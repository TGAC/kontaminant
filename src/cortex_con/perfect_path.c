/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * Development team: 
 *       R. Ramirez-Gonzalez (Ricardo.Ramirez-Gonzalez@bbsrc.ac.uk)
 *       R. Leggett (richard@leggettnet.org.uk)
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */
 
#include <string.h>
#include <stdint.h>
#ifdef THREADS
#include <pthread.h>
#endif
#include <assert.h>


#include <global.h>
#include <flags.h>
#include <binary_kmer.h>
#include <element.h>
#include <open_hash/hash_table.h>
#include <dB_graph.h>
#include <path.h>
#include <perfect_path.h>
#include <logger.h>

static int limit = MAX_PATH_LENGTH;

static pathStep *get_next_step(pathStep * current_step, pathStep * next_step, pathStep * reverse_step, dBGraph * db_graph)
{
    pathStep *step = db_graph_get_next_step(current_step, next_step, reverse_step, db_graph);

    assert(step != NULL);
    
    step->label = Undefined;
    Nucleotide n;
   
    // TODO: If this works, then there's possibly a bug to fix for the multi-colour situation
    if (next_step->node == NULL) {
        next_step->label = Undefined;
        return next_step;
    }
 
    if (db_node_has_precisely_one_edge_all_colours(next_step->node, next_step->orientation, &n)) {
        next_step->label = n;
    }

    if (!db_node_has_precisely_one_edge_all_colours(next_step->node, opposite_orientation(next_step->orientation), &n)) {
        next_step->label = Undefined;
    }	

    return next_step;
}


static void path_step_do_nothing(pathStep * ps)
{
    return;
}

static boolean clean_path(Path * p, dBGraph * db_graph)
{
    return !path_is_empty(p);
}

static boolean continue_traversing(pathStep * current_step,
                                   pathStep * next_step,
                                   pathStep * reverse_step, Path * temp_path,
                                   dBGraph * db_graph)
{
	pathStep first;
	
	boolean cont;
    cont = current_step->label != Undefined;
    //cont = cont && db_node_has_precisely_one_edge_all_colours(current_step->node, current_step->orientation, &n);
    path_get_step_at_index(0, &first, temp_path);
    int n_fwd, n_rev;
    /* We don't do these checks for the first node - in case it's a Y node */
    if(temp_path->length > 1) {
        if (db_node_check_for_any_flag(next_step->node, next_step->orientation == forward? VISITED_FORWARD:VISITED_REVERSE)) {
            cont = false;
        }
        
        /* Check for a cycle - as this is a perfect path, we only need to check the first node. If we come
           back in at one of the other nodes, then it will result in two edges in one orientation */
        
        if (path_step_equals_without_label(&first, current_step)) {
            path_add_stop_reason(LAST, PATH_FLAG_IS_CYCLE, temp_path);
            cont = false;
        }
        /* Now check for more than one edge in either direction */
        n_fwd = db_node_edges_count_all_colours(current_step->node, current_step->orientation);
        n_rev = db_node_edges_count_all_colours(current_step->node, opposite_orientation(current_step->orientation));
        
        if (n_fwd == 0) {
            path_add_stop_reason(LAST, PATH_FLAG_STOP_BLUNT_END, temp_path);
            cont = false;
        }
        if (n_fwd > 1) {
            path_add_stop_reason(LAST, PATH_FLAG_DIVERGING_PATHS, temp_path);
            cont = false;
        }
        if (n_rev > 1) {
            path_add_stop_reason(LAST, PATH_FLAG_CONVERGING_PATHS, temp_path);
            cont = false;
        }
        if (!path_has_space(temp_path)) {
            path_add_stop_reason(LAST, PATH_FLAG_LONGER_THAN_BUFFER, temp_path);
            cont = false;
        }      
    }
                                                  
    
    
	return cont;
}


pathStep *get_first_step_identity(pathStep * first_step, dBGraph * db_graph){
	return first_step;
}


static void store_last(Path * p, void * arg){
    pathStep * last = (pathStep *) arg;
    path_get_last_step_reverse(last, p);
}

static pathStep *get_first_step(pathStep * first_step, dBGraph * db_graph)
{
	//printf("|");
	WalkingFunctions wf;
	pathStep tmp_step;
	Nucleotide n;

	
    path_step_assign(&tmp_step, first_step);
	if(tmp_step.orientation == undefined){
		tmp_step.orientation = reverse; //This will allow, when reversing the walk, to walk in the forward direction, which "naturally" appends the bases to the end of the path. 
	}else{
		return first_step;
	}

	/* Check if forward node has one edge, if not that reverse node does.
       While we're at it, assign label to n. */
	if (db_node_has_precisely_one_edge_all_colours(tmp_step.node, tmp_step.orientation, &n)) {	
	} else if (db_node_has_precisely_one_edge_all_colours(tmp_step.node, opposite_orientation(tmp_step.orientation), &n)) {
        tmp_step.orientation = opposite_orientation(tmp_step.orientation);
    } else {
 
        return NULL;
    }    
	tmp_step.label = n;

    /* If it's a blunt node in the opposite orientation, then we can't get a better starting step */
    if (db_node_is_blunt_end_all_colours(tmp_step.node, opposite_orientation(tmp_step.orientation))) {
		path_step_assign(first_step, &tmp_step);
		return first_step;
	}
    
    //printf("#");
    /* Otherwise, do a perfect path walk to try and get a better start step */    
	//Path * temp_path = path_new(limit, db_graph->kmer_size);
	Path *temp_path = path_get_buffer_path();
	temp_path->id = -2;
	
	
	db_node_has_precisely_one_edge_all_colours(tmp_step.node, tmp_step.orientation, &tmp_step.label);
	wf.get_starting_step = &get_first_step_identity;
	wf.continue_backwards = &clean_path;
	wf.post_step_action = &path_step_do_nothing;
	wf.pre_step_action = &path_step_do_nothing;
	wf.step_action = &path_step_do_nothing;
	wf.continue_traversing = &continue_traversing;
	wf.get_next_step = &get_next_step;
    wf.output_callback = &path_do_nothing;
    
    wf.node_callbacks.used = 0;
	wf.step_actions.used = 0;
    wf.path_callbacks.used = 0 ;
    
	
    db_graph_add_path_callback_with_args(&wf, &store_last, (void *) first_step);
    
	
	if (DEBUG) {
		printf("[get_first_step]Orientation: %s\n", tmp_step.orientation == forward ? "forward" : "reverse");
		printf("[get_first_step] (%d):", path_get_edges_count(temp_path));
		path_step_print(first_step, db_graph->kmer_size, stdout);
		printf("\n");
    }
    
	db_graph_generic_walk(&tmp_step, temp_path, &wf, db_graph);
	
    path_free_buffer_path(temp_path);
	return first_step;
}

void perfect_path_base_callback(Path * p)
{
	pathStep p1, p2;
	path_get_last_step(&p1, p);
	path_get_step_at_index(0, &p2, p);
	
	//Remove the repeated node in a cycle
	if (path_is_cycle(p)) {
		
		if (path_step_equals
		    (&p1,
		     &p2)) {
				p1.label = Undefined;
				path_remove_last(p);
				path_add_node(&p1, p);
			}
	}
	
	int n_fwd = db_node_edges_count_all_colours(p2.node, p2.orientation);
    int n_rev = db_node_edges_count_all_colours(p2.node, opposite_orientation(p2.orientation));
    
    int n_fwd_f = db_node_edges_count_all_colours(p1.node, p1.orientation);
    int n_rev_f = db_node_edges_count_all_colours(p1.node, opposite_orientation(p1.orientation));
    
        
    if (n_fwd == 0) {
        path_add_stop_reason(FIRST, PATH_FLAG_STOP_BLUNT_END, p);
    }
    if (n_fwd > 1) {
        path_add_stop_reason(FIRST, PATH_FLAG_DIVERGING_PATHS, p);
    }
    
    if (n_rev > 1) {
        path_add_stop_reason(FIRST, PATH_FLAG_CONVERGING_PATHS, p);
    }
        
	if(n_rev  > 1 && n_fwd_f > 1 ){
		path_add_stop_reason(FIRST, PATH_FLAG_IS_DOUBLE_Y,p);
	}
    
    if(n_rev_f > 1){
        path_add_stop_reason(LAST, PATH_FLAG_CONVERGING_PATHS, p);
    }
    
    if(n_fwd_f > 1){
        path_add_stop_reason(LAST, PATH_FLAG_DIVERGING_PATHS, p);
    }
	
	
}

static void pre_step_action(pathStep * ps){
	
	if(ps->orientation == forward){
        db_node_action_set_flag(ps->node, VISITED_FORWARD);
	}else{
        db_node_action_set_flag(ps->node, VISITED_REVERSE);
	}
	
}

static void post_step_action(pathStep * ps){
	if(ps->orientation == forward){
        db_node_action_unset_flag(ps->node, VISITED_FORWARD);
		
	}else{
        db_node_action_unset_flag(ps->node, VISITED_REVERSE);
	}
}


WalkingFunctions * perfect_path_get_funtions(WalkingFunctions *
											 walking_functions)
{
	//walking_functions->find_first_node = &always_true;
	walking_functions->get_starting_step = &get_first_step;
	walking_functions->continue_backwards = &clean_path;
	walking_functions->post_step_action = &post_step_action;
	walking_functions->pre_step_action = &pre_step_action;
	
	walking_functions->continue_traversing = &continue_traversing;
	walking_functions->get_next_step = &get_next_step;
	walking_functions->output_callback = &perfect_path_base_callback;
	walking_functions->step_action = &path_step_do_nothing;
    
    //Clearing the array of functions to concatenate. 
    walking_functions->step_actions.used = 0;
    walking_functions->path_callbacks.used = 0;
    walking_functions->node_callbacks.used = 0;
	
	return walking_functions;
}





int perfect_path_get_path_with_callback(dBNode * node, Orientation orientation, void (*node_action) (dBNode * node), void (*path_action) (Path * path),dBGraph * db_graph)
{
	
	WalkingFunctions wf;
	perfect_path_get_funtions(&wf);
	
	pathStep first;
	first.node = node;
	//first.orientation = orientation == undefined? reverse: orientation;
	first.orientation = orientation;
	first.label = Undefined;
	
	if (orientation != undefined) {
		db_node_has_precisely_one_edge_all_colours(node, orientation, &first.label);
		wf.get_starting_step = &get_first_step_identity;	
	}
	
	//wf.find_first_node = &find_first;
	
	
    db_graph_add_node_action(&wf, node_action);
    db_graph_add_path_callback(&wf, path_action);
    
	
	Path *path = path_get_buffer_path();
	
	int ret = db_graph_generic_walk(&first, path, &wf, db_graph);
	path_free_buffer_path(path);
	
	return ret;
}


int perfect_path_get_path_with_callback_with_args(dBNode * node, Orientation orientation, void (*node_action) (dBNode * , void *), void * node_args,  void (*path_action) (Path * , void *), void * path_args,dBGraph * db_graph)
{
	
	WalkingFunctions wf;
	perfect_path_get_funtions(&wf);
	
	pathStep first;
	first.node = node;
	//first.orientation = orientation == undefined? reverse: orientation;
	first.orientation = orientation;
	first.label = Undefined;
	
	if (orientation != undefined) {
		db_node_has_precisely_one_edge_all_colours(node, orientation, &first.label);
		wf.get_starting_step = &get_first_step_identity;	
	}
	
	//wf.find_first_node = &find_first;
	
	
    db_graph_add_node_action_with_args (&wf, node_action, node_args);
    db_graph_add_path_callback_with_args (&wf, path_action, path_args);
    
	
	Path *path = path_get_buffer_path();
	
	int ret = db_graph_generic_walk(&first, path, &wf, db_graph);
	path_free_buffer_path(path);
	
	return ret;
}

int perfect_path_get_path_from_step_with_callback_with_args(pathStep  first,
												  void (*node_action) (dBNode * node, void *),void * node_args,
												  void (*path_action) (Path * path, void *), void * path_args,
												  dBGraph * db_graph)
{
	
	WalkingFunctions wf;
	perfect_path_get_funtions(&wf);
	
	
	
	
	wf.get_starting_step = &get_first_step_identity;
	
    db_graph_add_node_action_with_args (&wf, node_action, node_args);
    db_graph_add_path_callback_with_args (&wf, path_action, path_args);
	
	Path *path = path_get_buffer_path();
	
	int ret = db_graph_generic_walk(&first, path, &wf, db_graph);
	path_free_buffer_path(path);
	
	return ret;
}


static void pp_copy_path(Path * p, void * args) {
    Path * path = (Path *) args;    
    path_copy(path, p);
}
int perfect_path_get_path(dBNode * node, Orientation orientation,
						  void (*node_action) (dBNode * node),
						  dBGraph * db_graph, Path * path)
{
	
	
	
	perfect_path_get_path_with_callback_with_args(node, orientation, (void *) node_action, NULL,  &pp_copy_path,(void *)  path, db_graph);
	return path_get_edges_count(path);
	
}


typedef struct{
    long long count_kmers;
    dBGraph * db_graph;
    Path * path;
    int singleton_length;
    long long count_sing;
    long long count_rep;
    double graph_cov;
    FILE * fout_cov;
    FILE * fout;
    PathCounts  counts;
    long long count_nodes;
} perfect_path_print_supernodes_args;

static 	void print_supernode(dBNode * node, void * arg) {
    
    perfect_path_print_supernodes_args * args = (perfect_path_print_supernodes_args * ) arg;
    
    args->count_kmers++;
    if (db_node_check_flag_visited(node) == false && db_node_check_flag_not_pruned(node)) {
        
        perfect_path_get_path(node, undefined,
                              &db_node_action_set_flag_visited,
                              args->db_graph, args->path);
        
        if (path_is_singleton(args->singleton_length, args->path)) {
            args->count_sing++;
        }else if(path_is_repetitive(args->graph_cov, args->path)){
            args->count_rep++;
        }else{
            if (args->fout_cov != NULL) {
                path_to_coverage(args->path, args->fout_cov);
            }
#ifdef SOLID
            if(args->db_graph->output_base_space){
                path_to_base_space_fasta(args->path, args->fout);
            }else{
                path_to_fasta(args->path, args->fout);
            }
#else                
            path_to_fasta(args->path, args->fout);
#endif
            
            
            path_counts_add(args->path, &args->counts);
            args->count_nodes++;
            path_increase_id(args->path);
        }
        
        
    }
}

void perfect_path_print_paths(char *filename, int max_length, int singleton_length,
							  boolean with_coverages, dBGraph * db_graph)
{
	FILE * fout = NULL;
	FILE * fout_cov = NULL;
	int i;
	fout = fopen(filename, "w");
	
	if (with_coverages) {
		char filename_cov[strlen(filename) + 10];
		sprintf(filename_cov, "%s_cov", filename);
		fout_cov = fopen(filename_cov, "w");
	}
	limit = max_length;
	
	//Path *path = path_new(max_length, db_graph->kmer_size);
	//path->id=-1;
    
    perfect_path_print_supernodes_args ** args = calloc(db_graph->number_of_threads, sizeof(perfect_path_print_supernodes_args * ));
    
    for (i = 0; i < db_graph->number_of_threads; i++) {
        args[i] = calloc(1, sizeof(perfect_path_print_supernodes_args)) ;
        args[i]->db_graph = db_graph;
        args[i]->path = path_new(max_length, db_graph->kmer_size);
        args[i]->fout = fout;
        args[i]->fout_cov = fout_cov;//TODO: Make the printing function "thread safe"
    }
  
	//buffers = path_array_new(2);
	
	
	double graph_cov = db_graph_get_average_coverage(db_graph);
    log_and_screen_printf("Average coverage: %5.2f \n", graph_cov);
    
    
	
	hash_table_traverse_with_args(&print_supernode, (void ** ) args ,db_graph);
	
	
	log_and_screen_printf("%'d nodes visited [%'qd singletons, %'qd repetitive]\n", args[0]->count_nodes, args[0]->count_sing, args[0]->count_rep);//TODO: At some point we can make this multithreading
	path_counts_print_and_log(&args[0]->counts);
    
    
	for (i = 0; i < db_graph->number_of_threads; i++) {
        free(args[i]);
        path_destroy(args[i]->path);
    }
    free(args);
	
	
	fclose(fout);
	if (with_coverages) {
		fclose(fout_cov);
	}
}

//Gets all the paths from a given node. 
int perfect_path_get_all_paths_from(dBNode * node, Orientation orientation,
									 PathArray * pa, int limit, dBGraph * db_graph)
{
	
	int found = 0;
	pathStep ps;
	ps.orientation = orientation;
	ps.node = node;
	ps.flags = 0;
    ps.path = NULL;
    Nucleotide n;
    for (n = 0; n < 4; n++) {
		ps.label = n;
		path_reset(pa->paths[n]);
        path_set_limit(limit,pa->paths[n] );		
		if (db_node_edge_exist_any_colour(node, n, orientation)) {
			perfect_path_get_path_from_step_with_callback_with_args(ps,
														  NULL, NULL,
														  &pp_copy_path,pa->paths[n], 
														  db_graph);
			
		}

    }
	
	return found;
}
