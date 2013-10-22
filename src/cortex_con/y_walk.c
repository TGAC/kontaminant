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
#include <assert.h>
#include <pthread.h>

#include <global.h>
#include <flags.h>
#include <binary_kmer.h>
#include <element.h>
#include <open_hash/hash_table.h>
#include <dB_graph.h>
#include <path.h>
#include <perfect_path.h>
#include <logger.h>
#include <y_walk.h>

static pathStep *get_next_step(pathStep * current_step, pathStep * next_step,
                               pathStep * reverse_step, dBGraph * db_graph)
{
	pathStep *step = db_graph_get_next_step(current_step, next_step, reverse_step,
                                            db_graph);
	assert(step != NULL);
	
	if (step->node != NULL) {
		step->label = Undefined;
		
		
		Nucleotide n;
		db_node_has_precisely_one_edge(next_step->node, next_step->orientation, &n);
		
		if (db_node_edges_count(next_step->node, next_step->orientation) > 0) {
			next_step->label = n;
		}
	}
    
	return next_step;
}


/*
 * Returns false if the node has one and only one in the orientation given by the next_step.   
 * */
static boolean is_conflicted(pathStep * current_step, pathStep * next_step,
                             pathStep * reverse_step, dBGraph * db_graph)
{
    
	if (current_step->label == Undefined){//We are in the last node, therefore, it is the end of the walk. 
		return true;
	} else  {
		return false;
	}
}

static void path_step_do_nothing(pathStep * ps)
{
	return;
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

static boolean continue_searching(pathStep * current_step,
                                  pathStep * next_step,
                                  pathStep * reverse_step, Path * temp_path,
                                  dBGraph * db_graph)
{
	pathStep first;
	
	boolean cont = true;
	
	//check for cycles -- check if this can be done better
	if(temp_path->length > 0){
		path_get_step_at_index(0, &first, temp_path);
		cont = !path_step_equals_without_label(&first, next_step);
		//cont = !path_contains_step(next_step, temp_path);
        if (db_node_check_for_any_flag(next_step->node, next_step->orientation == forward? VISITED_FORWARD:VISITED_REVERSE)) {
            cont = false;
        }
	}
	//TODO: add the validations of VISITED_FORWARD and VISITED_REVERSE. 
    //TODO: Add the reasons to stop
	return !db_node_check_for_flag(current_step->node, Y_START) &&
    path_has_space(temp_path) &&
    !is_conflicted(current_step, next_step, reverse_step, db_graph) &&
    cont;
}


static boolean continue_traversing(pathStep * current_step,
                                   pathStep * next_step,
                                   pathStep * reverse_step, Path * temp_path,
                                   dBGraph * db_graph)
{
	//pathStep first;
	
	boolean cont = true;
    
	boolean next_visited = 
	(db_node_check_for_flag(next_step->node, VISITED_FORWARD) == true && next_step->orientation==forward )||
	(db_node_check_for_flag(next_step->node, VISITED_REVERSE) == true && next_step->orientation==reverse);
    
    if (next_visited && next_step->label != Undefined) {
        path_action_set_flag(temp_path, IS_CYCLE);
        path_add_stop_reason(LAST, PATH_FLAG_IS_CYCLE, temp_path);
        cont = false;
    }
    
    if (db_node_check_for_flag(current_step->node, Y_START)) {
        cont = false;
    }
    
    if (current_step->label == Undefined){
        path_add_stop_reason(LAST, PATH_FLAG_STOP_BLUNT_END, temp_path);
        cont = false;
    }
    
    if(!path_has_space(temp_path)){
        path_add_stop_reason(LAST, PATH_FLAG_LONGER_THAN_BUFFER, temp_path);
        cont = false;
    }
    return cont;
}



static void mark_from_end(dBNode * node, void * args) {
    WalkingFunctions  * wf  = ( WalkingFunctions  * ) args;
    pathStep first;
    first.node = node;
    first.orientation = undefined;
    first.path = NULL;
    Path *buffer = NULL;
    if (!db_node_check_for_any_flag (node, Y_START ) ) {	//We already marked this as an end-point
        int fwd_count = db_node_edges_count(node, forward);
        int rev_count = db_node_edges_count(node, reverse);
        if (fwd_count == 1 && rev_count > 1) {
            buffer = path_get_buffer_path();
            first.orientation = forward;
        } else if (rev_count == 1 && fwd_count > 1) {
            first.orientation = reverse;
        }
    }
    if (first.orientation != undefined) {
        buffer = path_get_buffer_path();
        path_reset(buffer);
        first.node = node;
        
        //just to get label of step
        db_node_has_precisely_one_edge(node, first.orientation, &first.label);
        
        db_graph_generic_walk(&first, buffer, wf, wf->db_graph);
        path_free_buffer_path(buffer);
       
    }
}


static void mark_double_y_callback(Path * p, void * args) {
    
    pathStep ps;
    
    long long * double_y_count  = (long long *) args;
    path_get_last_step(&ps, p);
    
    if (db_node_edges_count(ps.node, ps.orientation) > 1) {
        db_node_action_set_flag(path_last_node(p), Y_START);
        (*double_y_count)++;
    }
}
void mark_double_y(dBGraph * db_graph)
{
    
	WalkingFunctions wf;
	perfect_path_get_funtions(&wf);
	wf.continue_traversing = &continue_searching;
	wf.get_starting_step = &get_first_step_identity;
	wf.pre_step_action = &pre_step_action;
	wf.post_step_action = &post_step_action;
	wf.db_graph = db_graph;

    void * args[1];
    args[0] = &wf;
    
	long long double_y_count = 0;
	db_graph_add_path_callback_with_args(&wf, &mark_double_y_callback, (void * ) &double_y_count);
    
	log_and_screen_printf("Marking double Ys\n");
	hash_table_traverse_with_args(&mark_from_end, args, db_graph);
	log_and_screen_printf("%'qd  double Ys found\n", double_y_count);
    
}


/**
 * 
 * 
 * We have to think in a better way to make it reentrant if we want
 * to multithread this
 * 
 */
Path *y_walk_get_path(dBNode * node, Orientation orientation,
                      void (*node_action) (dBNode * node),
                      dBGraph * db_graph, boolean both_directions,Path * path)
{
    
	Path *buff1 = path_get_buffer_path();
	
	path_reset(buff1);
	path_reset(path);	//Make sure that the paths are clean
    
	pathStep first;
	boolean only_one_edge;
	first.node = node;
	first.orientation =  orientation;
	
	WalkingFunctions wf;
	perfect_path_get_funtions(&wf);
	
    db_graph_add_node_action(&wf, node_action);
	
	
	wf.continue_traversing = &continue_traversing;
	
	wf.get_next_step =  &get_next_step;
//	wf.step_action = &step_action;
	wf.pre_step_action = &pre_step_action;
	wf.post_step_action = &post_step_action;
	
    
	if(orientation != undefined){
		only_one_edge = db_node_has_precisely_one_edge(first.node, first.orientation, &first.label);
		wf.get_starting_step = &get_first_step_identity;
		
	}else {
        first.orientation = reverse;
        only_one_edge = db_node_has_precisely_one_edge(first.node, first.orientation, &first.label);
    }
    
    
	//we got the first side. We use path, since from here we will append
	//the second half of the path.
    
 //   db_graph_add_path_callback(<#WalkingFunctions *wf#>, <#void (*path_callback)(Path *)#>)
    if (both_directions) {
        db_graph_add_path_callback_with_args(&wf, &path_reverse, path);
    }else{
        path_r
    }
    
    
	void callback(Path * p) {
        // We have to copy the path, either in reverse or as it is, otherwise we lose it.
        if (both_directions == true) {
            path_reverse(p, path);
        } else {
            path_copy(path, p);
        }
        // printf("Ywalk first callback path: %s\n", p->seq);
	}
    
	//wf.output_callback = &callback;
    
	path_reset(buff1);
	if (db_node_edges_count(first.node, first.orientation) > 0 &&  !db_node_check_for_flag(first.node, Y_START)) {		
		db_graph_generic_walk(&first, buff1, &wf, db_graph);
	} else {
        if (both_directions == true) {
            first.label = Undefined;
            // RML: This wasn't reversing the orientation - shouldn't it?
            first.orientation = opposite_orientation(first.orientation);
            path_add_node(&first, path); //In this way, we can use the last node samelesly for the second round, without checking if it was able to walk in the first place
        }
	}
	
	if (both_directions == true){		
		//The second time, we just need to append
		void callback2(Path * p) {
			path_append(path, p);
		}
		
		wf.output_callback = &callback2;
		path_get_last_step(&first, path);//Starting over from the last step in the path. 
		db_node_has_precisely_one_edge(first.node, first.orientation, &first.label);
        
		if(db_node_edges_count(first.node, first.orientation) > 0 ){
			db_graph_generic_walk(&first, buff1, &wf, db_graph);
		}
        
	}
	
	path_free_buffer_path(buff1);
    
    return path;	
}

void y_walk_print_paths(char *filename, int max_length, int singleton_length, 
                        boolean with_coverages, boolean with_viz, dBGraph * db_graph)
{
    
	Path *path = path_get_buffer_path();	//We will try to use only this buffer path.
    
	path_reset(path);
    
	WalkingFunctions wf;
    
	FILE *fout;
	FILE *fout_cov;
	FILE *fout_viz;
	
	fout = fopen(filename, "w");
    
	if (with_coverages) {
		char filename_cov[strlen(filename) + 10];
		sprintf(filename_cov, "%s_cov", filename);
		fout_cov = fopen(filename_cov, "w");
	}
	
	if (with_viz) {
		char filename_cov[strlen(filename) + 10];
		sprintf(filename_cov, "%s.viz", filename);
		fout_viz = fopen(filename_cov, "w");
		path_graphviz_open_header(fout_viz);
	}
	
	int count_nodes = 0;
    
	//path->id=-1;
	long long count_kmers = 0;
	long long count_sing = 0;
    long long count_rep = 0;
    double graph_cov = db_graph_get_average_coverage(db_graph);
	PathCounts counts;
    path_counts_reset(&counts); 
    
	mark_double_y(db_graph);
	path_reset(path);
    
	perfect_path_get_funtions(&wf);
	wf.continue_traversing = &continue_traversing;
    //TODO: write a function that actually looks for where to start... 
	wf.get_starting_step = &get_first_step_identity;
    
	void print_supernode(dBNode * node) {
        
		count_kmers++;
		if (db_node_check_flag_visited(node) == false && db_node_check_flag_not_pruned(node)) {
            
			y_walk_get_path(node, undefined,
                            &db_node_action_set_flag_visited,
                            db_graph, true,path);
            
			if (path_is_singleton(singleton_length, path)) {
                count_sing++;
            }else if(path_is_repetitive(graph_cov, path)){
                count_rep++;
            }else {
				if (with_coverages) {
					path_to_coverage(path, fout_cov);
				}
				if(with_viz){
					path_graphviz_line(fout_viz, path);
				}
#ifdef SOLID
                if(db_graph->output_base_space){
					path_to_base_space_fasta(path, fout);
                }else{
                  	path_to_fasta(path, fout);   
                }
#else                
				path_to_fasta(path, fout);
#endif
				if (path->length == max_length) {
					log_and_screen_printf("contig length equals max length [%i] for node_%i\n", max_length, count_nodes);
				}
				path_counts_add(path, &counts);
				count_nodes++;
				path_increase_id(path);
			}
		}
	}
    
	hash_table_traverse(&print_supernode, db_graph);
	log_and_screen_printf("%'d nodes visited [%'qd singletons, %'qd repetitive]\n", count_nodes, count_sing);
	path_counts_print_and_log(&counts);
	
	fclose(fout);
	if (with_coverages) {
		fclose(fout_cov);
	}
	if (with_viz) {
		path_graphviz_close_header(fout_viz);
		fclose(fout_viz);
	}
    
    path_free_buffer_path(path);
    
	return;
}


