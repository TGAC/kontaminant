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
#include <branches.h>
#include <logger.h>


static boolean continue_traversing(pathStep * current_step,
								   pathStep * next_step,
								   pathStep * reverse_step, Path * temp_path,
								   dBGraph * db_graph)
{
	pathStep first;
	
	boolean cont = true;
    cont = current_step->label != Undefined;
    //cont = cont && db_node_has_precisely_one_edge(current_step->node, current_step->orientation, &n);
    path_get_step_at_index(0, &first, temp_path);
    int n_fwd, n_rev;
    /* We don't do these checks for the first node - in case it's a Y node */
    if(temp_path->length > 1) {
        if (db_node_check_for_any_flag(next_step->node, next_step->orientation == forward? VISITED_FORWARD:VISITED_REVERSE)) {
            cont = false;
        }
        
        if (path_step_equals_without_label(&first, current_step) || path_has_in_step(next_step, temp_path)) {
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
       
    }
    
    if(path_get_length(temp_path) >= path_get_limit(temp_path)){
        cont = false;
        path_add_stop_reason(LAST, PATH_FLAG_LONGER_THAN_BUFFER, temp_path);
    }
    
    if(temp_path->in_nodes_count > db_graph->max_double_y_complexity && temp_path->out_nodes_count > db_graph->max_double_y_complexity){
        cont = false;
        path_add_stop_reason(LAST, PATH_TOO_COMPEX, temp_path);
    }
    
    
	return cont;
}

static pathStep *get_next_step(pathStep * current_step, pathStep * next_step,
							   pathStep * reverse_step, dBGraph * db_graph)
{
    
	pathStep *step = db_graph_get_next_step(current_step, next_step, reverse_step, db_graph);
    
    assert(step != NULL);
    if (step->node != NULL) {
		step->label = Undefined;
		Nucleotide n;
		if (db_node_has_precisely_one_edge_all_colours
		    (next_step->node, next_step->orientation, &n)) {
			 next_step->label = n;
		}else{
            next_step->label = path_step_get_unvisited_edge_all_colours(next_step);
        }
        
	}
	
	return next_step;
}

pathStep *get_first_step_label(pathStep * first_step, dBGraph * db_graph){
    assert(first_step->node != NULL);
    if (first_step->label == Undefined) {
         first_step->label = path_step_get_unvisited_edge_all_colours(first_step);
    }
	return first_step;
}

static boolean clean_path(Path * p, dBGraph * db_graph){
	if(path_is_empty(p)){
        return false;
    }
    boolean clean = true;
    pathStep last;
    
    path_get_last_step(&last, p);
    if(path_step_has_unvisited_edge_all_colours(&last)){
        Nucleotide n = path_step_get_unvisited_edge_all_colours(&last);
        path_modfy_last_label(n, p);
        clean = false;
    }
    
	return clean;
}

WalkingFunctions * branches_get_funtions(WalkingFunctions * walking_functions){
    perfect_path_get_funtions(walking_functions);
    walking_functions->continue_backwards = &clean_path;
    walking_functions->get_next_step = &get_next_step;
    walking_functions->continue_traversing = &continue_traversing;
    walking_functions->get_starting_step = &get_first_step_label;
    return walking_functions;
}


/**
 * Returns all the paths up to certain length. All the paths come from the pathBuffer, you are responsible
 * of freeing them from the buffer with path_array_free_from_buffer after using the buffer. 
 */
PathArray * branches_get_all_paths_from(pathStep * first, int max_length, dBGraph * db_graph){
//    first->label = Cytosine;
    first->flags = 0;
    first->orientation = reverse;
    first->label = Undefined;
    PathArray * pa = path_array_new(4);// where is freeing the path array?
    Path * buff = path_get_buffer_path();
    path_set_limit(max_length, buff);
    
    WalkingFunctions wf;
    branches_get_funtions(&wf);
    
    
    void copy_path(Path * p){
        Path * tmp = path_get_buffer_path();
        path_copy(tmp, p);
        //path_to_fasta(p, stdout);
        path_array_add_path(tmp, pa);
    }
    wf.output_callback = &copy_path;
    
    
    db_graph_generic_walk(first, buff, &wf, db_graph);
    
    path_free_buffer_path(buff);
    return pa;
}

void branches_get_all_paths_from_with_callback(pathStep * first, int max_length, void (*path_action) (Path * path), dBGraph * db_graph){
    //    first->label = Cytosine;
    first->flags = 0;
    first->orientation = reverse;
    first->label = Undefined;
    Path * buff = path_get_buffer_path();
    path_set_limit(max_length, buff);
    
    WalkingFunctions wf;
    branches_get_funtions(&wf);
    
    
    wf.output_callback = path_action;
    
    
    db_graph_generic_walk(first, buff, &wf, db_graph);
    
    path_free_buffer_path(buff);
    return ;
}



/*
 * Overwrites the content of path. Returns the path if found, NULL otherwise. 
 */
Path * branches_get_path_between(pathStep * first, pathStep * last, int max_length, int max_depth, Path * path,  dBGraph * db_graph){
    //    first->label = Cytosine;
    first->flags = 0;
    first->orientation = reverse;
    first->label = Undefined;
    Path * tmp = NULL;
    Path * buff = path_get_buffer_path();
    path_set_limit(max_length, buff);
    boolean found = false;
    
    WalkingFunctions wf;
    branches_get_funtions(&wf);
    
    boolean ( * cont_old)(pathStep * current_step, pathStep * next_step, pathStep * reverse_step, Path * temp_path, dBGraph * db_graph) =  wf.continue_traversing;
    
    boolean continue_traversing2(pathStep * current_step, pathStep * next_step, pathStep * reverse_step, Path * temp_path, dBGraph * db_graph2){
        found = binary_kmer_comparison_operator(current_step->node->kmer, last->node->kmer);
        if(found){
            return false;
        }
        if (temp_path->out_nodes_count > max_depth) {
            return false;
        }
        return cont_old(current_step, next_step,reverse_step,  temp_path, db_graph2);
    }
    
    boolean (*continue_backwards_old)(Path * path, dBGraph * db_graph) = wf.continue_backwards;
    boolean continue_backwards2(Path * path2, dBGraph * db_graph2){
        if(found){
            if(path_get_length(path2) > 0){
                return true;
            }else{
                return false;
            }
        }else{
            return   continue_backwards_old(path2, db_graph2) ; //We stop searching once we have a match. 
        }
    }
    wf.continue_backwards = &continue_backwards2;
    wf.continue_traversing = &continue_traversing2;
    void copy_path(Path * p){
        if(found){
            tmp = path;
            path_copy(tmp, p);
        }
    }
    wf.output_callback = &copy_path;
    db_graph_generic_walk(first, buff, &wf, db_graph);
    path_free_buffer_path(buff);
    
    return tmp;
}

/*
 * Overwrites the content of path. Returns the path if found, NULL otherwise. 
 */
Path * branches_get_path_between_paths(Path * first_path, Path * last_path, int max_length, int max_depth, Path * path,  dBGraph * db_graph){
    //    first->label = Cytosine;
    pathStep first;
    pathStep first_last;

    pathStep second_first;
    pathStep second_last;
 //   pathStep current_step_in_read;
    
    path_get_step_at_index(0, &second_first, last_path);
    path_get_last_step(&second_last, last_path);
    path_get_step_at_index(0, &first, first_path);
    path_get_last_step(&first_last, first_path);
    
    Path * tmp = NULL;
    Path * buff = path_get_buffer_path();
    path_set_limit(max_length, buff);
    boolean found = false;
    boolean found_second_first = false;
    boolean found_second_last = false;
    boolean found_first_last = false;
//    int current_index;
//    int first_length = path_get_length(first_path);
    
    
    WalkingFunctions wf;
    branches_get_funtions(&wf);
    
    boolean ( * cont_old)(pathStep * current_step, pathStep * next_step, pathStep * reverse_step, Path * temp_path, dBGraph * db_graph) =  wf.continue_traversing;
    
    boolean continue_traversing2(pathStep * current_step, pathStep * next_step, pathStep * reverse_step, Path * temp_path, dBGraph * db_graph2){
        
       //int curr_len = path_get_length(temp_path);
        
        /*if(curr_len <= first_length){//TODO: make this more flexible, think on 454 homopolymers. 
            path_get_step_at_index(curr_len-1, &current_step_in_read, temp_path);
            if (current_step_in_read.node != current_step->node) {
                return false;
            }
        }*/
        
        if(!found_first_last){
            found_first_last = first_last.node == current_step->node;
        }
        
        if(!found_second_first){
            found_second_first = second_first.node == current_step->node;
        }
        if(!found_second_last){
            found_second_last = second_last.node == current_step->node;
        }
        
        
        
        //found = binary_kmer_comparison_operator(current_step->node->kmer, last.node->kmer);
        found = found_second_last && found_second_first;
        if(found){
            return false;
        }
        if (temp_path->out_nodes_count > max_depth) {
            return false;
        }
        return cont_old(current_step, next_step,reverse_step,  temp_path, db_graph2);
    }
    
    boolean (*continue_backwards_old)(Path * path, dBGraph * db_graph) = wf.continue_backwards;
    boolean continue_backwards2(Path * path2, dBGraph * db_graph2){
        if(found){
            if(path_get_length(path2) > 0){
                return true;
            }else{
                return false;
            }
        }else{
            return   continue_backwards_old(path2, db_graph2) ; //We stop searching once we have a match. 
        }
    }
    wf.continue_backwards = &continue_backwards2;
    wf.continue_traversing = &continue_traversing2;
    void copy_path(Path * p){
        if(found){
            tmp = path;
            path_copy(tmp, p);
        }
    }
    wf.output_callback = &copy_path;
    db_graph_generic_walk(&first, buff, &wf, db_graph);
    path_free_buffer_path(buff);
    
    return tmp;
}


