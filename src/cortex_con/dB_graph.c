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

#include <structs.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <assert.h>
#ifdef THREADS
#include <pthread.h>
#endif
#include <structs.h>

#include <global.h>
#include <flags.h>
#include <binary_kmer.h>
#include <element.h>
#include <open_hash/hash_table.h>
#include <dB_graph.h>
#include <cleaning.h>
#include <file_reader.h>
#include <path.h>
#include <perfect_path.h>
#include <logger.h>

#ifdef ENABLE_READ_PAIR
#include <binary_tree.h>
#include <read_pair.h>
#endif

#ifdef ENABLE_MARK_PAIR
#include <binary_tree.h>
#include <mark_pair.h>
#endif

/**
 * Returns the next node, and the reverse to come back.
 * This function DO NOT  tells which is the next node
 * label to traverse, since it depends on the algorithm
 * to decide which one to use. It, however, gives the 
 * reverse lable to be able to traverse back. 
 *  
 */
#ifndef SOLID2
pathStep *db_graph_get_next_step(pathStep * current_step, pathStep * next_step, pathStep * rev_step, dBGraph * db_graph)
{
    
	assert(current_step != NULL);
	assert(next_step != NULL);
	assert(rev_step != NULL);
	
	assert(current_step->node != NULL);
	assert(current_step->label != Undefined);
	
	next_step->node = NULL;
	rev_step->node = current_step->node;
    next_step->flags = 0;
    rev_step->flags = current_step->flags & PATH_STEP_MASK_VISITED;
    
	BinaryKmer local_copy_of_kmer;
	binary_kmer_assignment_operator(local_copy_of_kmer,	current_step->node->kmer);
    
	BinaryKmer tmp_kmer;
	//dBNode * next_node = NULL;
    
	// after the following line tmp_kmer and rev_kmer are pointing to the same B Kmer
	BinaryKmer *rev_kmer =
    binary_kmer_reverse_complement(&local_copy_of_kmer, db_graph->kmer_size, &tmp_kmer);
    
	if (current_step->orientation == reverse) {
		rev_step->label = binary_kmer_get_last_nucleotide(&local_copy_of_kmer);
		binary_kmer_assignment_operator(local_copy_of_kmer, *rev_kmer);
	} else {//TODO: This could be avoided by reversing just the first nucleotide, not requiring to reverse always. 
		rev_step->label = binary_kmer_get_last_nucleotide(rev_kmer);
	}
    
	binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&local_copy_of_kmer, current_step->label, db_graph->kmer_size);
    
	//get node from table
	next_step->node = hash_table_find(element_get_key(&local_copy_of_kmer, db_graph->kmer_size, &tmp_kmer), db_graph);
	rev_step->node = next_step->node;
    
	if (next_step->node != NULL) {
		next_step->orientation = db_node_get_orientation(&local_copy_of_kmer, next_step->node, db_graph->kmer_size);
		rev_step->orientation = opposite_orientation(next_step->orientation);
	}    
    //#ifdef __DEBUG
	else {
        //		if (DEBUG) {
        //
        char tmpseq[db_graph->kmer_size];
        printf("[db_graph_get_next_step] Cannot find %s so get a NULL node\n",
               binary_kmer_to_seq(&tmp_kmer, db_graph->kmer_size, tmpseq));
        //Commented by ricardo, to reduce the log as for the traversing
        //      algorithm relays on having this as null
        //		}
	}
    //#endif
    
	next_step->label = Undefined;
	return next_step;
}

pathStep *db_graph_get_next_step_with_reverse(pathStep * current_step, pathStep * next_step, pathStep * rev_step, dBGraph * db_graph)
{
    
	assert(current_step != NULL);
	assert(next_step != NULL);
	assert(rev_step != NULL);
	
	assert(current_step->node != NULL);
	assert(current_step->label != Undefined);
	
	next_step->node = NULL;
	rev_step->node = current_step->node;
    next_step->flags = 0;
    rev_step->flags = current_step->flags & PATH_STEP_MASK_VISITED;
    
	BinaryKmer local_copy_of_kmer;
	binary_kmer_assignment_operator(local_copy_of_kmer,	current_step->node->kmer);
    
	BinaryKmer tmp_kmer;
	//dBNode * next_node = NULL;
    
	// after the following line tmp_kmer and rev_kmer are pointing to the same B Kmer
	BinaryKmer *rev_kmer =
    binary_kmer_reverse_complement(&local_copy_of_kmer, db_graph->kmer_size, &tmp_kmer);
    
	if (current_step->orientation == reverse) {
		rev_step->label = binary_kmer_get_last_nucleotide(&local_copy_of_kmer);
		binary_kmer_assignment_operator(local_copy_of_kmer, *rev_kmer);
	} else {//TODO: This could be avoided by reversing just the first nucleotide, not requiring to reverse always. 
		rev_step->label = binary_kmer_get_last_nucleotide(rev_kmer);
	}
    
	binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&local_copy_of_kmer, current_step->label, db_graph->kmer_size);
    
	//get node from table
	next_step->node = hash_table_find(element_get_key(&local_copy_of_kmer, db_graph->kmer_size, &tmp_kmer), db_graph);
	rev_step->node = next_step->node;
    
	if (next_step->node != NULL) {
		next_step->orientation = db_node_get_orientation(&local_copy_of_kmer, next_step->node, db_graph->kmer_size);
		rev_step->orientation = opposite_orientation(next_step->orientation);
	}    
    //#ifdef __DEBUG
	else {
        //		if (DEBUG) {
        //
        char tmpseq[db_graph->kmer_size];
        printf("[db_graph_get_next_step] Cannot find %s so get a NULL node\n",
               binary_kmer_to_seq(&tmp_kmer, db_graph->kmer_size, tmpseq));
        //Commented by ricardo, to reduce the log as for the traversing
        //      algorithm relays on having this as null
        //		}
	}
    //#endif
    
	next_step->label = Undefined;
	return next_step;
}


#else
pathStep *db_graph_get_next_step(pathStep * current_step, pathStep * next_step, pathStep * rev_step, dBGraph * db_graph)
{
    
	assert(current_step != NULL);
	assert(next_step != NULL);
	assert(rev_step != NULL);
	
	assert(current_step->node != NULL);
	
	assert(current_step->label != Undefined);
	
	next_step->node = NULL;
	rev_step->node = current_step->node;
    
	BinaryKmer local_copy_of_kmer;
	binary_kmer_assignment_operator(local_copy_of_kmer,	current_step->node->kmer);
    
	BinaryKmer tmp_kmer;
	//dBNode * next_node = NULL;
    
	// after the following line tmp_kmer and rev_kmer are pointing to the same B Kmer
	BinaryKmer *rev_kmer =
    binary_kmer_reverse_complement(&local_copy_of_kmer, db_graph->kmer_size, &tmp_kmer);//On solid, this is only asigning. 
    
	if (current_step->orientation == reverse) {
		rev_step->label = binary_kmer_get_first_nucleotide(&local_copy_of_kmer, db_graph->kmer_size);
        
        binary_kmer_right_shift_and_insert_new_base_at_left_end(&local_copy_of_kmer, current_step->label, db_graph->kmer_size);
        //  binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&local_copy_of_kmer, current_step->label, db_graph->kmer_size);
	} else {
		rev_step->label = binary_kmer_get_last_nucleotide(rev_kmer);
        binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&local_copy_of_kmer, current_step->label, db_graph->kmer_size);
	}
    
	
    
	//get node from table
	next_step->node = hash_table_find(element_get_key(&local_copy_of_kmer, db_graph->kmer_size, &tmp_kmer), db_graph);
	rev_step->node = next_step->node;
    
	if (next_step->node != NULL) {
		next_step->orientation = current_step->orientation;
		rev_step->orientation = opposite_orientation(next_step->orientation);
	}    
    //#ifdef __DEBUG
	else {
		if (DEBUG) {
            
			char tmpseq[db_graph->kmer_size];
			printf("[db_graph_get_next_step] Cannot find %s so get a NULL node\n",
                   binary_kmer_to_seq(&tmp_kmer, db_graph->kmer_size, tmpseq));
			//Commented by ricardo, to reduce the log as for the traversing
			//      algorithm relays on having this as null
		}
	}
    //#endif
    
	next_step->label = Undefined;
	return next_step;
}



#endif

/**
 * This function resets the flags of the whole graph. A future improvement
 * should be to do this in a traditional look, using a single binary operation. 
 * That will reduce the number of jumps, and should be paralellizable and vectorizable. 
 */
void db_graph_reset_flags(dBGraph * db_graph)
{
	if (DEBUG) {
		printf("[db_graph_identify_branches] Cleaning flags...\n");
	}
	printf("Cleaning flags...\n");
#ifdef THREADS
	hash_table_threaded_traverse(&db_node_action_clear_flags, db_graph);
#else
	hash_table_traverse(&db_node_action_clear_flags, db_graph);
#endif
	printf("\n");
	fflush(stdout);
}

int db_graph_get_perfect_path_with_first_edge(pathStep * first_step, void (*node_action) (dBNode* node), Path * path, dBGraph * db_graph)
{
    
	dBNode *node = first_step->node;
    
	Nucleotide fst_nucleotide = first_step->label;
	//printf("First node\n");
	//printNode(first_step->node, db_graph->kmer_size);
	dBNode *current_node = node;
	Nucleotide nucleotide, nucleotide2;
	int length = 0;
	char tmp_seq[db_graph->kmer_size + 1];
	tmp_seq[db_graph->kmer_size] = '\0';
    
	//sanity check
	if (node == NULL) {
		printf
        ("[db_graph_get_perfect_path_with_first_edge] db_graph_get_perfect_path: can't pass a null node\n");
		exit(1);
	}
    
	path_reset(path);
    
	if (DEBUG) {
		printf
        ("[db_graph_get_perfect_path_with_first_edge] Node %i in path: %s\n",
         path->length,
         binary_kmer_to_seq(element_get_kmer(current_node),
                            db_graph->kmer_size, tmp_seq));
	}
	//first edge defined
	nucleotide = fst_nucleotide;
	boolean added = true;
    
	pathStep current_step, next_step, rev_step;
	path_step_assign(&next_step, first_step);
	do {
		path_step_assign(&current_step, &next_step);
		/*current_step.node = next_step.node;
         current_step.orientation = next_step.orientation;
         current_step.label = next_step.label; */
		//printNode(current_step.node, db_graph->kmer_size);
		added = path_add_node(&current_step, path);
		if (added) {
            
			node_action(current_step.node);
			//                      db_graph_get_next_node(pathStep * current_step,
			//                                      pathStep * next_step, pathStep * rev_step, dBGraph * db_graph);
			db_graph_get_next_step(&current_step, &next_step,
                                   &rev_step, db_graph);
            
			//sanity check
			if (next_step.node == NULL) {
				fprintf(stderr,
                        "[db_graph_get_perfect_path_with_first_edge] dB_graph: didnt find node in hash table: %s %c %s\n",
                        binary_kmer_to_seq(element_get_kmer
                                           (current_step.node),
                                           db_graph->kmer_size,
                                           tmp_seq),
                        binary_nucleotide_to_char
                        (current_step.label),
                        current_step.orientation ==
                        forward ? "forward" : "reverse");
				exit(1);
			}
            
		}
	}
	while (added && !path_is_cycle(path) &&	//loop
	       db_node_has_precisely_one_edge(next_step.node, opposite_orientation(next_step.orientation), &nucleotide2) &&	//multiple entries
	       db_node_has_precisely_one_edge(next_step.node, next_step.orientation, &next_step.label));	//has one next edge only
	if (current_step.node != NULL) {
		added = path_add_node(&next_step, path);
		if (added) {
			node_action(next_step.node);
		}
	}
	if (DEBUG) {
		if (next_step.node == NULL) {
			printf("\n[db_graph_get_perfect_path_with_first_edge] Next node is null! \n");
		} else {
			printf("\n[db_graph_get_perfect_path_with_first_edge] Last node in path: %s %i length: %i\n",
                   binary_kmer_to_seq(element_get_kmer(next_step.node), db_graph->kmer_size, tmp_seq),
                   db_node_get_edges(next_step.node),	//next_node->edges,
			       path->length);
		}
	}
    
	return length;
    
}

// An all colours version of db_graph_get_perfect_path_with_first_edge
int db_graph_get_perfect_path_with_first_edge_all_colours(pathStep * first_step, void (*node_action) (dBNode * node), Path * path, dBGraph * db_graph)
{
	dBNode *node = first_step->node;
    
	Nucleotide fst_nucleotide = first_step->label;
	dBNode *current_node = node;
	Nucleotide nucleotide, nucleotide2;
	char tmp_seq[db_graph->kmer_size + 1];
	tmp_seq[db_graph->kmer_size] = '\0';
    
	//sanity check
	if (node == NULL) {
		printf("[db_graph_get_perfect_path_with_first_edge_all_colours] db_graph_get_perfect_path: can't pass a null node\n");
		exit(1);
	}
    
	path_reset(path);
    
	if (DEBUG) {
		printf("[db_graph_get_perfect_path_with_first_edge_all_colours] Node %i in path: %s\n",
               path->length, binary_kmer_to_seq(element_get_kmer(current_node), db_graph->kmer_size, tmp_seq));
	}
	//first edge defined
	nucleotide = fst_nucleotide;
	boolean added = true;
    
	pathStep current_step, next_step, rev_step;
	path_step_assign(&next_step, first_step);
	do {
		path_step_assign(&current_step, &next_step);
		added = path_add_node(&current_step, path);
		if (added) {
			node_action(current_step.node);
			//                      db_graph_get_next_node(pathStep * current_step,
			//                                      pathStep * next_step, pathStep * rev_step, dBGraph * db_graph);
			db_graph_get_next_step(&current_step, &next_step, &rev_step, db_graph);
            
			//sanity check
			if (next_step.node == NULL) {
				fprintf(stderr, "[db_graph_get_perfect_path_with_first_edge_all_colours] dB_graph: didnt find node in hash table: %s %c %s\n",
                        binary_kmer_to_seq(element_get_kmer(current_step.node), db_graph->kmer_size, tmp_seq),
                        binary_nucleotide_to_char(current_step.label),
                        current_step.orientation == forward ? "forward" : "reverse");
				exit(1);
			}
		} else {
			// If not added (we ran out of space), then remove the last node and make it's label undefined...
			pathStep ps;
			path_get_last_step(&ps, path);
			ps.label = Undefined;
			path_remove_last(path);
			path_add_node(&ps, path);
			if (DEBUG) {
				printf("[db_graph_get_perfect_path_with_first_edge_all_colours] Trying to correct last label failed.\n");
			}
		}
	}
	while (added && !path_is_cycle(path) &&	//loop
	       db_node_has_precisely_one_edge_all_colours(next_step.node, opposite_orientation(next_step.orientation), &nucleotide2) &&	//multiple entries
	       db_node_has_precisely_one_edge_all_colours(next_step.node, next_step.orientation, &next_step.label));	//has one next edge only
    
	if (current_step.node != NULL) {
		if (!(db_node_has_precisely_one_edge_all_colours(next_step.node, opposite_orientation(next_step.orientation), &nucleotide2) &&
              db_node_has_precisely_one_edge_all_colours(next_step.node,	next_step.orientation, &next_step.label)))
		{
			next_step.label = Undefined;
			if (DEBUG) {
				printf("[db_graph_get_perfect_path_with_first_edge_all_colours] Changing last label to 'N'.\n");
			}
		}
		added = path_add_node(&next_step, path);
		if (added) {
			node_action(next_step.node);
		}
	}
    
	if (DEBUG) {
		if (next_step.node == NULL) {
			printf("[db_graph_get_perfect_path_with_first_edge_all_colours] Next node is null! \n");
		} else {
			printf("[db_graph_get_perfect_path_with_first_edge_all_colours] Last node in path: %s %i length: %i\n",
                   binary_kmer_to_seq(element_get_kmer(next_step.node), db_graph->kmer_size, tmp_seq), db_node_get_edges_all_colours(next_step.node), path->length);
		}
	}
    
	return path_get_edges_count(path);
    
}


void db_graph_print_status(dBGraph * db_graph)
{
	log_and_screen_printf("dBGraph:\n unique kmers: %'lld\n", db_graph->unique_kmers);
	log_and_screen_printf(" Capacity: %'lld \n", (db_graph->bucket_size * db_graph->number_buckets));
	float cap = (float)db_graph->unique_kmers /(float) (db_graph->bucket_size * db_graph->number_buckets) * 100;
	log_and_screen_printf(" Occupied: %f%%\n", cap);
}


int db_graph_get_perfect_path(dBNode * node, Orientation orientation, void (*node_action) (dBNode * node), dBGraph * db_graph, Path * path)
{
    
	//sanity check
	if (node == NULL) {
		printf("[db_graph_get_perfect_path] can't pass a null node\n");
		exit(1);
	}
	if (path == NULL) {
		printf("[db_graph_get_perfect_path] can't pass a null path\n");
		exit(1);
	}
	if (db_graph == NULL) {
		printf("[db_graph_get_perfect_path] can't pass a null db_graph\n");
		exit(1);
	}
	perfect_path_get_path(node, orientation, node_action, db_graph, path);
    
	return path_get_edges_count(path);
    //return  path_get_length(path);
}



// limit is the max length
// min_coverage, max_coverage and avg_coveragte refer to the internal nodes
int db_graph_supernode(dBNode * node, void (*node_action) (dBNode * node), Path * path, dBGraph * db_graph)
{
	return db_graph_get_perfect_path(node, undefined, node_action, db_graph, path);
    
}

void printNode(dBNode * dbn, short int kmerSize)
{
	if (dbn != NULL) {
		char seq1[kmerSize];
		printf("MemAdd:\%p\n", dbn);
		printf("Node:\t%s \n",
		       binary_kmer_to_seq(&(dbn->kmer), kmerSize, seq1));
		printf("Flags:\t%x \n", dbn->flags);
		printf("Edges:\t%x \n", db_node_get_edges(dbn));
        int i;
        for (i = 0; i < NUMBER_OF_COLOURS; i++) {
            printf("Col %d:\t%d \n", i, dbn->coverage[i]);
        }
	} else {
		printf("NULL node");
	}
	return;
}


void db_graph_add_node_action(WalkingFunctions * wf, void (*node_action)(dBNode * node)){
    
    assert(wf);
    
    if (node_action == NULL) {
        return;
    }
    
    if(wf->node_callbacks.used >= MAX_STACKED_FUNCTIONS){
        fprintf(stderr, "Trying to overload more functions than we have capacity.");
        assert(wf->node_callbacks.used >= MAX_STACKED_FUNCTIONS);
        exit(-1);
    }
   
    wf->node_callbacks.callback[wf->node_callbacks.used] = NULL;
    
    wf->node_callbacks.callback[wf->node_callbacks.used++] = (void (*) ) node_action; //Make sure the action is not in the array yet.  
    
    
}

void db_graph_add_path_callback(WalkingFunctions * wf, void (*path_callback)(Path * path)){
    
    assert(wf != NULL);
    if (path_callback == NULL) {
        return;
    }
    
    if(wf->path_callbacks.used >= MAX_STACKED_FUNCTIONS){
        fprintf(stderr, "Trying to overload more functions than we have capacity.");
        assert(wf->path_callbacks.used >= MAX_STACKED_FUNCTIONS);
        exit(-1);
    }
    
    wf->path_callbacks.args[wf->path_callbacks.used] = NULL;
    wf->path_callbacks.callback[wf->path_callbacks.used++] = (void (*) ) path_callback; 
   
}

boolean db_graph_remove_path_callback(WalkingFunctions * wf, void * funct){
    int used_orig = wf->path_callbacks.used;
    int removed = 0;
    int i;
    
    for (i = 0; i < used_orig; i++) {
        if(wf->path_callbacks.callback[i] == funct) {
            removed++;
        }
        if (i + removed < MAX_STACKED_FUNCTIONS) {
            wf->path_callbacks.args[i] = wf->path_callbacks.args[i + removed];
            wf->path_callbacks.callback[i] = wf->path_callbacks.callback[i + removed];

        }else{
            wf->path_callbacks.args[i] = NULL;
            wf->path_callbacks.callback[i] = NULL;
        }
    }
    wf->path_callbacks.used -= removed;
    assert( wf->path_callbacks.used >= 0);
    return  removed > 0;
}

void db_graph_add_step_action(WalkingFunctions * wf, void (*step_action)(pathStep * ps)){
    
    assert(wf != NULL);
    if (step_action == NULL) {
        return;
    }
    
    if(wf->step_actions.used >= MAX_STACKED_FUNCTIONS){
        fprintf(stderr, "Trying to overload more functions than we have capacity.");
        assert(wf->step_actions.used >= MAX_STACKED_FUNCTIONS);
        exit(-1);
    }
    
    wf->step_actions.callback[wf->step_actions.used] = NULL;
    wf->step_actions.callback[wf->step_actions.used++] = (void (*) ) step_action; //Make sure the action is not in the array yet.  
    
}

void db_graph_add_node_action_with_args(WalkingFunctions * wf, void (*node_action)(dBNode * node, void * arg), void * args){
    
    assert(wf);
    
    if (node_action == NULL) {
        return;
    }
    
    if(wf->node_callbacks.used >= MAX_STACKED_FUNCTIONS){
        fprintf(stderr, "Trying to overload more functions than we have capacity.");
        assert(wf->node_callbacks.used >= MAX_STACKED_FUNCTIONS);
        exit(-1);
    }
    
    wf->node_callbacks.callback[wf->node_callbacks.used] = args;
    wf->node_callbacks.callback[wf->node_callbacks.used++] = node_action; //Make sure the action is not in the array yet.  
    
    
}

void db_graph_add_path_callback_with_args(WalkingFunctions * wf, void (*path_callback)(Path * path, void * arg), void * args){
    
    assert(wf != NULL);
    if (path_callback == NULL) {
        return;
    }
    
    if(wf->path_callbacks.used >= MAX_STACKED_FUNCTIONS){
        fprintf(stderr, "Trying to overload more functions than we have capacity.");
        assert(wf->path_callbacks.used >= MAX_STACKED_FUNCTIONS);
        exit(-1);
    }
    
    wf->path_callbacks.args[wf->path_callbacks.used] = args;
    wf->path_callbacks.callback[wf->path_callbacks.used++] = path_callback; 
    
}

void db_graph_add_step_action_with_args(WalkingFunctions * wf, void (*step_action)(pathStep * ps, void * arg), void * args){
    
    assert(wf != NULL);
    if (step_action == NULL) {
        return;
    }
    
    if(wf->step_actions.used >= MAX_STACKED_FUNCTIONS){
        fprintf(stderr, "Trying to overload more functions than we have capacity.");
        assert(wf->step_actions.used >= MAX_STACKED_FUNCTIONS);
        exit(-1);
    }
    
    wf->step_actions.callback[wf->step_actions.used] = args;
    wf->step_actions.callback[wf->step_actions.used++] = step_action; //Make sure the action is not in the array yet.  
    
}


static void execute_path_callbacks(Path * p, PathCallbackArray * callbacks){
    int i;
    for (i = 0; i < callbacks->used; i++) {
        void  (* f)() = callbacks->callback[i];
        if (callbacks->args[i] == NULL) {
            f(p);
        }else{
            f(p, callbacks->args[i]);
        }
    }  
}


static void execute_path_step_callbacks(pathStep * p, PathStepActionCallbackArray * callbacks){
    int i;
    for (i = 0; i < callbacks->used; i++) {
        void  (* f)() = callbacks->callback[i];
        if (callbacks->args[i] == NULL) {
            f(p);
        }else{
            f(p, callbacks->args[i]);
        }
    }  
}

static void execute_node_callbacks(dBNode * n, NodeActionCallbackArray * callbacks){
    int i;
    for (i = 0; i < callbacks->used; i++) {
        void  (* f)() = callbacks->callback[i];
        if (callbacks->args[i] == NULL) {
            f(n);
        }else{
            f(n, callbacks->args[i]);
        }
    }  
}

// clip a tip in the graph (the tip starts in node)
// limit is max length for tip
// node_action is applied to all the elements in the tip
// returns the length of the tip (0 means no length)
int db_graph_db_node_clip_tip(dBNode * node, int limit, void (*node_action) (dBNode * node), dBGraph * db_graph)
{
	int length_tip = 0;
    
	length_tip = db_graph_db_node_clip_tip_with_orientation(node, forward, limit, node_action, db_graph);
    
	if (length_tip == 0) {
        //printf("PLEASE NOTE: I do get here!\n");
		length_tip = db_graph_db_node_clip_tip_with_orientation(node, reverse, limit, node_action, db_graph);
	}
    
	return length_tip;
}

struct print_supernode_args{
    long long count_kmers ;
	long long count_sing ;
    int count_nodes;
    int max_length;
    boolean with_coverages;
    dBGraph * db_graph;
    Path * path;
    FILE *fout;
	FILE *fout_cov;
};

static 	void print_supernode(dBNode * node, void * v) {
    struct print_supernode_args * args = v;
    args->count_kmers++;
    if (db_node_check_flag_visited(node) == false) {
        int length = db_graph_supernode(node,
                                        &db_node_action_set_flag_visited,
                                        args->path, args->db_graph);
        
        if (length > 0) {
            if (args->with_coverages) {
                path_to_coverage(args->path, args->fout_cov);
            }
            
            path_to_fasta(args->path, args->fout);
            path_increase_id(args->path);
            if (length == args->max_length) {
                printf("contig length equals max length [%i] for node_%i\n", args->max_length, args->count_nodes);
            }
            args->count_nodes++;
            path_reset(args->path);
            
        } 
    }
}  

void db_graph_print_supernodes(char *filename, int max_length, boolean with_coverages, dBGraph * db_graph)
{
	
    struct print_supernode_args  args[1];
    
    args[0].max_length = max_length;
	args[0].fout = fopen(filename, "w");
    
	if (with_coverages) {
		char filename_cov[strlen(filename) + 10];
		sprintf(filename_cov, "%s_cov", filename);
		args[0].fout_cov = fopen(filename_cov, "w");
	}
    
	args[0].count_nodes = 0;
    
	args[0].path = path_new(max_length, db_graph->kmer_size);
	if (!args[0].path) {
		fprintf(stderr, "\n[db_graph_print_supernodes] Can't get memory for new path.\n\n");
		exit(-1);
	}
    
	hash_table_traverse_with_args(&print_supernode, (void**) &args, db_graph);
	log_and_screen_printf("%'qd nodes visited [%'qd singletons]\n", args[1].count_kmers, args[1].count_sing);
    
	path_destroy(args[0].path);
	fclose(args[1].fout);
	if (with_coverages) {
		fclose(args[1].fout_cov);
	}
}

typedef struct {
    long long tmp_cov[NUMBER_OF_COLOURS];
    dBGraph * db_graph;
} db_graph_stats_args;

static void get_cov(dBNode * node, void * dbsa){
    db_graph_stats_args * args = (db_graph_stats_args * ) dbsa;
    int i;
    boolean all_common;
    long long curr_cov;
    
    if(db_node_check_flag_not_pruned(node)){
        all_common = true;
        for(i = 0; i < NUMBER_OF_COLOURS; i++){
            curr_cov = element_get_coverage_by_colour(node, i);
            if(curr_cov){
                args->tmp_cov[i] += curr_cov;
                args->db_graph->colour_kmers[i]++;
            }else{
                all_common = false;
            }
        }
        if(all_common){
            args->db_graph->common_kmers_in_all_colours ++;
        }
    } 
}


void db_graph_calculate_stats(dBGraph * db_graph){
    //TODO: This should be easy to multithread... 
    
    if(db_graph->calculated == true){
        return;
    }
    double avg_cov = 0; 
    //	long long tmp_cov[NUMBER_OF_COLOURS];
    short i;
    //   long long curr_cov;
    db_graph_stats_args ** tmp_args = calloc(1, sizeof(db_graph_stats_args *));
    tmp_args[0] = calloc(1, sizeof(db_graph_stats_args));
    
    for(i = 0; i < NUMBER_OF_COLOURS; i++){
        
        tmp_args[0]->tmp_cov[i] = 0;
        db_graph->colour_kmers[i] = 0;
    }
    tmp_args[0]->db_graph = db_graph;
    
	log_and_screen_printf("Calculating graph stats...\n");
    
    db_graph->common_kmers_in_all_colours = 0;
    
    hash_table_traverse_with_args(&get_cov,(void **) tmp_args, db_graph);
    for(i = 0; i < NUMBER_OF_COLOURS; i++){
	    avg_cov = (double)tmp_args[0]->tmp_cov[i] / (double)db_graph->colour_kmers[i];
	    db_graph->average_coverage[i] = avg_cov;
	    
    }
    db_graph->calculated = true;
    free(tmp_args[0]);
    free(tmp_args);
    
}

//This just gets the average colour of the first colour
double db_graph_get_average_coverage(dBGraph * db_graph){
	if (db_graph->calculated == false) {//We already have calculated the stats, there is no need 
        db_graph_calculate_stats(db_graph);
    }
	
	return db_graph->average_coverage[0];
}

double db_graph_get_average_coverage_by_colour(short colour, dBGraph * db_graph){
	if (db_graph->calculated == false) {//We already have calculated the stats, there is no need 
        db_graph_calculate_stats(db_graph);
    }
	
	return db_graph->average_coverage[colour];
}


struct print_coverage_args{
    long long max_cov;
	long long *count;
};

static 	void f_max_cov(dBNode * node, void * v) {
    struct print_coverage_args * args = v;
    long long tmp_cov = element_get_coverage_all_colours(node);
    if (args->max_cov < tmp_cov) {
        args->max_cov = tmp_cov;
    }
}


static void sum_cov(dBNode * node, void * v) {
    struct print_coverage_args * args = v;
    long long tmp_cov = element_get_coverage_all_colours(node);
    args->count[tmp_cov - 1]++;
}

void db_graph_print_coverage(FILE * out, dBGraph * db_graph)
{
    
    struct print_coverage_args args[1];
	args[0].max_cov = 0;
	
    
	fprintf(stdout, "Getting all the coverages...\n");
    
	hash_table_traverse_with_args(&f_max_cov, (void **)&args ,db_graph);
    
	args[0].count = calloc(args[0].max_cov, sizeof(long long));
    
	if(args[0].count == NULL){
		fprintf(stderr, "[db_graph_print_coverage] WARNING! Unable to allocate memory to count the coverages (%lli)\n", args[0].max_cov);
		return;
	}
    
	
	hash_table_traverse_with_args(&sum_cov, (void **)&args, db_graph);
	fprintf(out, "Kmer\tCoverage\n");
    long long tmp_cov;
	for (tmp_cov = 0; tmp_cov < args[0].max_cov; tmp_cov++) {
		fprintf(out, "%lld\t%lld\n", tmp_cov + 1, args[0].count[tmp_cov]);
	}
    
	fprintf(stdout, "done \n");
	free(args[0].count);
}

struct kmer_cov_args{
    char *tmp_kmer;
    FILE * out; 
    dBGraph * db_graph;
};

static void print_cov(dBNode * node, void * v) {
    struct kmer_cov_args * args = v;
    
    binary_kmer_to_seq(element_get_kmer(node), args->db_graph->kmer_size, args->tmp_kmer);
    fprintf(args->out, "%s\t%d\n", args->tmp_kmer, element_get_coverage_all_colours(node));
}

void db_graph_print_kmer_coverage(FILE * out, dBGraph * db_graph)
{
    
    
    struct kmer_cov_args args[1];
	args[0].tmp_kmer = calloc(db_graph->kmer_size + 1, sizeof(char));
    args[0].db_graph = db_graph;
    args[0].out = out;
    
	fprintf(stdout, "Printing kmer and coverage..");
	fprintf(out, "Kmer\tCoverage\n");
    
	hash_table_traverse_with_args(&print_cov, (void **) &args, db_graph);
	fprintf(stdout, "done \n");
    
	free(args[0].tmp_kmer);
}

struct dump_binary_args {
    long long count;
    boolean(*condition) (dBNode * node);
    FILE * fout;
    dBGraph * db_graph;
    short colour;
};

//routine to dump graph
static void print_node_binary(dBNode * node, void * v) {
    struct dump_binary_args *  args = v;
    if (args->condition(node)) {
        args->count++;
        db_node_print_binary(args->fout, node, args->db_graph->kmer_size);
    }
}

// Should combine with below function
void db_graph_dump_binary(char *filename, boolean(*condition) (dBNode * node), dBGraph * db_graph)
{
	FILE *fout;		//binary output
	fout = fopen(filename, "w");
	if (fout == NULL) {
		fprintf(stderr, "cannot open %s", filename);
		exit(1);
	}
    
	int mean_read_len=0;//Mario - you can plumb this in if you want. See cortex_var/core/graph_info.c, and how the GraphInfo object is used in my file_reader * cprtex_var.c, 
	// for how I did it.
	long long total_seq=0;
	print_binary_signature(fout,db_graph->kmer_size,1, mean_read_len, total_seq);
    
    
	struct dump_binary_args * args[1];
    struct dump_binary_args arg;
    arg.db_graph = db_graph;
    arg.fout = fout;
    arg.condition = condition;
    arg.count = 0;
    args[0] = &arg;
    
	hash_table_traverse_with_args(&print_node_binary, (void **)&args,db_graph);
	fclose(fout);
    
	log_and_screen_printf("%'lld kmers dumped\n", arg.count);
}

//routine to dump graph
static void print_node_binary_by_colour(dBNode * node, void * v) {
    struct dump_binary_args *  args = v;
    if ((args->condition(node)) && (node->coverage[args->colour] > 0)) {
        args->count++;
        db_node_print_binary_by_colour(args->fout, node, args->colour, args->db_graph->kmer_size);
    }
}

// Should combine with above function
void db_graph_dump_binary_by_colour(char *filename, boolean(*condition) (dBNode * node), short colour, dBGraph * db_graph)
{
	FILE *fout;		//binary output
	fout = fopen(filename, "w");
	if (fout == NULL) {
		fprintf(stderr, "cannot open %s", filename);
		exit(1);
	}
    
	int mean_read_len=0;//Mario - you can plumb this in if you want. See cortex_var/core/graph_info.c, and how the GraphInfo object is used in my file_reader * cprtex_var.c, 
	// for how I did it.
	long long total_seq=0;
	print_binary_signature(fout,db_graph->kmer_size,1, mean_read_len, total_seq);
    
    
	long long count = 0;

    struct dump_binary_args args[1];

    args[0].db_graph = db_graph;
    args[0].fout = fout;
    args[0].condition = condition;
    args[0].count = 0;
    args[0].colour = colour;
	hash_table_traverse_with_args(&print_node_binary_by_colour, (void **) args,db_graph);
	fclose(fout);
    
	log_and_screen_printf("%qd kmers dumped\n", count);
}





int db_graph_db_node_clip_tip_with_orientation(dBNode * node, Orientation orientation, int limit, void (*node_action) (dBNode * node), dBGraph * db_graph)
{
	Nucleotide nucleotide, reverse_nucleotide;
	int length = 0;
	int i;
	dBNode *nodes[limit];
	Orientation next_orientation;
	dBNode *next_node;
	char seq[db_graph->kmer_size + 1];
    
	//starting in a blunt end also prevents full loops
	if (db_node_is_blunt_end_all_colours(node, opposite_orientation(orientation))) {
		boolean join_main_trunk = false;
        
		while (db_node_has_precisely_one_edge_all_colours(node, orientation, &nucleotide)) {
			nodes[length] = node;
            
			next_node = db_graph_get_next_node(node, orientation, &next_orientation, nucleotide, &reverse_nucleotide, db_graph);
            
			if (next_node == NULL) {
				printf("dB_graph_db_node_clip_tip_with_orientation: didnt find node in hash table: %s\n",
                       binary_kmer_to_seq(element_get_kmer(node), db_graph->kmer_size, seq));
				exit(1);
			}
            
			length++;
            
			if (length > limit) {
				break;
			}
			//want to stop when we join another trunk
			if (!db_node_has_precisely_one_edge_all_colours(next_node, opposite_orientation(next_orientation), &nucleotide) || !db_node_has_precisely_one_edge_all_colours(next_node, next_orientation, &nucleotide)) {
				join_main_trunk = true;
				break;
			}
			//keep track of node
            
			node = next_node;
			orientation = next_orientation;
		}
        
		if (!join_main_trunk) {
			length = 0;
		} else {	//clear edges and mark nodes as pruned
			for (i = 0; i < length; i++) {
                
				if (DEBUG) {
					printf("CLIPPING node: %s\n", binary_kmer_to_seq(element_get_kmer(nodes[i]), db_graph->kmer_size, seq));
				}
                
				node_action(nodes[i]);
                //perhaps we want to move this inside the node action?
				//(we already did). db_node_reset_edges_all_colours(nodes[i]);
			}
            
			if (DEBUG) {
				printf("RESET %c BACK\n", binary_nucleotide_to_char(reverse_nucleotide));
			}
			db_node_reset_edge_all_colours(next_node, opposite_orientation(next_orientation), reverse_nucleotide);//Why is this necessary?  
		}
        
	}
    
	return length;
}


#ifndef SOLID
//it doesn't check that it is a valid arrow -- it just assumes the arrow is fine
dBNode *db_graph_get_next_node(dBNode * current_node, Orientation current_orientation, Orientation * next_orientation,
                               Nucleotide edge, Nucleotide * reverse_edge, dBGraph * db_graph)
{
	if (edge == Undefined) {
		return NULL;	//If it is undefined, the next node is null...
	}
    
	BinaryKmer local_copy_of_kmer;
	binary_kmer_assignment_operator(local_copy_of_kmer, current_node->kmer);
    
	BinaryKmer tmp_kmer;
	dBNode *next_node = NULL;
    
	// after the following line tmp_kmer and rev_kmer are pointing to the same B Kmer
	BinaryKmer *rev_kmer = binary_kmer_reverse_complement(&local_copy_of_kmer, db_graph->kmer_size, &tmp_kmer);
    
	if (current_orientation == reverse) {
		*reverse_edge = binary_kmer_get_last_nucleotide(&local_copy_of_kmer);
		binary_kmer_assignment_operator(local_copy_of_kmer, *rev_kmer);
	} else {
		*reverse_edge = binary_kmer_get_last_nucleotide(rev_kmer);
	}
    
	binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&local_copy_of_kmer, edge, db_graph->kmer_size);
    
	//get node from table
	next_node = hash_table_find(element_get_key(&local_copy_of_kmer, db_graph->kmer_size, &tmp_kmer), db_graph);
    
	if (next_node != NULL) {
		*next_orientation = db_node_get_orientation(&local_copy_of_kmer, next_node, db_graph->kmer_size);
	}
    
	return next_node;
}
#else

dBNode *db_graph_get_next_node(dBNode * current_node, Orientation current_orientation,Orientation * next_orientation, Nucleotide edge, Nucleotide * reverse_edge, dBGraph * db_graph){
    if (edge == Undefined) {
		return NULL;	//If it is undefined, the next node is null...
	}
	BinaryKmer local_copy_of_kmer;
	binary_kmer_assignment_operator(local_copy_of_kmer, current_node->kmer);
    
	BinaryKmer tmp_kmer;
	dBNode *next_node = NULL;
    
	// after the following line tmp_kmer and rev_kmer are pointing to the same B Kmer
	BinaryKmer *rev_kmer = binary_kmer_reverse_complement(&local_copy_of_kmer, db_graph->kmer_size, &tmp_kmer);
    
	if (current_orientation == reverse) {        
        *reverse_edge = binary_kmer_get_first_nucleotide(&local_copy_of_kmer, db_graph->kmer_size);
        binary_kmer_right_shift_and_insert_new_base_at_left_end(&local_copy_of_kmer, edge, db_graph->kmer_size);
    } else {
        *reverse_edge = binary_kmer_get_last_nucleotide(rev_kmer);
        binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&local_copy_of_kmer, edge, db_graph->kmer_size);
	}
    
	
    
	//get node from table
	next_node = hash_table_find(element_get_key(&local_copy_of_kmer, db_graph->kmer_size, &tmp_kmer), db_graph);
    
	if (next_node != NULL) {
		*next_orientation = current_orientation;
	}
    
	return next_node;
}

#endif
int db_graph_generic_walk(pathStep * first_step, Path * path, WalkingFunctions * functions, dBGraph * db_graph)
{
	//dBNode * node = first_step->node;
	//sanity check
	if (first_step == NULL) {
		printf("[db_graph_generic_walk] can't pass a null node\n");
        assert(0);
		exit(1);
	}
    
	if (path == NULL) {
		printf("[db_graph_generic_walk] can't pass a null path\n");
        assert(0);
		exit(1);
	}
    
	if (functions == NULL) {
		printf("[db_graph_generic_walk] can't pass  null functions\n");
        assert(0);
		exit(1);
	}
    
	if (db_graph == NULL) {
		printf("[db_graph_generic_walk] can't pass  null db_graph\n");
        assert(0);
		exit(1);
	}
    
	pathStep current_step, next_step, rev_step;
	
	path_step_assign(&next_step, first_step);
    next_step.path = path; //This way, on all the assignments we keep the pointer to the path that was sent to the function originally. 
	functions->get_starting_step(&next_step, db_graph);
    
	boolean  try = true;
	int count = 0;
	boolean walked, added;
    
	do {
		walked = false;
        do {            
            added = path_add_node(&next_step, path);
            path_step_assign(&current_step, &next_step);
            try = false;
            functions->pre_step_action(&current_step);
            
            added = false;            
            if (current_step.label != Undefined) {
                functions->get_next_step(&current_step,&next_step, &rev_step,db_graph);
            }
            
            functions->step_action(&current_step);
            execute_path_step_callbacks(&current_step,&functions->step_actions);
            execute_node_callbacks(current_step.node, &functions->node_callbacks);
        }
        while (functions->continue_traversing(&current_step, &next_step, &rev_step, path,db_graph));
		
		if (path->length > 0) {
			walked = true;
		}
        
		if (walked) {
			count++;
			functions->output_callback(path);
            execute_path_callbacks(path, &functions->path_callbacks);
            
		}
		do {
			pathStep ps;
            
			path_get_last_step(&ps, path);
			if(ps.node != NULL){
                functions->post_step_action(&ps);
                path_remove_last(path);
            }
		}while (functions->continue_backwards(path, db_graph));
        if(path_get_length(path) > 0){//This is to enable the backtracking. 
            path_get_last_step(&current_step, path);//The continue_backwards should fix the last step
            functions->get_next_step(&current_step,&next_step, &rev_step,db_graph);
            
        }
	}while (path_get_length(path) > 0);
	return count;
}


struct print_graph_viz_args{
    dBGraph * db_graph;
    FILE * fp;
};


static void print_graphviz_colours(Element * node, void * v) {
    struct print_graph_viz_args * args = v;
    if (node != NULL) {
        BinaryKmer tmp, co;
        short kmer_size = args->db_graph->kmer_size;
        binary_kmer_assignment_operator(tmp, node->kmer);
        binary_kmer_reverse_complement(&tmp, kmer_size, &co);
        char seq[kmer_size], seqNext[kmer_size], seq1[kmer_size];
        binary_kmer_to_seq(&tmp, kmer_size, seq1);
        char *print = db_node_check_for_any_flag(node, STARTING_FORWARD  | BRANCH_NODE_FORWARD | BRANCH_NODE_REVERSE | END_NODE_FORWARD | END_NODE_REVERSE | X_NODE) ? "ellipse" : "ellipse";
        
#ifdef ENABLE_BUBBLEPARSE
        char *node_colour = ((node->coverage[0] == 0) && (node->coverage[1] > 0)) ? "orange":(((node->coverage[0] > 0) && (node->coverage[1] == 0)) ? "green" : "black");
#else
        char *node_colour = "black";
#endif
        if (db_node_check_for_any_flag(node, BRANCH_NODE_FORWARD | BRANCH_NODE_REVERSE | X_NODE)) {
            print = "circle";
        }
        
        char reverse_seq1[kmer_size + 1];
        char label_string[512];
        
        seq_reverse_complement(seq1, kmer_size, reverse_seq1);
        //sprintf(label_string, "\"%s\\n%s\\nC0=%d C1=%d\"", seq1, reverse_seq1, node->coverage[0], node->coverage[1]);
        
#ifdef ENABLE_BUBBLEPARSE
        sprintf(label_string,
                "<<table border=\"0\" cellpadding=\"0\" cellspacing=\"0\">"
                "<tr><td><font color=\"blue\">%s</font></td></tr>"
                "<tr><td><font color=\"red\">%s</font></td></tr>"
                "<tr><td><font color=\"green\">%d</font> <font color=\"orange\">%d</font></td></tr>"
                "</table>>", seq1, reverse_seq1,
                node->coverage[0], node->coverage[1]);
#else
        sprintf(label_string,
                "<<table border=\"0\" cellpadding=\"0\" cellspacing=\"0\">"
                "<tr><td><font color=\"blue\">%s</font></td></tr>"
                "<tr><td><font color=\"red\">%s</font></td></tr>"
                "<tr><td><font color=\"black\">%d</font> </td></tr>"
#ifdef ENABLE_READ_PAIR
                "<tr><td><font color=\"blue\">0x%llX</font></td></tr>"
                "<tr><td><font color=\"red\"> 0x%llX </font></td></tr>"
#endif				
                
#ifdef SOLID
                "<tr><td><font color=\"blue\">%c</font></td></tr>"
                "<tr><td><font color=\"red\">%c</font></td></tr>"
                
#endif
                "</table>>", seq1, reverse_seq1,
                node->coverage[0]
#ifdef ENABLE_READ_PAIR
                , (long long unsigned int)db_node_get_signature(0,0,node),
                (long long unsigned int)db_node_get_signature(1,0,node)
#endif	
#ifdef SOLID
                ,binary_nucleotide_base_space_to_char(db_node_get_starting_base(forward, node))
                ,binary_nucleotide_base_space_to_char(db_node_get_starting_base(reverse, node))
#endif
                
                );
#endif
        fprintf(args->fp, "%s [label=%s, shape=%s, color=%s]\n", seq1, label_string, print, node_colour);
        
        binary_kmer_to_seq(&tmp, kmer_size, seq);
        binary_kmer_left_shift(&tmp, 2, kmer_size);
        Orientation o = forward;
        //db_node_print_to_stream(node, args->db_graph, stderr);

        for (o = 0; o <=1; o++) {

            Edges ea[NUMBER_OF_COLOURS];
            int i;
            
            for (i = 0; i < NUMBER_OF_COLOURS; i++) {
                ea[i] = db_node_get_edges_for_orientation_by_colour(node, o, i);
            }
            
            Nucleotide n;
            for (n = 0; n < 4; n++) {
                BinaryKmer bk;
                Key k = &bk;
                boolean colour0edge = (((ea[0] >> n) & 1) == 1);
                boolean colour1edge = false;
               // fprintf(stderr, "%c %d %d %d\n", binary_nucleotide_to_char(n), o, colour0edge, colour1edge );
#ifdef ENABLE_BUBBLEPARSE
                colour1edge = (((ea[1] >> n) & 1) == 1);
#endif
                if (colour0edge || colour1edge) {
                    char *labelcolour = (colour0edge && colour1edge) ? "black":(colour0edge ? "green":"orange");
                    binary_kmer_modify_base(&tmp, n, kmer_size, 0);
                    binary_kmer_to_seq(element_get_key(&tmp, kmer_size, k), kmer_size, seqNext);
                                        fprintf(args->fp,
                            "%s -> %s [ label=<<font color=\"%s\">%c</font>>, color=\"%s\"];\n",
                            seq, seqNext, labelcolour,
                            binary_nucleotide_to_char(n),
                            o == forward ? "blue" : "red");
                }
            }
            
            binary_kmer_assignment_operator(tmp, co);
            binary_kmer_left_shift(&tmp, 2, kmer_size);
        }
        
    }
}

void db_graph_write_graphviz_file(char *filename, dBGraph * db_graph)
{
	FILE *fp;
    
	
    
	if (db_graph->unique_kmers < 450) {
		if (filename) {
			fp = fopen(filename, "w");
		} else {
			fp = stderr;
		}
		if (fp) {
			fprintf(fp, "digraph finite_state_machine {\n");
			fprintf(fp, "\trankdir=LR;\n");
			fprintf(fp, "\tsize=\"12,8\";\n");
			fprintf(fp, "\tfontsize=18;\n");
			fprintf(fp, "\tnode [shape = circle];\n");
            struct print_graph_viz_args ** args = calloc(1, sizeof(struct print_graph_viz_args*));
            args[0] = calloc(1, sizeof(struct print_graph_viz_args));
            args[0]->db_graph = db_graph;
            args[0]->fp = fp;
			hash_table_traverse_with_args(&print_graphviz_colours, (void **)args, db_graph);
            free(args[0]);
            free(args);
			fprintf(fp, "}\n");
			fclose(fp);
		}
		printf("\nWritten graphviz file %s.\n\n", filename);
	} else {
		printf("\nNote: Too many kmers to output graphviz file.\n\n");
	}
}

// Build array of neighbouring nodes
int db_graph_get_neighbouring_nodes_all_colours(dBNode* start, Orientation orientation, dBNode* neighbours[], Nucleotide labels[], dBGraph * db_graph) {
    int num_neighbours = 0;
    Nucleotide nucleotide;
    for (nucleotide = Adenine; nucleotide < Undefined; 
         nucleotide++) {
        if (db_node_edge_exist_any_colour(start, nucleotide, orientation)) {
            pathStep current_step, next_step, rev_step;
            current_step.node = start;
            current_step.label = nucleotide;
            current_step.orientation = orientation;
            db_graph_get_next_step(&current_step, &next_step, &rev_step, db_graph);
            
            if (next_step.node != NULL) {
                labels[num_neighbours] = nucleotide;
                neighbours[num_neighbours] = next_step.node;
                num_neighbours++;
            }
        }
    }
    
    
    return num_neighbours;
}

Nucleotide db_graph_get_best_next_step_nucleotide(dBNode * from, dBNode * previous,  Orientation orientation, dBGraph * db_graph ){
    assert(from != NULL);
    assert(previous != NULL);
    Nucleotide best = Undefined;
    
#if defined (ENABLE_MARK_PAIR) || defined (ENABLE_READ_PAIR)
    dBNode * neighbours[4];
    Nucleotide labels[4];
    int found = db_graph_get_neighbouring_nodes_all_colours(from, orientation, neighbours, labels, db_graph);
    int i;
#endif
#ifdef ENABLE_READ_PAIR    
    int j;
    //TODO: add all the read pairs, at the moment I will assume that we have coverage in the first pair. 
    
    
    // boolean valid = false;
    int  best_matches = 0, temp_matches = 0, second_best_matches = 0;
    ReadPairSignature rps = 0;
    for(j = 0; j < NUMBER_OF_SIGNATURES; j++){
        rps |= from->signatures[j];//TODO: make this line cleaner 
    }
    if (read_pair_count_bits(rps)  > 2 * (READ_PAIR_LENGTH / 3 ) ) {
        return Undefined;
    }
    
    for (i = 0; i < found; i++) {
        temp_matches = 0;
        for(j = 0; j < NUMBER_OF_SIGNATURES; j++){
            temp_matches += read_pair_count_common_bits(previous->signatures[j], neighbours[i]->signatures[j]);//TODO: make this line cleaner 
        }
        if(temp_matches > db_graph->rpda->pair[0]->min_pair_coverage && temp_matches > best_matches){ //Hardcoded 2 matches, may be a variable? Make this even more flexible... 
            second_best_matches = best_matches;
            best_matches = temp_matches;
            best = labels[i];
        }else if(temp_matches == best_matches){
            best = Undefined;
        }
    }
    if(best_matches - second_best_matches < db_graph->rpda->pair[0]->min_pair_coverage ){
        best = Undefined;
    }
    
#elif ENABLE_MARK_PAIR
    uint32_t conn_count = 0 ;
    uint32_t curr_count;
    boolean valid = true;
    uint32_t supernode_prev = db_node_get_supernode(previous);
    uint32_t supernode_curr = 0;
    if (supernode_prev == 0) {
        valid = false;
    }
    for (i = 0; i < found && valid; i++) {
        supernode_curr = db_node_get_supernode(neighbours[i]);                                            
        if (supernode_curr != 0) {
            curr_count = mark_pair_connections_between(supernode_prev, supernode_curr, db_graph->supernode_link);
            if (conn_count && curr_count) {
                valid = false;  
            }
            if (conn_count == 0 && curr_count ) {
                conn_count = curr_count;
            }
        }
    }
    if (valid == false) {
        best = Undefined;
    }
#else
    //TODO: Make this with the coverage. 
    
#endif
    
    return best;
}



void  db_node_print_to_stream(dBNode * node, dBGraph * graph, FILE * f){
    binary_kmer_to_stream(&node->kmer, graph->kmer_size, f);    
    
    Orientation o;
    Nucleotide n;
  
    
    for ( o=forward; o<=reverse; o++) {
        
        for (n = Adenine; n < Undefined; n++) {
            if (db_node_edge_exist(node, n, o)) {
                fprintf(f, "%c",binary_nucleotide_to_char(n));
            }else{
                fprintf(f, " ");
            }
        }
        fprintf(f, "-");
    }
    fprintf(f, "\n");
    
}


#ifdef DEBUG_CLEANUP
void db_graph_cleanup_graph(dBGraph * db_graph)
{
	Orientation orientation;
    
	void clean_node(dBNode * node) {
		void check_edge(Nucleotide nucleotide) {
			// Check if edge exists
			if (db_node_edge_exist_any_colour(node, nucleotide, orientation)) {
				// Now see if kmer it leads to exists
				pathStep current_step, next_step, rev_step;
				current_step.node = node;
				current_step.label = nucleotide;
				current_step.orientation = orientation;
				db_graph_get_next_step(&current_step, &next_step, &rev_step, db_graph);
				if (next_step.node == NULL) {
					printf("Missing node\n");
					db_node_reset_edge_all_colours(node, orientation, nucleotide);
				}
			}
		}
        
		orientation = forward;
		nucleotide_iterator(&check_edge);
        
		orientation = reverse;
		nucleotide_iterator(&check_edge);
	}
    
	printf("\nTidying up graph...\n");
	hash_table_traverse(&clean_node, db_graph);
}
#endif



