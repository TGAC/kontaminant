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
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <limits.h>
#include <assert.h>
#ifdef THREADS
#include <pthread.h>
#endif

#include <global.h>
#include <flags.h>
#include <nucleotide.h>
#include <seq.h>
#include <binary_kmer.h>
#include <element.h>
#include <path.h>
#include <math.h>
#include <logger.h>

Path *path_new(int max_length, short kmer_size)
{
	Path *path = calloc(1, sizeof(Path));
	if(path == NULL){
			fprintf(stderr, "[path_new]Unable to allocate path\n");
			exit(-1);
	}
	path->nodes = calloc(max_length, sizeof(dBNode *));
	path->orientations = calloc(max_length, sizeof(Orientation));
	path->labels = calloc(max_length + 1, sizeof(Nucleotide));
	path->seq = calloc(max_length + 1 + kmer_size, sizeof(char));
	path->in_nodes = calloc(PATH_MAX_IN_NODES, sizeof(int));
    path->step_flags = calloc(max_length, sizeof(Flags));
    path->max_virtual_length = max_length;
	if(path->nodes == NULL){
		fprintf(stderr, "[path_new]Unable to allocate nodes\n");
		exit(-1);
	}
	
	if(path->orientations == NULL){
		fprintf(stderr, "[path_new]Unable to allocate orientations\n");
		exit(-1);
	}
	
	if(path->labels == NULL){
		fprintf(stderr, "[path_new]Unable to allocate labels\n");
		exit(-1);
	}
	
	if(path->seq == NULL){
		fprintf(stderr, "[path_new]Unable to allocate seq\n");
		exit(-1);
	}
	
	if(path->in_nodes == NULL){
		fprintf(stderr, "[path_new]Unable to allocate in_nodes\n");
		exit(-1);
	}
	
	path->in_nodes_capacity = PATH_MAX_IN_NODES;
    path->in_nodes_count = 0;
	path->max_length = max_length;
	flags_action_unset_flag(IS_CYCLE, &(path->flags));
	flags_action_set_flag(PRINT_FIRST, &(path->flags));
	path->kmer_size = kmer_size;
	path->header = 0;
	path->used = false;

#ifdef ENABLE_MARK_PAIR
    path->supernodes_total_count = 0;
    path->supernodes = NULL;
    path->supernodes_count = NULL;
#endif
	
    if ((!path->nodes) || (!path->orientations) || (!path->labels)
	    || (!path->seq)) {
		free(path);
		path = 0;
	}

	return path;
}

void path_destroy(Path * path)
{
	if (path == NULL) {
		fprintf(stderr, "[path_destroy] The path is null");
		exit(-1);
	}
	free(path->nodes);
	free(path->orientations);
	free(path->labels);
	free(path->seq);
	free(path->in_nodes);
    free(path->step_flags);

	if (path->header) {
		free(path->header);
	}

	free(path);

}

void path_increase_id(Path * path)
{
	path->id++;
}

//We can only set the limit to new paths! 
void path_set_limit(int limit, Path * p){
    assert(limit <= p->max_length);
    assert(p->length == 0);
    p->max_virtual_length = limit;
}

int path_get_limit(Path * p){
    return p->max_virtual_length;
}

void path_reset(Path * path)
{
	if (path == NULL) {
		fprintf(stderr, "[path_reset] The path is null");
		exit(-1);
	}
	while (path->length > 0) {
		path_remove_last(path);
	}

	path->seq[0] = '\0';
	path->in_nodes_count = 0;
    path->out_nodes_count = 0;
	path->length = 0;
    path->new_nodes = 0;
    path->max_virtual_length = path->max_length;
	flags_action_clear_flags(&(path->flags));
	flags_action_set_flag(PRINT_FIRST, &(path->flags));
	flags_action_set_flag(NEW_PATH, &(path->flags));
	
	path_clean_stop_reason(path);
	
	int i;
		
	for(i = 0; i < path->in_nodes_capacity; i++){
		path->in_nodes[i] = 0;
	}
#ifdef ENABLE_MARK_PAIR
    path->supernodes_total_count = 0;
#endif
}

/**
 * To get a natural order of two steps. 
 * At the moment we are using just the address to compare, is
 * cheaper than comparing the kmers. If the nodes are the same, 
 * we use the orientations to distinguish them, and if Use with care. 
 */
int path_step_compare(pathStep * a, pathStep * b){
	if (a->node < b->node) {
		return -1;
	}else if (a->node > b->node) {
		return 1;
	}else if (a->orientation < b->orientation) {
		return -1;
	}else if (a->orientation > b->orientation) {
		return 1;
	}else if (a->label < b->label) {
		return -1;
	}else if (a->label > b->label) {
		return 1;
	}else {
		return 0;
	}

}

boolean path_step_has_unvisited_edge_all_colours(pathStep * step){
    Edges e = db_node_get_edges_for_orientation_all_colours(step->node, step->orientation);
    Edges visited = step->flags & PATH_STEP_ALL_VISITED; //We relay on being using the 4 least significan bit.
    return e != visited;
}



Nucleotide path_step_get_unvisited_edge_all_colours(pathStep * step){
    Edges e = db_node_get_edges_for_orientation_all_colours(step->node, step->orientation);
    Edges visited = step->flags & PATH_STEP_ALL_VISITED; //We relay on being using the 4 least significan bit.
    Nucleotide n, nucleotide = Undefined; 
    Edges edges = ~visited & e;
    for (n = 0; n < 4 && nucleotide == Undefined; n++) {
		if ((edges & 1) == 1) {
			nucleotide = n;
		}
		edges >>= 1;
	}
    
    return nucleotide;
}

boolean path_step_equals(pathStep * step, pathStep * other)
{
	return step->node == other->node
	    && step->orientation == other->orientation
	    && step->label == other->label;
}

// Compare steps, without comparing labels
boolean path_step_equals_without_label(pathStep * step, pathStep * other)
{
	return step->node == other->node
	    && step->orientation == other->orientation;
}

void path_modfy_last_label(Nucleotide n, Path * p){
    int last_index = path_get_length(p)-1;
    char last_c = binary_nucleotide_to_char(n);
    p->labels[last_index] = n;
    p->seq[last_index] = last_c;
    p->step_flags[last_index] |= binary_nucleotide_to_edge(n);
    
}

boolean unlabelled_path_step_equals(pathStep * step, pathStep * other)
{
	return step->node == other->node
	    && step->orientation == other->orientation;
}


boolean path_contains_step(pathStep * step, Path * path)
{    
    pathStep tmp;
    int i;
    int len = path_get_length(path);
    for (i = 0; i < len; i++) {
        path_get_step_at_index(i, &tmp, path);
        if (path_step_equals(step, &tmp)) {
			return true;
		}
    }
	return false;
}

boolean path_contains_step_without_label(pathStep * step, Path * path)
{
    pathStep tmp;
    int i;
    int len = path_get_length(path);
    for (i = 0; i < len; i++) {
        path_get_step_at_index(i, &tmp, path);
        if (path_step_equals_without_label(step, &tmp)) {
			return true;
		}
    }
	return false;
}


boolean path_is_first_equals_to_last(Path * path)
{
	if(path->length < 2)
		return false;
	
	pathStep first, last;
	path_get_last_step(&last, path);
	path_get_step_at_index(0, &first, path);
	return path_step_equals_without_label(&first, &last);
}

boolean path_contains_step_from_index(int index, pathStep * step, Path * path)
{
	pathStep tmp;
    int i;
    int len = path_get_length(path);
    for (i = index; i < len; i++) {
        path_get_step_at_index(i, &tmp, path);
        if (path_step_equals(step, &tmp)) {
			return true;
		}
    }
	return false;
}

boolean path_contains_step_from_index_without_label(int index, pathStep * step, Path * path)
{
    pathStep tmp;
    int i;
    int len = path_get_length(path);
    for (i = index; i < len; i++) {
        path_get_step_at_index(i, &tmp, path);
        if (path_step_equals_without_label(step, &tmp)) {
			return true;
		}
    }
	return false;
}

boolean path_has_space(Path * path)
{
	if (path->length >= path->max_virtual_length) {
		return false;
	}
	return true;

}

void path_add_in_step(pathStep * step, Path * path)
{
	assert(path != NULL);
	assert(step != NULL);
	assert(step->node != NULL);
	assert(path->in_nodes != NULL);
	assert(path->in_nodes_count < path->in_nodes_capacity);
	
	path->in_nodes[path->in_nodes_count++] = path->length;
	
	if (path->in_nodes_count == path->in_nodes_capacity) {
		int new_capacity = path->in_nodes_capacity + PATH_IN_NODES_CAPACITY_INCREASE;
		int * in_tmp_ptr = realloc(path->in_nodes, new_capacity * sizeof(int));
        int i;
		if (in_tmp_ptr == NULL){
			fprintf(stderr, "out of memory while increasing the number  to %d of in paths  in path: \n", 
                        new_capacity);
			path_to_fasta(path, stderr);
			exit(1);
        }
        
        path->in_nodes = in_tmp_ptr;
        
        for (i=path->in_nodes_count; i<new_capacity; i++) {
            path->in_nodes[i] = 0;
        }
        
        path->in_nodes_capacity = new_capacity;
	}
	
	assert(path->in_nodes_count <= path->in_nodes_capacity);
}

boolean path_has_in_step(pathStep * step, Path * path){
	boolean found = false;
	pathStep tmp_step;
	int i;
	for(i = 0; i < path->in_nodes_count && !found; i++ ){
			path_get_step_at_index(path->in_nodes[i], &tmp_step, path);
			if(path_step_equals_without_label(&tmp_step, step)){
				found = true;
			}
		}
		
	return found;
}

int path_index_of_last_in_node(Path * path){
    int i = path->in_nodes_count; 
    if( i !=  0){
        return path->in_nodes[i-1];
    }
    return -1;

}

int path_get_first_in_node_after(int pos, Path * path){
    assert(pos >= 0);
    int i, arr_pos = -1;
    for (i = path->in_nodes_count-1; i > 0; i--) {
        if(path->in_nodes[i] > pos){
            arr_pos = i;
        }
    }
    return arr_pos != -1? path->in_nodes[arr_pos]: -1;
}


boolean path_add_node(pathStep * step, Path * path)
{
	assert(path!=NULL);
	assert(step != NULL);
	assert(step->node != NULL);
	
	if (step->node == NULL) {
		fprintf(stderr, "[path_add_node] The node is null\n");
		//exit(-1);
		return false;	//If you try to add a null node, you dont add it. Currently never happens because the assert... 
	}
	if (path->length == path->max_virtual_length) {
		//log_and_screen_printf("[path_add_node] The node_%qd has reached the limit[%d]\n",
		//	     path->id, path->max_virtual_length);
        //path_to_fasta(path, <#FILE *fout#>)
		return false;
	}

	dBNode *node = step->node;
	if (DEBUG) {
		char tmp_seq[path->kmer_size + 1];
		tmp_seq[path->kmer_size] = '\0';
		printf
		    ("[path_add_node] Node %i in path(%lli): %s add %c (%s).\n",
		     path->length, path->id,
		     binary_kmer_to_seq(element_get_kmer(node), path->kmer_size,
					tmp_seq),
		     binary_nucleotide_to_char(step->label),
		     step->orientation == reverse ? "reverse" : "forward");
	}
	Orientation orientation = step->orientation;
	Nucleotide nucleotide = step->label;
	pathStep first;
	if(path->length > 0){
		path_get_step_at_index(0, &first, path);
		if(path_step_equals_without_label(&first, step)){ //TODO this is not enough for a cycle as the cycle doesn't need to close on first node (it is also disregarding the orientation)
			flags_action_set_flag(IS_CYCLE, &(path->flags));
		}
		if(path_has_in_step(step, path)){
			flags_action_set_flag(IS_CYCLE, &(path->flags));
		}
		
		
	}
    int edges_in = db_node_edges_count_all_colours(step->node,  opposite_orientation(step->orientation));
    if(edges_in > 1){
        path_add_in_step(step, path);
    }
    if(!db_node_check_flag_visited(step->node)){
        path->new_nodes++;
    }
    
    int edges_out = db_node_edges_count_all_colours(step->node, step->orientation);
    if(edges_out > 1){
        path->out_nodes_count++;
    }
    
	/*flags_action_set_flag(FIND_AGAIN, &(path->flags)); */
	flags_action_unset_flag(NEW_PATH | EMPTY_PATH, &(path->flags));

	path->nodes[path->length] = node;
	path->orientations[path->length] = orientation;
	path->labels[path->length] = nucleotide;
    if(nucleotide != Undefined){
        path->step_flags[path->length] = step->flags | binary_nucleotide_to_edge(nucleotide);//We add the current nucleotide to the flags. The step may come with already marked nucleotides. 
	}
    path->seq[path->length] = nucleotide == Undefined ? '\0'
	    : binary_nucleotide_to_char(nucleotide);
	path->seq[path->length + 1] = '\0';
	
#ifdef ENABLE_MARK_PAIR
    path_add_supernode(step->node->supernode, path);
#endif
    
	path->length++;

	return true;

}

void path_step_print(pathStep * step, int kmer_size, FILE * f)
{
	char tmp_seq[kmer_size + 1];
	tmp_seq[kmer_size] = '\0';
	fprintf(f, "pathStep: %s %c (%s).",
		binary_kmer_to_seq(element_get_kmer(step->node), kmer_size,
				   tmp_seq),
		binary_nucleotide_to_char(step->label),
		step->orientation == reverse ? "reverse" : "forward");
}

Nucleotide path_last_nucleotide(Path * path)
{
	return (path->length <= 0) ? Undefined : path->labels[path->length - 1];
}

Flags path_last_flags(Path * path)
{
	return (path->length <= 0) ? 0 : path->step_flags[path->length - 1];
}

static void compute_label(dBNode * node, Orientation o, char *label)
{
	int i = 0;
	Nucleotide n;
    for (n = Adenine; n < Undefined; n++) {
        if (db_node_edge_exist(node, n, o)) {
			label[i] = binary_nucleotide_to_char(n);
			i++;
		}
    }
    
	label[i] = '\0';
}

static void check_print_first(Path * path)
{
	boolean ignoreFirst = false;
#ifndef SHORT_FLAGS
	Flags fs = db_node_get_flags(path->nodes[0], PRINT_FORWARD | PRINT_REVERSE);
	// Include the first node if:
	// - first node orientation is forward and PRINT_FORWARD is set
	// - first node orientation is reverse and PRINT_REVERSE is set	
	ignoreFirst = (path->orientations[0] == forward) ? ((fs == PRINT_FORWARD) ? false : true) : ((fs == PRINT_REVERSE) ? false : true);
#endif
	if (!ignoreFirst) {
		flags_action_set_flag(PRINT_FIRST, &(path->flags));
	} else {
		flags_action_unset_flag(PRINT_FIRST, &(path->flags));
	}
}


/*
 * Returns a path array with all the perfect paths from the path. If the path is a perfect path, it returns an array containing a copy of the path. 
 * WARNING! dont forget to use the free from buffer functions, and destroy the path array 
 */
PathArray * path_split_in_perfect_paths(Path * p){
    assert(p != NULL);
    assert(path_get_length(p) > 0);
    PathArray * pa =    path_array_new(1);
    int path_length = path_get_length(p);
    int i, count_fwd, count_rev;
    Path * current_path = path_get_buffer_path();
    path_array_add_path(current_path, pa);
    pathStep tmp;
    i = 0;
    path_get_step_at_index(i, &tmp, p);
    count_rev = db_node_edges_count_all_colours(tmp.node, reverse);
    path_add_node(&tmp, current_path);  
    if(count_rev > 1){
        current_path = path_get_buffer_path();
        path_add_node(&tmp, current_path);    
        path_array_add_path(current_path, pa);
    }
    for(i = 1; i < path_length; i++){
        path_get_step_at_index(i, &tmp, p);
        count_fwd =db_node_edges_count_all_colours(tmp.node, forward);
        count_rev =db_node_edges_count_all_colours(tmp.node, reverse);
        path_add_node(&tmp, current_path);
        if(count_fwd > 1 || count_rev > 1 ){
            current_path = path_get_buffer_path();
            path_add_node(&tmp, current_path);    
            path_array_add_path(current_path, pa);
        }
    }
    return pa;
}
//Merges the first array into the second one, and frees the original array, but it doesnt destroy the paths inside. Careful with the leaks
void path_array_merge(PathArray ** from, PathArray * to){
    PathArray * from_p = *from;
    int paths = path_array_get_number_of_paths(from_p);
    int i;
    for(i = 0; i < paths; i++){
         path_array_add_path(path_array_get(i, from_p), to);
    }
    path_array_destroy_struct(from);
   
    
}

void path_to_fasta_debug(Path * path, FILE * fout)
{
	short kmer_size = path->kmer_size;
	//int length = path->length;
    
	// Sanity checking
	if (path == NULL) {
		printf("[path_to_fasta] trying to print a null Path\n");
		exit(-1);
	}
    
	if (fout == NULL) {
		printf("[path_to_fasta] trying to print to a null FILE\n");
		exit(-1);
	}
    
	if (DEBUG) {
		printf("[path_to_fasta] About to print a path\n");
	}
	// Set value of PRINT_FIRST (whether first node included)
	//check_print_first(path);
    
	// Get coverage statistics from path
	double avg_coverage;
	int min_coverage;
	int max_coverage;
	path_get_statistics(&avg_coverage, &min_coverage, &max_coverage, path);
    
	// Get orientation of first and last node
	Orientation fst_orientation;
	//if (flags_check_for_flag(PRINT_FIRST, &(path->flags))) {
    fst_orientation = path->orientations[0];
	//} else {
	//	fst_orientation = path->orientations[1];
	//}
	Orientation lst_orientation = path->orientations[path->length];
    
	// Get the first node - this will be nodes[0] if PRINT_FIRST is
	// specified, or nodes[1] otherwise.
	dBNode *fst_node;
	if (flags_check_for_flag(PRINT_FIRST, &(path->flags))) {
		if (path->length == 0) {
			printf("[path_to_fasta] Trying to print an empty path[1]!\n");
			return;
		}
		fst_node = path->nodes[0];
	} else {
		if (path->length < 2) {
			fprintf(stderr, "[path_to_fasta] Trying to print an empty path[2]!\n");
			return;
		}
		fst_node = path->nodes[1];
	}
    
	// Get the last node
	dBNode *lst_node = path->nodes[path->length - 1];
    
	// Make a set of labels for first and last nodes which list the
	// acceptable forward and reverse path labels
	char fst_f[5], fst_r[5], lst_f[5], lst_r[5];
	compute_label(fst_node, forward, fst_f);
	compute_label(fst_node, reverse, fst_r);
	compute_label(lst_node, forward, lst_f);
	compute_label(lst_node, reverse, lst_r);
    
	// Places to store first and last kmer sequence
	char fst_seq[kmer_size + 1], lst_seq[kmer_size + 1];
	fst_seq[kmer_size] = '\0';
	lst_seq[kmer_size] = '\0';
	BinaryKmer tmp_kmer;
    
	// Get the first kmer sequence (or complement)
	BinaryKmer fst_kmer;
	binary_kmer_assignment_operator(fst_kmer, *(element_get_kmer(fst_node)));
	if (fst_orientation == reverse) {
		binary_kmer_reverse_complement(&fst_kmer, kmer_size, &tmp_kmer);
		binary_kmer_assignment_operator(fst_kmer, tmp_kmer);
	}
	binary_kmer_to_seq(&fst_kmer, kmer_size, fst_seq);
	if (DEBUG) {
		printf("[path_to_fasta] First kmer: %s\n", fst_seq);
    }
    
	// Get the last kmer sequence (or complement)
	BinaryKmer lst_kmer;
	binary_kmer_assignment_operator(lst_kmer, *(element_get_kmer(lst_node)));
	if (lst_orientation == reverse) {
		binary_kmer_reverse_complement(&lst_kmer, kmer_size, &tmp_kmer);
		binary_kmer_assignment_operator(lst_kmer, tmp_kmer);
	}
	binary_kmer_to_seq(&lst_kmer, kmer_size, lst_seq);
    
	// Output to file
	/*fprintf(fout,
            ">node_%qd length:%i average_coverage:%.2f min_coverage:%i max_coverage:%i fst_coverage:%i fst_kmer:%s fst_r:%s fst_f:%s lst_coverage:%i lst_kmer:%s lst_r:%s lst_f:%s\n",
            path->id,
            (flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? length +
             kmer_size : length + kmer_size - 1), avg_coverage,
            min_coverage, max_coverage,
            element_get_coverage_all_colours(fst_node), fst_seq,
            (fst_orientation == forward ? fst_r : fst_f),
            (fst_orientation == forward ? fst_f : fst_r),
            element_get_coverage_all_colours(lst_node), lst_seq,
            (lst_orientation == forward ? lst_r : lst_f),
            (lst_orientation == forward ? lst_f : lst_r));
    */
	binary_kmer_to_seq(&fst_kmer, flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? kmer_size : kmer_size - 1, fst_seq);
    
	int i, current;
	for(i = 0, current = 1; i < path->kmer_size; i++, current++) {
		fprintf(fout, "%c", fst_seq[i]);
		if(current % PATH_FASTA_LINE == 0) {
			fprintf(fout, "\n");	
		}
	}
	
	size_t len = strlen(path->seq);
	for(i = 0; i < len; i++, current++){
        if (path->step_flags[i] & PRINT_LABEL_AS_N) {
            fprintf(fout, "N");
        } else if (path->step_flags[i] & PRINT_LABEL_LOWERCASE) {
            fprintf(fout, "%c", tolower(path->seq[i]));
        } else {
            fprintf(fout, "%c",  path->seq[i]);
        }
		if(current % PATH_FASTA_LINE == 0){
			fprintf(fout, "\n");	
		}
	}
	fprintf(fout, "\n");
	
	fflush(fout);
}

void path_to_fasta(Path * path, FILE * fout)
{
	short kmer_size = path->kmer_size;
	int length = path->length;
    int max_length = path->max_length;
	// Sanity checking
	if (path == NULL) {
		fprintf(stderr,	"[path_to_fasta] trying to print a null Path\n");
		exit(-1);
	}

	if (fout == NULL) {
		fprintf(stderr,
			"[path_to_fasta] trying to print to a null FILE\n");
		exit(-1);
	}

	if (DEBUG) {
		printf("[path_to_fasta] About to print a path\n");
	}
    
    if (length == max_length) {
        log_and_screen_printf
        ("contig length equals max length [%i] for node_%i\n",
         max_length, path->id);
    }
	// Set value of PRINT_FIRST (whether first node included)
	//check_print_first(path);

	// Get coverage statistics from path
	double avg_coverage;
	int min_coverage;
	int max_coverage;
	path_get_statistics(&avg_coverage, &min_coverage, &max_coverage, path);

	// Get orientation of first and last node
	Orientation fst_orientation;

    fst_orientation = path->orientations[0];
	
	Orientation lst_orientation = path->orientations[path->length];

	// Get the first node - this will be nodes[0] if PRINT_FIRST is
	// specified, or nodes[1] otherwise.
	dBNode *fst_node;
	if (flags_check_for_flag(PRINT_FIRST, &(path->flags))) {
		if (path->length == 0) {
			fprintf(stderr, "[path_to_fasta] Trying to print an empty path[1]!\n");
			return;
		}
		fst_node = path->nodes[0];
	} else {
		if (path->length < 2) {
			fprintf(stderr,	"[path_to_fasta] Trying to print an empty path[2]!\n");
			return;
		}
		fst_node = path->nodes[1];
	}

	// Get the last node
	dBNode *lst_node = path->nodes[path->length - 1];

	// Make a set of labels for first and last nodes which list the
	// acceptable forward and reverse path labels
	char fst_f[5], fst_r[5], lst_f[5], lst_r[5];
	compute_label(fst_node, forward, fst_f);
	compute_label(fst_node, reverse, fst_r);
	compute_label(lst_node, forward, lst_f);
	compute_label(lst_node, reverse, lst_r);

	// Places to store first and last kmer sequence
	char fst_seq[kmer_size + 1], lst_seq[kmer_size + 1];
	fst_seq[kmer_size] = '\0';
	lst_seq[kmer_size] = '\0';
	BinaryKmer tmp_kmer;

	// Get the first kmer sequence (or complement)
	BinaryKmer fst_kmer;
	binary_kmer_assignment_operator(fst_kmer, *(element_get_kmer(fst_node)));
	if (fst_orientation == reverse) {
		binary_kmer_reverse_complement(&fst_kmer, kmer_size, &tmp_kmer);
		binary_kmer_assignment_operator(fst_kmer, tmp_kmer);
	}
	binary_kmer_to_seq(&fst_kmer, kmer_size, fst_seq);
	if (DEBUG) {
		printf("[path_to_fasta] First kmer: %s\n", fst_seq);
    }
    
	// Get the last kmer sequence (or complement)
	BinaryKmer lst_kmer;
	binary_kmer_assignment_operator(lst_kmer, *(element_get_kmer(lst_node)));
	if (lst_orientation == reverse) {
		binary_kmer_reverse_complement(&lst_kmer, kmer_size, &tmp_kmer);
		binary_kmer_assignment_operator(lst_kmer, tmp_kmer);
	}
	binary_kmer_to_seq(&lst_kmer, kmer_size, lst_seq);

	// Output to file
	fprintf(fout,
            ">node_%qd length:%i average_coverage:%.2f min_coverage:%i max_coverage:%i fst_coverage:%i fst_r:%s fst_f:%s lst_coverage:%i lst_r:%s lst_f:%s\n",
            path->id,
            (flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? length + kmer_size : length + kmer_size - 1), avg_coverage,
            min_coverage,
            max_coverage,
            element_get_coverage_all_colours(fst_node),
            (fst_orientation == forward ? fst_r : fst_f),
            (fst_orientation == forward ? fst_f : fst_r),
            element_get_coverage_all_colours(lst_node),
            (lst_orientation == forward ? lst_r : lst_f),
            (lst_orientation == forward ? lst_f : lst_r));

	binary_kmer_to_seq(&fst_kmer, flags_check_for_flag(PRINT_FIRST,	&(path->flags)) ? kmer_size : kmer_size - 1, fst_seq);

	int i, current= 1;
    
//If we are doing solid, we print in color space. 
#ifdef SOLID 

    fprintf(fout, "%c",  binary_nucleotide_base_space_to_char(path_get_first_base(path)));

#endif
    
	for(i = 0, current = 1 ; i < path->kmer_size; i++, current++) {
		fprintf(fout, "%c", fst_seq[i]);
		if(current % PATH_FASTA_LINE == 0) {
			fprintf(fout, "\n");	
		}
	}
	
	size_t len = strlen(path->seq);
	for(i = 0; i < len; i++, current++) {
        if (path->step_flags[i] & PRINT_LABEL_AS_N) {
            fprintf(fout, "N");
        } else if (path->step_flags[i] & PRINT_LABEL_LOWERCASE) {
            fprintf(fout, "%c", tolower(path->seq[i]));
        } else {
            fprintf(fout, "%c",  path->seq[i]);
        }
        
		if(current % PATH_FASTA_LINE == 0) {
			fprintf(fout, "\n");	
		}
	}
	fprintf(fout, "\n");
	
	fflush(fout);
}

// Cloned from path_to_fasta with changes for file colour coding
// NOTE: Has not been updated for PRINT_LABEL_AS_N flag
void path_to_fasta_colour(Path * path, FILE * fout, char *id)
{
	short kmer_size = path->kmer_size;

	// Sanity checking
	if (path == NULL) {
		fprintf(stderr,	"[path_to_fasta_colour] trying to print a null Path\n");
		exit(-1);
	}

	if (fout == NULL) {
		fprintf(stderr,	"[path_to_fasta_colour] trying to print to a null FILE\n");
		exit(-1);
	}

	if (DEBUG) {
		printf("[path_to_fasta_colour] About to print a path\n");
	}
	// Set value of PRINT_FIRST (whether first node included)
	//check_print_first(path);

	// Get coverage statistics from path
	double avg_coverage;
	int min_coverage;
	int max_coverage;
	path_get_statistics(&avg_coverage, &min_coverage, &max_coverage, path);

	// Get orientation of first and last node
	Orientation fst_orientation;
	if (flags_check_for_flag(PRINT_FIRST, &(path->flags))) {
		fst_orientation = path->orientations[0];
	} else {
		fst_orientation = path->orientations[1];
	}
	Orientation lst_orientation = path->orientations[path->length];

	// Get the first node - this will be nodes[0] if PRINT_FIRST is
	// specified, or nodes[1] otherwise.
	dBNode *fst_node;
	if (flags_check_for_flag(PRINT_FIRST, &(path->flags))) {
		if (path->length == 0) {
			fprintf(stderr,	"[path_to_fasta_colour] Trying to print an empty path[1]!\n");
			return;
		}
		fst_node = path->nodes[0];
	} else {
		if (path->length < 2) {
			fprintf(stderr,	"[path_to_fasta_colour] Trying to print an empty path[2]!\n");
			return;
		}
		fst_node = path->nodes[1];
	}

	// Get the last node
	dBNode *lst_node = path->nodes[path->length - 1];

	// Make a set of labels for first and last nodes which list the
	// acceptable forward and reverse path labels
	char fst_f[5], fst_r[5], lst_f[5], lst_r[5];
	compute_label(fst_node, forward, fst_f);
	compute_label(fst_node, reverse, fst_r);
	compute_label(lst_node, forward, lst_f);
	compute_label(lst_node, reverse, lst_r);

	// Places to store first and last kmer sequence
	char fst_seq[kmer_size + 1], lst_seq[kmer_size + 1];
	fst_seq[kmer_size] = '\0';
	lst_seq[kmer_size] = '\0';
	BinaryKmer tmp_kmer;

	// Get the first kmer sequence (or complement)
	BinaryKmer fst_kmer;
	binary_kmer_assignment_operator(fst_kmer, *(element_get_kmer(fst_node)));
	if (fst_orientation == reverse) {
		binary_kmer_reverse_complement(&fst_kmer, kmer_size, &tmp_kmer);
		binary_kmer_assignment_operator(fst_kmer, tmp_kmer);
	}
	binary_kmer_to_seq(&fst_kmer, kmer_size, fst_seq);

	if (DEBUG) {
		printf("[path_to_fasta_colour] First kmer: %s\n", fst_seq);
    }
    
	// Get the last kmer sequence (or complement)
	BinaryKmer lst_kmer;
	binary_kmer_assignment_operator(lst_kmer, *(element_get_kmer(lst_node)));
	if (lst_orientation == reverse) {
		binary_kmer_reverse_complement(&lst_kmer, kmer_size, &tmp_kmer);
		binary_kmer_assignment_operator(lst_kmer, tmp_kmer);
	}
	binary_kmer_to_seq(&lst_kmer, kmer_size, lst_seq);

	// Output to file
	fprintf(fout,
            ">%s length:%i %s average_coverage:%.2f min_coverage:%i max_coverage:%i fst_coverage:%i fst_r:%s fst_f:%s lst_coverage:%i lst_r:%s lst_f:%s\n",
            id,
            (flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? ((int)strlen(path->seq) + kmer_size) : ((int)strlen(path->seq) + kmer_size - 1)),
            path->header,
            avg_coverage,
            min_coverage,
            max_coverage,
            element_get_coverage_all_colours(fst_node),
            (fst_orientation == forward ? fst_r : fst_f),
            (fst_orientation == forward ? fst_f : fst_r),
            element_get_coverage_all_colours(lst_node),
            (lst_orientation == forward ? lst_r : lst_f),
            (lst_orientation == forward ? lst_f : lst_r));

	fprintf(fout, "%s", binary_kmer_to_seq(&fst_kmer, kmer_size, fst_seq));
	fprintf(fout, "%s\n", flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? path->seq : path->seq + 1);
		
    fflush(fout);

log_and_screen_printf("path_to_fasta: PRINT_FIRST %i path->seq length %i path->length %i\n", flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? 1:0, strlen(path->seq), path->length);
}

void path_to_coverage(Path * path, FILE * fout)
{
	if (path == NULL) {
		fprintf(stderr,	"[path_to_coverage] trying to print a null Path");
		exit(-1);
	}

	if (fout == NULL) {
		fprintf(stderr,	"[path_to_coverage] trying to print to a null FILE");
		exit(-1);
	}

	check_print_first(path);
	
    int i = flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? 0 : 1;

	fprintf(fout, ">node_%qd \n", path->id);

	for (; i < path->kmer_size - 1; i++) {
		fprintf(fout, "%i ", element_get_coverage_all_colours(path->nodes[flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? 0 : 1]));
	}

	for (i = 0; i < path->length; i++) {
		fprintf(fout, "%i ", element_get_coverage_all_colours(path->nodes[i]));
	}
    
	fprintf(fout, "\n");
}

// Clone of path_to_coverage for handling colour coded files
void path_to_coverage_colour(Path * path, FILE * fout, char *id, short colour)
{
	short kmer_size = path->kmer_size;

	// Sanity checking
	if (path == NULL) {
		fprintf(stderr,	"[path_to_coverage_colour] trying to print a null Path");
		exit(-1);
	}

	if (fout == NULL) {
		fprintf(stderr,	"[path_to_coverage_colour] trying to print to a null FILE");
		exit(-1);
	}
	// Set value of PRINT_FIRST (whether first node included)
	//check_print_first(path);
	int i = flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? 0 : 1;

	// If colour 0, write the label
	if (colour == 0) {
		fprintf(fout, ">%s\n", id);
	}
	// Write first kmer coverage
	for (; i < kmer_size - 1; i++) {
		fprintf(fout, "%i ", element_get_coverage_by_colour(path->nodes[flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? 0 :	1], colour));
	}

	// Write path coverage
	i = flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? 0 : 1;
	for (; i < path->length; i++) {
		fprintf(fout, "%i ", element_get_coverage_by_colour(path->nodes[i], colour));
	}

	fprintf(fout, "\n");

log_and_screen_printf("path_to_coverage: PRINT_FIRST %i path->length %i\n", flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? 1:0, path->length);
}

void path_iterator_from_index(int index, void (*step_action) (pathStep * step), Path * path)
{
	int i = 0;
	pathStep ps;

	if (index >= path->length) {
		printf("[path_iterator_from_index] Error: index (%d) greater than path length (%d)\n", index, path->length);
		exit(1);
	}

	for (i = index; i < path->length; i++) {
		ps.node = path->nodes[i];
		ps.label = path->labels[i];
		ps.orientation = path->orientations[i];
		step_action(&ps);
	}
}

void path_inner_iterator(void (*step_action) (pathStep * step), Path * path)
{

	int i;
	pathStep ps;

	for (i = 1; i < path->length - 1; i++) {
		ps.node = path->nodes[i];
		ps.label = path->labels[i];
		ps.orientation = path->orientations[i];
		step_action(&ps);
	}
}


void path_iterator(void (*step_action) (pathStep * step), Path * path)
{
	int i = 0;
	pathStep ps;
	for (i = 0; i < path->length; i++) {
		ps.node = path->nodes[i];
		ps.label = path->labels[i];
		ps.orientation = path->orientations[i];
		step_action(&ps);
	}
}

void path_iterator_with_args(void (*step_action) (pathStep * , void *),void * args,  Path * path)
{
	int i = 0;
	pathStep ps;
	for (i = 0; i < path->length; i++) {
		ps.node = path->nodes[i];
		ps.label = path->labels[i];
		ps.orientation = path->orientations[i];
		step_action(&ps, args);
	}
}

void path_inner_iterator_with_args(void (*step_action) (pathStep * step, void*),void * args,  Path * path)
{
    
	int i;
	pathStep ps;
    
	for (i = 1; i < path->length - 1; i++) {
		ps.node = path->nodes[i];
		ps.label = path->labels[i];
		ps.orientation = path->orientations[i];
		step_action(&ps,args);
	}
}


boolean path_is_singleton(int length, Path * path){
    boolean sing = false;
    double avg_cov;
    int min_cov, max_cov;
    pathStep first;
    pathStep last;
    path_get_step_at_index(0, &first, path);
    path_get_last_step(&last, path);
    
    if(path_get_length(path) < length){
        if(path_is_blunt(forward, path ) && path_is_blunt(reverse, path)){
            sing = true;
        }else if(path_get_length(path) == 1 && (path_is_blunt(forward, path ) || path_is_blunt(reverse, path))){
            sing = true;
        }else if(first.node == last.node){
            sing = true;
        }else{
            path_get_statistics(&avg_cov, &min_cov, &max_cov, path);
            if (avg_cov  < 2) {
                sing = true;
            }
        }
        
    }
    
    return sing;
}


boolean path_is_repetitive(double graph_cov, Path * p)
{
    double avg_cov = 0;
    int min_cov = 0;
    int max_cov = 0;
    boolean rep= false;
    if(path_has_stop_reason(FIRST, PATH_FLAG_DIVERGING_PATHS, p) &&path_has_stop_reason(LAST, PATH_FLAG_DIVERGING_PATHS, p) ){
        if(path_get_length(p)  < p->kmer_size * 2){
            rep = true;
        }else{
            path_get_statistics(&avg_cov, &min_cov, &max_cov, p);
            if(avg_cov > (graph_cov * 3) ){
                rep = true;
            }
        }

    }
    
    if(path_get_length(p)  < p->kmer_size * 2){
        if(avg_cov > (graph_cov * 3) ){
            rep = true;
        }
    }
    
    return rep;

}

void path_iterator_reverse(void (*step_action) (pathStep * step), Path * path)
{
	int i = 0;
	pathStep ps;
	for (i = path->length - 1; i >= 0; i--) {
		ps.node = path->nodes[i];
		ps.label = path->labels[i];
		ps.orientation = path->orientations[i];
		step_action(&ps);
	}
}

void path_iterator_with_index(void (*step_action) (int index, pathStep * step), Path * path)
{
	int i = 0;
	pathStep ps;
	for (i = 0; i < path->length; i++) {
	
		ps.node = path->nodes[i];
	
		ps.label = path->labels[i];
	
		ps.orientation = path->orientations[i];
	
		step_action(i,&ps);
	
	}
}



/**
 *
 *Gets the index of the step in the path, or -1 if it
 * is not found.  
 */
int path_index_of_step(pathStep * step, Path * path)
{

//TODO fill it.         
	return -1;
}

pathStep *path_get_last_step_reverse(pathStep * step, Path * path)
{
	return path_get_step_reverse(step, path, path->length - 1);
}

// Reverse a step at a given index
pathStep *path_get_step_reverse(pathStep * step, Path * path, int index)
{
	assert (path != NULL) ;
	assert(step!=NULL);
	//assert(path->length != 1);
	assert(path->length >=1);
	assert(index < path->length);
	
	step->node = path->nodes[index];
	step->orientation = opposite_orientation(path->orientations[index]);    
    step->flags = path->step_flags[index];
    
    if (path->length == 1) {
        step->label = Undefined;
    } else {    
        if (index == 0) {
            step->label = Undefined;
        } else {
            BinaryKmer second_to_last_kmer;
            binary_kmer_assignment_operator(second_to_last_kmer,
                            path->nodes[index - 1]->kmer);

            if (path->orientations[index - 1] == reverse) {
                step->label =
                    binary_kmer_get_last_nucleotide
                    (&second_to_last_kmer);
            } else {
                BinaryKmer tmp_kmer;
                BinaryKmer *rev_kmer =
                    binary_kmer_reverse_complement(&second_to_last_kmer,
                                   path->kmer_size,
                                   &tmp_kmer);
                step->label = binary_kmer_get_last_nucleotide(rev_kmer);
            }
        }
    }
    step->path = path;
	return step;
}

pathStep *path_get_step_at_index(int index, pathStep * step, Path * path)
{
	assert(path != NULL);
	assert(step != NULL);
	assert(index < path->length);
	if (path == NULL) {
		fprintf(stderr,
			"[path_get_step_at_index] passing a  null Path");
		exit(-1);
	}
	if (step == NULL) {
		fprintf(stderr,
			"[path_get_step_at_index] passing a  pathStep ");
		exit(-1);
	}
	if (index >= path->length) {
		fprintf(stderr,
			"[path_get_step_at_index] The queried index (%d) is greater than the path length (%d)",
			index, path->length);
		exit(-1);
	}
	step->node = path->nodes[index];
	step->orientation = path->orientations[index];
	step->label = path->labels[index];
    step->flags = path->step_flags[index];
    step->path = path;
	return step;
}

char *path_get_seq(char *tmp, Path * path)
{
	strcpy(tmp, path->seq);
	return tmp;
}

void path_graphviz_open_header(FILE * f){
	fprintf(f, "graph{\n"
	"graph[page=\"8.5,11\",size=\"7.5,7\",ratio=fill,center=1];\n"
	"fontsize=18;\n"
	"node [shape = circle];\n"
	"edge [dir=none];\n");

	
}

void path_graphviz_close_header(FILE * f){
	fprintf(f, "}\n");

}

void path_graphviz_line(FILE * f, Path * p){
	
	if(p == NULL || p->length == 0)
		return;
	
	pathStep  first ;
	pathStep last;
	
	path_get_step_at_index(0, &first, p);
	path_get_last_step(&last, p);
	
	short kmer_size = p->kmer_size;
	char  seq_f[kmer_size + 1];
	char  seq_l[kmer_size + 1];
	BinaryKmer bk;
	Key k = &bk;
	BinaryKmer * tmp = element_get_kmer(first.node) ;
	binary_kmer_to_seq(element_get_key(tmp ,kmer_size, k), kmer_size, seq_f);
	
	tmp = element_get_kmer(last.node) ;
	binary_kmer_to_seq(element_get_key(tmp ,kmer_size, k), kmer_size, seq_l);
	double size = ((double)p->length/(100000));
	//size =  size;
	//double weight  = fabs(log(1/size));
	fprintf(f, "%lld [ shape=point, color=red];\n", p->id);
	fprintf(f, "%s [ shape=point, color=blue];\n",seq_f);
	fprintf(f, "%s [ shape=point, color=blue];\n",seq_l);
	fprintf(f, "%s -- %lld [len=\"%f\"];\n", seq_f, p->id, size);
	fprintf(f, "%lld -- %s [len=\"%f\"];\n", p->id, seq_l, size);
}


static void path_print_contig_with_details(FILE * fout, Path * p){
    
    pathStep  ps;
    int i;
    fprintf(fout, ">node_%qd \n", p->id);
    for(i=0; i<path_get_length(p); i++){
        path_get_step_at_index(i, &ps, p);
     //   fprintf(fout, );
    }

}

void path_get_statistics(double *avg_coverage, int *min_coverage, int *max_coverage, Path * path)
{

	/**
	 * TODO: validate if the path is empty... think about singletons....
	 */
	int i = flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? 0 : 1;
	*max_coverage = 0;
	*min_coverage = INT_MAX;
	int sum_coverage = 0;

	for (; i < path->length; i++) {	//Calculate the return values for the current path.

		int coverage = element_get_coverage_all_colours(path->nodes[i]);
		sum_coverage += coverage;
		*max_coverage =
		    (*max_coverage < coverage) ? coverage : *max_coverage;
		*min_coverage =
		    (*min_coverage > coverage) ? coverage : *min_coverage;

	}
	int length = path_get_nodes_count(path);
	*avg_coverage = (double)sum_coverage / (double)(length);

	if (*min_coverage == INT_MAX) {
		*min_coverage = 0;
	}

}

dBNode *path_last_node(Path * path)
{
	return path->length > 0 ? path->nodes[path->length - 1] : NULL;
}

int path_get_length(Path * path){
	assert(path!=NULL);
	return path->length;
}

pathStep *path_get_last_step(pathStep * ps, Path * path)
{
    assert(path->length != 0);
	ps->node = path_last_node(path);
	ps->orientation = path_last_orientation(path);
	ps->label = path_last_nucleotide(path);
    ps->flags = path_last_flags(path);
	return ps;
}

boolean path_is_empty(Path * path)
{
	//return flags_check_for_flag(EMPTY_PATH, &(path->flags));
	assert(path != NULL);
	return path->length == 0;
}

/**
 * Returns if a path is blunt. 
 * If orientation is forward, it checks if the last node is blunt
 * If the orientation is reverse, it checks if the first node is blunt. 
 */ 
boolean path_is_blunt(Orientation o, Path * p){
    if (path_get_length(p) > 0) {
        pathStep ps;
        
        if(o == forward){
            path_get_last_step(&ps, p);
        }else if(o == reverse){
            path_get_step_reverse(&ps, p, 0);
        }else{
            fprintf(stderr, "Invalid orientation at path_is_blunt\n");
            exit(-1);
        }
        
        return db_node_is_blunt_end(ps.node, ps.orientation);
    }else{
        return  false;
    }
	
}

void path_add_stop_reason(PathEnd o, Flags f, Path * path){
    assert(o == FIRST || o == LAST);
	if(o == FIRST){
		flags_action_set_flag(f, &(path->stop_reasons_first));
	}else if(o == LAST){
		flags_action_set_flag(f, &(path->stop_reasons_last));
	}
}

boolean path_has_stop_reason(PathEnd o, Flags f, Path * path){
    assert(o == FIRST || o == LAST);
	return o == FIRST?flags_check_for_flag(f, &(path->stop_reasons_first)): flags_check_for_flag(f, &(path->stop_reasons_last));
}

boolean path_has_any_stop_reason(PathEnd o, Flags f, Path * path){
    assert(o == FIRST || o == LAST);
	return o == FIRST?flags_check_for_any_flag(f, &(path->stop_reasons_first)): flags_check_for_flag(f, &(path->stop_reasons_last));
}

void path_clean_stop_reason(Path * path){
	
	path->stop_reasons_first = 0;
	path->stop_reasons_last = 0;
		//flags_action_clear_flags(&(path->stop_reasons_first));
	//flags_action_clear_flags(&(path->stop_reasons_last));
}



boolean path_is_cycle(Path * path)
{
	
	/*if(path->len < 2){
		return false;	
	}
	if(!flags_check_for_flag(IS_CYCLE, &(path->flags))){
		int len = path->length;
		int i,j;
		pathStep psi, psj;
		
		for(i = 0; i < len; i++){
			path_get_step_at_index(i, &psi, path);
		
			for(j = i+i; j < len; j++){
				path_get_step_at_index(j, &psj, path);
				
				if(path_step_equals_without_label(&psi, &psj)){
					i = len;
					j = len;
					flags_action_set_flag(IS_CYCLE, &(path->flags));
					//if (DEBUG) {
						fprintf(stdout,
								"[path_add_node] The node_%qd is cycle!\n",
								path->id);
					//}
				}
			}
		}
			
	}*/
	return flags_check_for_flag(IS_CYCLE, &(path->flags));
}

Orientation path_last_orientation(Path * path)
{
	if(path->length == 0)
		return undefined;
	return path->orientations[path->length - 1];
}

int path_free_spaces(Path * path)
{
	return path->max_length - path->length;
}

void path_remove_last(Path * path)
{
	if (path->length == 0) {
		printf
		    ("[path_remove_last] Trying to remove node from the node_%qd which is empty\n",
		     path->id);
		flags_action_set_flag(ERROR_PATH, &(path->flags));
		flags_action_unset_flag(FIND_AGAIN, &(path->flags));
        assert(path->length > 0);
		return;		//TODO: use the flags to tell the caller the invalid state
	}

	if (DEBUG) {
		char tmp_seq[path->kmer_size + 1];
		tmp_seq[path->kmer_size] = 0;	//'\0';
		printf
		    ("[path_remove_last] Removing node %i in path(%lld): %s %c\n",
		     path->length - 1, path->id,
		     binary_kmer_to_seq(element_get_kmer(path_last_node(path)),
					path->kmer_size, tmp_seq),
		     path->seq[path->length - 1] ==
		     0 ? 'N' : path->seq[path->length - 1]);
		//printNode(current_node, db_graph->kmer_size);
	}
    
    if(db_node_edges_count_all_colours(path_last_node(path), path_last_orientation(path)) > 1){
        path->out_nodes_count--;
        assert(path->out_nodes_count >=0);
    }
    
	path->length--;
#ifdef ENABLE_MARK_PAIR
    path_add_supernode(path->nodes[path->length]->supernode, path);
#endif
	path->seq[path->length] = '\0';
	path->nodes[path->length] = NULL;
	path->orientations[path->length] = 0;
	path->labels[path->length] = 0;
    path->step_flags[path->length] = 0;
	if ((path->in_nodes > 0) && (path->in_nodes_count > 0)) {
		if (path->in_nodes[path->in_nodes_count-1] == path->length){
			path->in_nodes[path->in_nodes_count-1] = 0;
            path->in_nodes_count--;
		}
	}
	flags_action_set_flag(FIND_AGAIN, &(path->flags));
	if (path->length == 0) {
		//flags_action_clear_flags(&(path->flags));
		//flags_action_set_flag(PRINT_FIRST, &(path->flags));
		//flags_action_set_flag(NEW_PATH, &(path->flags));
		flags_action_set_flag(EMPTY_PATH, &(path->flags));
		flags_action_unset_flag(FIND_AGAIN, &(path->flags));
	}
}

void path_do_nothing(Path * p)
{
}

boolean path_to_retry(Path * path)
{

	if (DEBUG) {
		printf("[path_to_retry]Flags: %x\n", path->flags);
	}
	return flags_check_for_any_flag(FIND_AGAIN | NEW_PATH, &(path->flags))
	    && !flags_check_for_any_flag(ERROR_PATH, &(path->flags));
}

void path_step_assign(pathStep * to, pathStep * from)
{
	to->label = from->label;
	to->node = from->node;
	to->orientation = from->orientation;
    to->path = from->path;
    to->flags = from->flags;
}

PathArray *path_array_new(short number_of_paths)
{
	PathArray *pa = calloc(1, sizeof(PathArray));
#ifdef THREADS
    pthread_mutex_init(&pa->mutex, NULL);
#endif
	if (pa == NULL) {
		fprintf(stderr,
			"[path_array_new] Not enough memory to allocate the PathArray");
		exit(-1);
	}
	pa->number_of_paths = 0;
	pa->capacity = number_of_paths;
	pa->paths = calloc(number_of_paths, sizeof(Path *));
	if (pa->paths == NULL) {
		fprintf(stderr,
			"[path_array_new] Not enough memory to allocate the Paths for the PathArray ");
	}
	return pa;
}

int path_array_get_number_of_paths(PathArray * pa){
    return pa->number_of_paths;
}

void path_array_destroy(PathArray * pa)
{
    
    //TODO: design a cleaver lock. 
    while (pa->number_of_paths > 0) {
		pa->number_of_paths--;
		path_destroy(pa->paths[pa->number_of_paths]);

	}
	free(pa->paths);
	free(pa);
}

void path_array_destroy_struct (PathArray ** pa){
    free((*pa)->paths);
    free((*pa));
    (*pa)=NULL;
}

//WARNING: This is not thread safe, the method calling it must be. 
static void  path_array_double_capacity(PathArray * pa){

    int new_capacity = pa->capacity * 2;
    Path ** new_array = realloc(pa->paths, new_capacity * sizeof(Path * ));
    if(new_array == NULL){
        fprintf(stderr, "[path_array_double_capacity] Unable to double the size of the PathArray (size %d)", new_capacity);
        assert(new_array != NULL);
        exit(-1);
    }
    pa->capacity = new_capacity;
    pa->paths = new_array;
    
    
    
}

boolean path_array_add_path(Path * p, PathArray * pa)
{
	if (pa == NULL) {
		fprintf(stderr, "[path_array_add_path] PathArray is null");
		exit(-1);
	}

	if (p == NULL) {
		fprintf(stderr, "[path_array_add_path] Path is null");
		exit(-1);
	}

	if (pa->number_of_paths >= pa->capacity) {
		path_array_double_capacity(pa);
	}
	pa->paths[pa->number_of_paths] = p;
	pa->number_of_paths++;
	return true;
}

void path_action_clear_flags(Path * node)
{
	flags_action_clear_flags(&(node->flags));
}

void path_action_set_flag(Path * node, Flags f)
{
	flags_action_set_flag(f, &(node->flags));

}

void path_action_unset_flag(Path * node, Flags f)
{
	flags_action_unset_flag(f, &(node->flags));
}

Flags path_get_flags(Path * node, Flags f)
{
	return node->flags & f;
}

boolean path_check_for_flag(Path * node, Flags flag)
{
	return flags_check_for_flag(flag, &(node->flags));
    
}

boolean path_check_for_any_flag(Path * path, Flags flag)
{
	return flags_check_for_any_flag(flag, &(path->flags));
}

int path_percentage_new_nodes(Path * path){
    return (100 * path->new_nodes)/path->length;
    
}

void path_reverse(Path * source, Path * destination)
{
	pathStep new_step;
	int i;

	for (i = source->length - 1; i >= 0; i--) {
		path_get_step_reverse(&new_step, source, i);
        new_step.flags = new_step.flags & PATH_STEP_MASK_VISITED; //Because we are reversing, we cant keep track of all the paths which were already visited in this walk
		path_add_node(&new_step, destination);
        // Update step flags - could be more elegant
        //destination->step_flags[destination->length-1] = source->step_flags[i];
	}
	
}

int path_get_index_of_last_in_node(Path * p){
    int count = p->in_nodes_count;
    int index = 0;
    if (count) {
        index = p->in_nodes[count-1];
    }
    return index;

}

boolean path_append(Path * destination, Path * source){
	int i;
	pathStep new_step;
	pathStep first_step;
	pathStep last_step;
	boolean success = true;

	if (source->length == 0) {
		fprintf(stderr, "[path_append] The source path is empty!\n");
		//exit(-1);
        assert(source->length != 0);
		return false;
	}

	if ((destination->length > 0)
	    &&
	    (unlabelled_path_step_equals
	     (path_get_step_at_index(0, &first_step, source),
	      path_get_last_step(&last_step, destination)))) {
		path_remove_last(destination);
		if (DEBUG) {
			printf("[path_append] Removing last step.\n");
		}
	} else {
		if (DEBUG) {
			printf("[path_append] No need to remove last step.\n");
		}
	}

	for (i = 0; i < source->length; i++) {
		new_step.node = source->nodes[i];
		new_step.orientation = source->orientations[i];
		new_step.label = source->labels[i];
        new_step.flags = source->step_flags[i];
		if (!path_add_node(&new_step, destination)) {
			success = false;
			break;
		} /*else {
            // Update step flags - could be more elegant
            destination->step_flags[destination->length - 1] = source->step_flags[i];
        }*/
	}

	return success;
}

int path_get_nodes_count(Path * path){
	return path->length;
}

int path_get_edges_count(Path * path){
	return path->length - 1;
}

boolean paths_equal(Path * path_a, Path * path_b){
	int i;
	boolean paths_equal = true;

	if ((!path_a) || (!path_b))
		return false;

	if (path_a->length != path_b->length)
		return false;

	for (i = 0; i < path_a->length; i++) {
		if ((path_a->nodes[i] != path_b->nodes[i]) ||
		    (path_a->labels[i] != path_b->labels[i]) ||
		    (path_a->orientations[i] != path_b->orientations[i])) {
			paths_equal = false;
			break;
		}
	}

	return paths_equal;
}

static PathArray *path_buffers = NULL;
//static PathArray *short_path_buffers = NULL;

/**
 * Makes a copy of the "from" path into the "to" path. 
 * WARNING! The to path is cleared before the copy. 
 */
void path_copy(Path * to, Path * from)
{
	path_reset(to);
    to->max_virtual_length = from->max_virtual_length;
	to->flags = from->flags;
	to->stop_reasons_first = from->stop_reasons_first;
	to->stop_reasons_last = from->stop_reasons_last;
    if (from->length > 0) {
       path_append(to, from);
    }
	
}

#ifdef THREADS
static pthread_once_t path_array_init_once = PTHREAD_ONCE_INIT;
#endif
static void path_buffers_init(){
    path_buffers = path_array_new(MAX_PATH_BUFFERS);
}

void path_array_initialise_buffers(short kmer_size)
{
#ifdef THREADS	
    pthread_once(&path_array_init_once, path_buffers_init);
#else
    if (path_buffers == NULL) {
        path_buffers_init();
    }
#endif
	if (path_buffers->kmer_size == 0) {
		path_buffers->kmer_size = kmer_size;
        
//		path_buffers_short = path_array_new(MAX_PATH_BUFFERS);
//		for (i = 0; i < MAX_PATH_BUFFERS; i++) {
//			tmp = path_new(MAX_PATH_LENGTH, kmer_size);
//			tmp->id = i;
//			path_array_add_path(tmp, path_buffers);
			
//		}
	}
}

void path_array_destroy_buffers(){
	Path *tmp;
	int i;
	if (path_buffers != NULL){
		for (i = 0; i < path_buffers->number_of_paths; i++) {
			tmp = path_array_get(i, path_buffers);
			path_destroy(tmp);
		}
		free(path_buffers);
		path_buffers = NULL;
	}
}

PathArray *path_array_get_from_buffer_with_size(short size)
{
	PathArray *pa = path_array_new(size);
	int i;
	Path * tmp = NULL;
	for (i = 0; i < size; i++) {
		tmp = path_get_buffer_path();
		path_array_add_path(tmp, pa);
		tmp->id = i;
	}
    if (tmp != NULL) {
        pa->kmer_size= tmp->kmer_size;

    }
	return pa;
}

Path * path_array_get(int path, PathArray *pa){
	assert(path < pa->number_of_paths);
	return pa->paths[path];
}

void path_array_free_from_buffer(PathArray * pa)
{
	while (pa->number_of_paths > 0) {
		pa->number_of_paths--;
		path_free_buffer_path(pa->paths[pa->number_of_paths]);

	}
	free(pa->paths);
	free(pa);

}

void path_array_to_fasta(FILE * f, PathArray * pa){
	int i;
	for (i = 0; i<pa->number_of_paths; i++) {
		if (path_get_length(path_array_get(i,pa)) > 0) {
			path_to_fasta(path_array_get(i, pa), f);
		}
		
	}
}

Path *path_get_buffer_path(){
	//TODO make this thread safe
	assert(path_buffers !=NULL);
    
#ifdef THREADS
    pthread_mutex_lock(&path_buffers->mutex);
#endif
    
	Path *tmp = NULL, *found = NULL;
	int i;
	for (i = 0; i < path_buffers->number_of_paths && found == NULL; i++) {
		tmp = path_buffers->paths[i];
		if (!tmp->used) {
			found = tmp;
			
		}
	}
    tmp = NULL;
    assert(i < MAX_PATH_BUFFERS); //TODO: make this a growing array. 
    if(found == NULL){
    
        tmp = path_new(MAX_PATH_LENGTH, path_buffers->kmer_size);
        tmp->id = i;
        path_array_add_path(tmp, path_buffers);
        found = path_buffers->paths[i];
        
        //This is the old logic, that was in the initializer. 
        //		path_buffers_short = path_array_new(MAX_PATH_BUFFERS);
        //		for (i = 0; i < MAX_PATH_BUFFERS; i++) {
        //			tmp = path_new(MAX_PATH_LENGTH, kmer_size);
        //			tmp->id = i;
        //			path_array_add_path(tmp, path_buffers);
        
        //		}
    }
    
	assert(found != NULL);
    found->used = true;
#ifdef THREADS
    pthread_mutex_unlock(&path_buffers->mutex);
#endif
	return found;
}

void path_free_buffer_path(Path * path)
{
#ifdef THREADS
    pthread_mutex_lock(&path_buffers->mutex);
#endif
	path_reset(path);
	path->used = false;
    
#ifdef THREADS
    pthread_mutex_unlock(&path_buffers->mutex);
#endif
}

void path_step_mark_as_uncertain(int i, Path * path, boolean as_n) {
    if (as_n) {
        path->step_flags[i] |= PRINT_LABEL_AS_N;
    } else {
        path->step_flags[i] |= PRINT_LABEL_LOWERCASE;
    }
}

boolean is_step_marked_as_uncertain(int i, Path * path) {
	boolean r = (path->step_flags[i] & (PRINT_LABEL_AS_N | PRINT_LABEL_LOWERCASE)) > 0; 

	return r;    
}

static void step_mark_visited(pathStep * ps) {
    db_node_action_set_flag_visited(ps->node);
}

void path_mark_as_visited(Path* path) {
	path_iterator(&step_mark_visited, path);
}

// TO DO: This needs to be done properley!
void path_pairs_to_fasta(PathArray* pa, int distances[], FILE* fout) {
    int i, j;
    int total_length = 0;
    int kmer_size;

	// Sanity checking
	if (pa == NULL) {
		fprintf(stderr,
                "[path_pairs_to_fasta] trying to print a null Path\n");
		exit(-1);
	}
    
    if (pa->number_of_paths < 2) {
		fprintf(stderr,
                "[path_pairs_to_fasta] trying to print less than one path\n");
		exit(-1);
    }
    
	if (fout == NULL) {
		fprintf(stderr,
                "[path_pairs_to_fasta] trying to print to a null FILE\n");
		exit(-1);
	}
    
    kmer_size = pa->paths[0]->kmer_size;
    
    for (i=0; i<pa->number_of_paths; i++) {
        total_length += strlen(pa->paths[i]->seq) + kmer_size;
        if (i < (pa->number_of_paths-1)) {
            total_length += distances[i];
        }
    }
    
    int current = 1;
    fprintf(fout, ">rpnode_%qd length:%i\n", pa->paths[0]->id, total_length);
    for (i=0; i<pa->number_of_paths; i++) {
        Path *path = pa->paths[i];
        dBNode* fst_node = path->nodes[0];
        BinaryKmer fst_kmer;
        char fst_seq[kmer_size+1];
        binary_kmer_assignment_operator(fst_kmer, *(element_get_kmer(fst_node)));
        binary_kmer_to_seq(&fst_kmer, kmer_size, fst_seq);
        
        // Print first kmer
        for (j=0; j<kmer_size; j++, current++) {
            fprintf(fout, "%c", fst_seq[j]);
            if(current % PATH_FASTA_LINE == 0){
                fprintf(fout, "\n");	
            }            
        }
        
        // Print rest
        for (j=0; j<strlen(path->seq); j++, current++) {
            if (path->step_flags[i] & PRINT_LABEL_AS_N) {
                fprintf(fout, "N");
            } else if (path->step_flags[i] & PRINT_LABEL_LOWERCASE) {
                fprintf(fout, "%c", tolower(path->seq[i]));
            } else {
                fprintf(fout, "%c",  path->seq[i]);
            }
            if(current % PATH_FASTA_LINE == 0){
                fprintf(fout, "\n");	
            }
        }
        fprintf(fout, "\n");
        
        // Print Ns
        if (i < (pa->number_of_paths-1)) {
            for (j=0; j<distances[i]; j++, current++) {
                if(current % PATH_FASTA_LINE == 0){
                    fprintf(fout, "\n");	
                }
                fprintf(fout, "N");                
            }
        }
    }
}


void path_counts_reset(PathCounts * pc){
	pc->blunt_ends = 0;
	pc->converging_paths = 0;
	pc->diverging_paths = 0;
	pc->is_double_y = 0;
	pc->is_cycle = 0;
	pc->longer_than_buffer = 0;
	 
	pc-> minimum_double_y = 0; 
	pc-> total_double_y_lenght = 0;
}

void path_counts_print_and_log(PathCounts * pc){
	log_and_screen_printf("Blunt Ends\t%'lld\n", pc->blunt_ends);
	log_and_screen_printf("Converging paths\t%'lld\n", pc->converging_paths);
	log_and_screen_printf("Diverging paths\t%'lld\n", pc->diverging_paths);
	log_and_screen_printf("Is cycle\t%'lld\n", pc->is_cycle);
	log_and_screen_printf("Is double y\t%'lld\n", pc->is_double_y);
	log_and_screen_printf("Longer tha buffer\t%'lld\n", pc->longer_than_buffer);
}

void path_counts_add(Path * p, PathCounts * pc){
	
	
	if(path_has_stop_reason(LAST, PATH_FLAG_STOP_BLUNT_END,  p)){
		pc->blunt_ends++;
	}
	if(path_has_stop_reason(LAST, PATH_FLAG_CONVERGING_PATHS, p)){
		pc->converging_paths++;
	}
	if(path_has_stop_reason(LAST, PATH_FLAG_DIVERGING_PATHS, p)){
		pc->diverging_paths++;
	}
	
	if(path_has_stop_reason(LAST, PATH_FLAG_IS_CYCLE, p)){
		pc->is_cycle++;
	}
	
	if(path_has_stop_reason(LAST, PATH_FLAG_IS_DOUBLE_Y, p)){
		pc->is_double_y++;
	}
	
	if(path_has_stop_reason(FIRST, PATH_FLAG_STOP_BLUNT_END,  p)){
		pc->blunt_ends++;
	}
	if(path_has_stop_reason(FIRST, PATH_FLAG_CONVERGING_PATHS, p)){
		pc->converging_paths++;
	}
	if(path_has_stop_reason(FIRST, PATH_FLAG_DIVERGING_PATHS, p)){
		pc->diverging_paths++;
	}
	
	if(path_has_stop_reason(FIRST, PATH_FLAG_IS_CYCLE, p)){
		pc->is_cycle++;
	}
	
	if(path_has_stop_reason(FIRST, PATH_FLAG_IS_DOUBLE_Y, p)){
		pc->is_double_y++;
	}
}

#ifdef SOLID

static NucleotideBaseSpace path_get_first_base_from_position(int pos, Path * p){
    int i = pos;
    NucleotideBaseSpace nuc = Undef; 
    Nucleotide trans = Undefined;
    pathStep ps;
    path_get_step_at_index(i, &ps, p);
    if(i == 0){
       return db_node_get_starting_base(ps.orientation, ps.node); //If this is the first kmer, we dont need to go back. 
    }
    
    nuc = db_node_get_starting_base(ps.orientation, ps.node);
    i--;
    for(; i >= 0 && nuc != Undef; i--){
        
        
        
        path_get_step_at_index(i, &ps, p);
        if(ps.orientation == forward){
            trans = binary_kmer_get_first_nucleotide(&(ps.node->kmer), p->kmer_size);
        }else{
            trans = binary_kmer_get_last_nucleotide(&(ps.node->kmer));
        }
        //printf("transforming: ");
        //path_step_print(&ps, p->kmer_size, stdout);
        //printf(" %c %c> ", binary_nucleotide_base_space_to_char(nuc), binary_nucleotide_to_char(trans));
        nuc = binary_nucleotide_base_space_get_next_base(nuc, trans);
        
        //printf("%c\n", binary_nucleotide_base_space_to_char(nuc));

    }
       
    return nuc;
}

NucleotideBaseSpace path_get_first_base(Path * p){
	int i; 
	pathStep ps;

    int count[4];
    int max = 0;
	NucleotideBaseSpace nuc = Undef; 
    
   // Nucleotide trans = Undefined;
	//path_get_step_at_index(0, &ps, p);
   // nuc = db_node_get_starting_base(ps.orientation, ps.node);
   // printf("sarched... first %c, of length %d\n",  binary_nucleotide_base_space_to_char(nuc), p->length);
    
    for(i= 0; i < 4; i++){
        count[i]=0;
    }
    
    for(i = 0 ; /*nuc != Undef && */ i < p->length ; i++){
        
        path_get_step_at_index(i, &ps, p);
	    //nuc = db_node_get_starting_base(opposite_orientation(ps.orientation), ps.node);
        nuc = db_node_get_starting_base(ps.orientation, ps.node);
        if(nuc != Undef){
          //  printf("\n searching..\n"); 
           // path_step_print(&ps, p->kmer_size, stdout);
            nuc = path_get_first_base_from_position(i, p);
          //  printf("sarched... %d %c\n",i,  binary_nucleotide_base_space_to_char(nuc));
            count[nuc]++;
        }
    }
   
    
    
    for(i= 0; i < 4; i++){
        //printf("%c:%d\n",  binary_nucleotide_base_space_to_char(i), count[i]);
        if(count[i] > max){
            max = count[i];
            nuc = i;
        }
    }
    
    if(max == 0){
        nuc = Undef;
    }
    //printf("sarched... found %c\n",  binary_nucleotide_base_space_to_char(nuc));
    
    
    return nuc;
        
}
//This only makes sense when the run is all SOLiD
void path_to_base_space_fasta(Path * path, FILE * fout)
{
	short kmer_size = path->kmer_size;
	int length = path->length;

    
	// Sanity checking
	assert(path != NULL);
    assert(fout != NULL);
	
    NucleotideBaseSpace nuc = path_get_first_base(path);
    
    if(nuc == Undef){
        path_to_fasta(path, fout);
        return;//We will at least get the colour sequence.
    }
    
	// Get coverage statistics from path
	double avg_coverage;
	int min_coverage;
	int max_coverage;
	path_get_statistics(&avg_coverage, &min_coverage, &max_coverage, path);

	// Get orientation of first and last node
	Orientation fst_orientation;

    fst_orientation = path->orientations[0];
	
	Orientation lst_orientation = path->orientations[path->length];

	
	dBNode *fst_node;
	
    if (path->length == 0) {
        fprintf(stderr, "[path_to_fasta] Trying to print an empty path[1]!\n");
        return;
    }
    fst_node = path->nodes[0];
	

	// Get the last node
	dBNode *lst_node = path->nodes[path->length - 1];

	// Make a set of labels for first and last nodes which list the
	// acceptable forward and reverse path labels
	char fst_f[5], fst_r[5], lst_f[5], lst_r[5];
	compute_label(fst_node, forward, fst_f);
	compute_label(fst_node, reverse, fst_r);
	compute_label(lst_node, forward, lst_f);
	compute_label(lst_node, reverse, lst_r);

	// Places to store first and last kmer sequence
	char fst_seq[kmer_size + 1], lst_seq[kmer_size + 1];
	fst_seq[kmer_size] = '\0';
	lst_seq[kmer_size] = '\0';
	BinaryKmer tmp_kmer;

	// Get the first kmer sequence (or complement)
	BinaryKmer fst_kmer;
	binary_kmer_assignment_operator(fst_kmer, *(element_get_kmer(fst_node)));
	if (fst_orientation == reverse) {
		binary_kmer_reverse_complement(&fst_kmer, kmer_size, &tmp_kmer);
		binary_kmer_assignment_operator(fst_kmer, tmp_kmer);
	}
	binary_kmer_to_seq(&fst_kmer, kmer_size, fst_seq);
	if (DEBUG) {
		printf("[path_to_fasta] First kmer: %s\n", fst_seq);
    }
    
	// Get the last kmer sequence (or complement)
	BinaryKmer lst_kmer;
	binary_kmer_assignment_operator(lst_kmer, *(element_get_kmer(lst_node)));
	if (lst_orientation == reverse) {
		binary_kmer_reverse_complement(&lst_kmer, kmer_size, &tmp_kmer);
		binary_kmer_assignment_operator(lst_kmer, tmp_kmer);
	}
	binary_kmer_to_seq(&lst_kmer, kmer_size, lst_seq);

	// Output to file
	fprintf(fout,
            ">translated_node_%qd length:%i average_coverage:%.2f min_coverage:%i max_coverage:%i fst_coverage:%i fst_r:%s fst_f:%s lst_coverage:%i lst_r:%s lst_f:%s\n",
            path->id,
            (flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? length + kmer_size : length + kmer_size - 1), avg_coverage,
            min_coverage,
            max_coverage,
            element_get_coverage_all_colours(fst_node),
            (fst_orientation == forward ? fst_r : fst_f),
            (fst_orientation == forward ? fst_f : fst_r),
            element_get_coverage_all_colours(lst_node),
            (lst_orientation == forward ? lst_r : lst_f),
            (lst_orientation == forward ? lst_f : lst_r));

	
	int i, current= 1;
    
    
    fprintf(fout, "%c", binary_nucleotide_base_space_to_char( nuc));
    
//binary_kmer_to_seq(&fst_kmer, kmer_size, fst_seq);

    nuc = binary_kmer_to_base_seq(&fst_kmer, forward, nuc, fst_seq, path->kmer_size);
    
	for(i = 0, current ; i < path->kmer_size; i++, current++) {
		fprintf(fout, "%c", fst_seq[i]);
		if(current % PATH_FASTA_LINE == 0) {
			fprintf(fout, "\n");	
		}
	}
	
	size_t len = strlen(path->seq);
	char nuc_c;
    Nucleotide trans;
    for(i = 0; i < len; i++, current++) {
        
        trans = path->labels[i];
        //printf("\n PATH: %c %c > ", binary_nucleotide_base_space_to_char(nuc), binary_nucleotide_to_char(trans));
        nuc = binary_nucleotide_base_space_get_next_base(nuc, trans);
        nuc_c = binary_nucleotide_base_space_to_char(nuc);
        
        //printf("%c", nuc_c);
        
        
        fprintf(fout, "%c", nuc_c);
		if(current % PATH_FASTA_LINE == 0) {
			fprintf(fout, "\n");	
		}
	}
	fprintf(fout, "\n");
	
	fflush(fout);
}
#endif

#ifdef ENABLE_MARK_PAIR

static void init_supernodes_array(Path * path){
    path->supernodes = calloc(PATH_DEFAULT_SUPERNODES, sizeof(uint32_t));
    path->supernodes_count = calloc(PATH_DEFAULT_SUPERNODES, sizeof(uint32_t));
    if (path->supernodes_count == NULL) {
        fprintf(stderr, "Unable to allocate memory for supernodes count\n");
        assert(0);
        exit(-1);
    }
    if (path->supernodes == NULL) {
        fprintf(stderr, "Unable to allocate memory for supernodes\n");
        assert(0);
        exit(-1);
    }
    path->supernodes_capacity = PATH_DEFAULT_SUPERNODES;
    path->last_supernode = -1;
    path->supernodes_total_count = 0;
}

static void grow_supernodes_array(Path * path){
    unsigned long new_size = path->supernodes_capacity + PATH_SUPERNODES_CAPACITY_INCREASE;
    uint32_t * tmp_arr = realloc(path->supernodes, new_size * sizeof(uint32_t));
    if (tmp_arr == NULL) {
        fprintf(stderr, "Unable to reallocate memory for supernodes \n");
        assert(0);
        exit(-1);
    }
    path->supernodes = tmp_arr;
    tmp_arr = realloc(path->supernodes_count, new_size * sizeof(uint32_t));
    if (tmp_arr == NULL) {
        fprintf(stderr, "Unable to reallocate memory for supernodes count\n");
        assert(0);
        exit(-1);
    }
    path->supernodes_count = tmp_arr;
}

static int find_index_of_supernode(uint32_t supernode, Path * path){
    int i;
    int index = -1;
    if (path->supernodes == 0) {
        return index;
    }
    if (path->last_supernode > path->supernodes_total_count) {
        return index;
    }
    if (path->supernodes[path->last_supernode] == supernode) {
        index = path->last_supernode;
    }
    
    for (i = 0; i < path->supernodes_total_count && index == -1; i++) {
        if (path->supernodes[i] == supernode) {
            index = i;
        }
    }
    return index;
}

void path_add_supernode(uint32_t id,  Path * path){
    if (path->supernodes_capacity == 0) {
        init_supernodes_array(path);
    }
    if (path->supernodes_total_count == path->supernodes_capacity) {
        grow_supernodes_array(path);
    }
    int index = find_index_of_supernode(id, path);
    if (index == -1) {
        index = path->supernodes_total_count;
        path->supernodes[index] = id;
        path->supernodes_count[index] = 0;//force initialization of the counter. 
        path->supernodes_total_count++;
    }
    path->last_supernode = index;
    path->supernodes_count[index]++;
}

void path_remove_supernode(uint32_t id, Path * path){
    int index = find_index_of_supernode(id, path);
    path->supernodes_count[index]--;
    while (path->supernodes_count[path->supernodes_total_count - 1] == 0) {
        path->supernodes_total_count--;
        index = path->supernodes_total_count;
        path->supernodes[index] = 0;
        path->last_supernode = index - 1;
    }
}


#endif
