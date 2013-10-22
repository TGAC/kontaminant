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

#ifdef ENABLE_READ_PAIR

#include <limits.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>

#include <global.h>
#include <flags.h>
#include <nucleotide.h>
#include <binary_kmer.h>
#include <seq.h>
#include <element.h>
#include <open_hash/hash_table.h>
#include <dB_graph.h>
#include <path.h>
#include <perfect_path.h>
#include <y_walk.h>
#include <branches.h>
#include <binary_tree.h>
#include <read_pair.h>
#include <logger.h>

// Function prototypes
int read_pair_find_mate(Path* start_path, Path* mate_path, ReadPairDescriptorArray* rpda, dBGraph* db_graph);
void print_labels_within_insert_size(Path* path, ReadPairDescriptorArray* rpda, ReadPairDescriptor* rpd);
char* print_binary(ReadPairSignature n, char* s);
ReadPairJump* read_pair_find_jump(ReadPairDescriptorArray* rpda, dBNode* node);

#ifdef READ_PAIR_DEBUG_GRAPH
void output_bitfield_data(Path* path, ReadPairDescriptorArray *rpda, ReadPairDescriptor* rpd);
#endif

ReadPairDescriptorArray * new_read_pair_descriptor_array(char capacity, int walk_distance, int max_coverage,
                                                         int min_bits, int start_length, int max_paths,
                                                         int min_kmers, boolean stack_as_n)
{
	ReadPairDescriptorArray * ret = calloc(1, sizeof(ReadPairDescriptorArray));
	ret->capacity = capacity;
	ret->number_of_pairs = 0;
    ret->walk_distance = walk_distance;
    ret->maximum_coverage = max_coverage;
    ret->minimum_bits = min_bits;
    ret->minimum_start_length = start_length;
    ret->maximum_paths = max_paths;
    ret->minimum_kmers = min_kmers;
    ret->print_stack_as_n = stack_as_n;
    ret->minimum_jump_insert_size = 2000; //TODO: make this a parameter. 
	ret->pair = calloc(capacity, sizeof(ReadPairDescriptor*));
    
    printf("Read pair array created:\n");
    printf("         Walk distance: %i\n", ret->walk_distance);
    printf("      Maximum coverage: %i\n", ret->maximum_coverage);
    printf("              Min bits: %i\n", ret->minimum_bits);
    printf("  Minimum start length: %i\n", ret->minimum_start_length);
    printf("  Maximum search paths: %i\n", ret->maximum_paths);
    printf("         Minimum kmers: %i\n", ret->minimum_kmers);
	
    return ret;
}

void destroy_read_pair_descriptor_array(ReadPairDescriptorArray ** rpda)
{
	free((*rpda)->pair);
	free(*rpda);
	(*rpda) = NULL;
}

void read_pair_descriptor_array_add_pair(char first, char second, int colour, int insert, int read_length, int tolerance, int misses, int  min_pair_coverage, int max_pair_coverage, int depth,  double (* score_paths)(Path * current, Path * extending, ReadPairDescriptor *rpda ),
                                         ReadPairDescriptorArray * rpda)
{
	assert(rpda->number_of_pairs < rpda->capacity);
	
	ReadPairDescriptor * rpd= calloc(1, sizeof(ReadPairDescriptor));
	
    if (insert <= read_length) {
        printf("Error: Insert size must be greater than read length\n");
        exit(-1);
    }
    
	//rpda->pair[rpda->number_of_pairs];
	rpd->first = first;
	rpd->second = second;
	rpd->colour = colour;
    rpd->allowed_misses = misses;
    rpd->insert_size = insert;
    rpd->max_read_length = read_length+1;
	rpd->tolerance = tolerance;
	rpd->score_paths = score_paths;
	rpd->distance = insert;
    rpd->min_pair_coverage = min_pair_coverage;
    rpd->max_pair_coverage = max_pair_coverage;
	rpd->mate = true; //TODO: make sure we use it and we know what we are doing!
	rpd->used_pairs = 0;
    rpd->insert_distance_sum = 0;
    rpd->depth = depth;
	//if(insert >= rpda->minimum_jump_insert_size){
	//	rpd->jumps = binary_tree_new(100000, &read_pair_jump_compare);//It is a big size, since we expect a considerable amount of jumps. We have to assess if this is good enough. Besides, the BT should grow dinamically. 
	//}else{
    rpd->jumps = NULL;
	//}
	rpda->pair[rpda->number_of_pairs] = rpd;
	rpda->number_of_pairs++;
	
    printf("Added new read pair descriptor:\n");
    printf("         Insert size: %i\n", insert);
    printf("         Read length: %i\n", read_length);
    printf("          Inset size: %i\n", insert);
    printf("           Tolerance: %i\n", tolerance);
    printf("      Allowed misses: %i\n", misses);
    
	//if (rpd->distance < rpda->minimum_distance) {
	//	rpda->minimum_distance = rpd->distance;
	//}
}

int read_pair_count_bits(ReadPairSignature n) {
    int count = 0;
    int number_of_bits = sizeof(READ_PAIR_DATATYPE) * 8;
    int i;
    
    for (i=0; i<number_of_bits; i++) {
        count += (n & 1);
        n >>= 1;
    }
    
    return count;
}

int read_pair_count_common_bits(ReadPairSignature left, ReadPairSignature right)
{
	//int count = __builtin_popcount (left & right); //Only works with GCC! http://gurmeetsingh.wordpress.com/2008/08/05/fast-bit-counting-routines/	
	//fprintf(stdout, "%x,%x,%d\n", left, right, count);
    int count = 0;
    int number_of_bits = sizeof(READ_PAIR_DATATYPE) * 8;
    int i;
    
    for (i=0; i<number_of_bits; i++) {
        if (left & right & 1) {
            count++;
        }
        
        left  >>= 1;
        right >>= 1;
    }
    
	return count;
}


double simple_iterated_score_paths(Path * left, Path * right, ReadPairDescriptor * rpd)
{
	int distance = rpd->distance;
	
    //	int right_index =  path_get_length(left) - 1;
	
	int l_length = path_get_length(left);
	int r_length = path_get_length(right);	
	int left_offset = path_get_length(left) - 2 - distance;
	if (left_offset < 0) {
		return  0.0; // Will reduce to the half the probability, if we cant check how many are congruent.
	}
	if(path_get_length(right) < 1){
		return 0.0; //To short to relay without doing a deep search
	}
	double total_score = 0.0 , tmp_score;
	
	int i_left = left_offset, i_right = 1, total = 0;
	pathStep step_left, step_right; 
	
	
	ReadPairSignature rps_l;
	ReadPairSignature rps_r;
	int left_cov, right_cov, max_cov, hits;
	//printf("%d %d %d %d %d %d\n", i_left, l_length, i_right, r_length, total, distance);
	while(i_left < l_length && i_right < r_length /*&& total < distance*/){//Tha last option is to stop searching once the 
		
		
		path_get_step_at_index(i_left, &step_left, left);
		path_get_step_at_index(i_right, &step_right, right);
        
		right_cov =element_get_coverage_by_colour(step_right.node, rpd->colour);
		left_cov  =element_get_coverage_by_colour(step_left.node, rpd->colour);
		//max_cov = right_cov > left_cov?right_cov:left_cov;
		max_cov = right_cov + left_cov;
		rps_r = db_node_get_signature(rpd->second, rpd->colour, step_right.node) | db_node_get_signature(rpd->first, rpd->colour, step_right.node);
		rps_l = db_node_get_signature(rpd->second, rpd->colour, step_left.node) | db_node_get_signature(rpd->first, rpd->colour, step_left.node); //For this simple algorithm, we merge both of the signatures. 
		
		hits = read_pair_count_common_bits(rps_l, rps_r);
		
		tmp_score = (double)hits/(double) max_cov;
		//printf("\ttmp_score: %f\n", tmp_score);
		total_score += tmp_score;
		
		i_right++;
		i_left++;
		total++;
		
	}
	//printf("\ttotal_score: %f\n", total_score);
	if(total > 0)
		total_score /= (double) total;
	//printf("\t\ttotal_score: %f\n", total_score);
	/*
     ReadPairSignature right_signature1 = db_node_get_signature(rpd->second, rpd->colour, step_right.node);
     ReadPairSignature left_signature1 = db_node_get_signature(rpd->first, rpd->colour, step_left.node);
     */
    //	int count1 = read_pair_count_common_bits(left_signature1, right_signature1);
	
	/*ReadPairSignature right_signature2 = right_signature1 | db_node_get_signature(rpd->first, rpd->colour, step_right.node);
     ReadPairSignature left_signature2  = left_signature1 | db_node_get_signature(rpd->second, rpd->colour, step_left.node);
     
     int count1 = read_pair_count_common_bits(left_signature1, right_signature1);
     int count2 = read_pair_count_common_bits(left_signature2, right_signature2);
     
     double score1 = (double)(count1)/(double)min_cov;
     double score2 = (double)(count2)/(double)min_cov;
     double max_score =  score1 > score2? score1:score2;
     //double max_score = score1;
     printf("cov %d %d %d %f\n", right_cov, left_cov, min_cov, max_score);
     return max_score > 0.10? max_score:-1;
     
     if (score1 > 0 && score2 > 0) {
     return score1 * score2;
     }else {
     return max_score;
     }
     */
	return total_score;
}

/**
 * Over simplistic alogirithm, just takes in account the first node of each path, and the node on the corresponding distance. 
 */ 
double simple_score_paths(Path * left, Path * right, ReadPairDescriptor * rpd)
{
	int distance = rpd->distance;
    //	int right_index =  path_get_length(left) - 1;
	
	int left_offset = path_get_length(left) - 2 - distance;
	if (left_offset < 0) {
		return  0.5; // Will reduce to the half the probability, if we cant check how many are congruent.
	}
	if(path_get_length(right) < 1){
		return 0.0; //To short to relay without doing a deep search
	}
	
	pathStep step_left, step_right; 
	path_get_step_at_index(left_offset, &step_left, left);
	path_get_step_at_index(1, &step_right, right);
	
	int right_cov =element_get_coverage_by_colour(step_right.node, rpd->colour);
	int left_cov  =element_get_coverage_by_colour(step_left.node, rpd->colour);
	int min_cov = right_cov > left_cov?right_cov:left_cov;
	
	
	
	ReadPairSignature right_signature1 = db_node_get_signature(rpd->second, rpd->colour, step_right.node);
	ReadPairSignature left_signature1 = db_node_get_signature(rpd->first, rpd->colour, step_left.node);
	
    //	int count1 = read_pair_count_common_bits(left_signature1, right_signature1);
	
	ReadPairSignature right_signature2 =  db_node_get_signature(rpd->first, rpd->colour, step_right.node);
	ReadPairSignature left_signature2  =  db_node_get_signature(rpd->second, rpd->colour, step_left.node);
	
	int count1 = read_pair_count_common_bits(left_signature1, right_signature1);
	int count2 = read_pair_count_common_bits(left_signature2, right_signature2);
	
	double score1 = (double)(count1)/(double)min_cov;
	double score2 = (double)(count2)/(double)min_cov;
	double max_score =  score1 > score2? score1:score2;
	//double max_score = score1;
	//printf("cov %d %d %d %f\n", right_cov, left_cov, min_cov, max_score);
	return max_score > 0.10? max_score:-1;
	
	if (score1 > 0 && score2 > 0) {
		return score1 * score2;
	}else {
		return max_score;
	}
}

/**
 * It gives an score between 0 and 1, where 1 is a perfect match in the flags
 * and 0 none of the flags are common. In principle, with this method we can 
 * multiply several scores obtained. It can be seen as a probability. If the 
 * sequence it is to short, the application should stop, so validate before sending
 * to this function that the length is enough!
 * This also asumes that the 
 * 
 */
double paired_signatures_score_paths(Path * left, Path * right, ReadPairDescriptor * rpd)
{	
	assert(left != right);
	assert(left->length >= rpd->distance);
	assert(right->length >= rpd->distance);
	
	pathStep left_step, right_step;
	path_get_last_step(&left_step, left);
	path_get_step_at_index(0, &right_step, right);
	
	//printf("left:\t");
	//path_step_print(&left_step, left->kmer_size, stdout);
	//printf("\nright:\t");
	//path_step_print(&right_step, left->kmer_size, stdout);
	//printf("\n");
	
	assert(left_step.node == right_step.node);
	assert(left_step.orientation == right_step.orientation);
	
	int index_left_1  = rpd->first;
	int index_right_1 = rpd->second;
	
	int index_left_2  = rpd->second;
	int index_right_2 = rpd->first;
	
	
	int tolerance = rpd->tolerance;
	int i, j;
	
	
	int current_index_left, current_index_right;
	int length_l = path_get_length(left);
	int length_r = path_get_length(right);
	int tmp_score_1;
	int tmp_score_2;
	int tmp_sig_l_1;
	int tmp_sig_r_1;
	int tmp_sig_l_2;
	int tmp_sig_r_2;
	int total_score = 0,
	score_1, score_2;
	int total_cov = 0;
    int tmp_cov_l, tmp_cov_r, tmp_tmp_cov;
	
	for(i = 0; i < rpd->distance; i++){
		tmp_cov_r = tmp_cov_l = tmp_tmp_cov = 0;
		current_index_right = rpd->distance -  i;
		score_1 = score_2 = 0;
		
		for(j = 1 - tolerance; j < tolerance; j++){
			current_index_left =  length_l  - i + j;
			path_get_step_at_index(current_index_left, &left_step, left); 
			tmp_sig_l_1 = db_node_get_signature(index_left_1, rpd->colour, left_step.node);
			tmp_sig_l_2 = db_node_get_signature(index_left_2, rpd->colour, left_step.node);
			tmp_cov_l = element_get_coverage_by_colour(right_step.node, rpd->colour);
			if (current_index_right >= 0 && current_index_right < length_r) {
				path_get_step_at_index(current_index_right, &right_step, left);				
				tmp_sig_r_1 = db_node_get_signature(index_right_1, rpd->colour, right_step.node);
				tmp_sig_r_2 = db_node_get_signature(index_right_2, rpd->colour, right_step.node);
				tmp_score_1 = read_pair_count_common_bits(tmp_sig_l_1, tmp_sig_r_1);
				tmp_score_2 = read_pair_count_common_bits(tmp_sig_l_2, tmp_sig_r_2);
				
				if(tmp_score_1 > score_1){
					score_1 = tmp_score_1;
					
					
					tmp_tmp_cov = element_get_coverage_by_colour(right_step.node, rpd->colour);
					if(tmp_tmp_cov > tmp_cov_r){
						tmp_cov_r = tmp_tmp_cov;
					}
				}
				
				if(tmp_score_2 > score_2){
					score_2 = tmp_score_2;
					
					tmp_tmp_cov = element_get_coverage_by_colour(right_step.node, rpd->colour);
					if(tmp_tmp_cov > tmp_cov_r){
						tmp_cov_r = tmp_tmp_cov;
					}
				}
				
			}
		}
		
		total_score += score_1 + score_2;
		total_cov += tmp_cov_l + tmp_cov_r;
	}
	
	fprintf(stderr, "Comparing paths (total score: %d, total cov: %d\n", total_score, total_cov);
	//path_to_fasta(left, stderr);
	//path_to_fasta(right, stderr);
	
	return (double) total_score / (double) total_cov;
} 

Path * read_pair_select_path( Path * current,  Orientation orientation, PathArray * possible_paths, ReadPairDescriptorArray *rpda)
{
	assert(orientation == forward || orientation == reverse);
	
	//if(current->length < read_pair_get_minimum_length(rpda))//We don't have a "long" path to start with... 
	//		return current;
	
	Path * tmp = path_get_buffer_path();
	
	if(orientation == reverse){
		path_reverse(current, tmp);
	}else{
		path_copy(tmp, current);
	}
	
	Path * left = 0;
	Path * right = 0;
	int i,j;	
	
	//double * scores = calloc(rpda->capacity, sizeof(double));
	
	if(orientation == reverse){
		left = tmp;
		// right is set below
	}else {
		right = tmp;
		// left is set below
	}
	
	double max_score=  0;
	int index_max = -1;
	double tmp_score;
	ReadPairDescriptor * rpd;
	for(i = 0; i < possible_paths->number_of_paths; i++){
		if(orientation == reverse){
			right = path_array_get(i, possible_paths);
		}else {
			left = path_array_get(i, possible_paths);
		}
		
		
		
		if (path_get_length(left) > 0 && path_get_length(right) > 0) {
			//printf("left\n");
			//path_to_fasta(left, stdout);
			//printf("right\n");
			//path_to_fasta(right, stdout);
			tmp_score = 1;
			for(j=0; j<rpda->number_of_pairs;j++){
				
				rpd =  rpda->pair[j];
				tmp_score *= rpd->score_paths(left, right, rpd);//TODO: Get the score and select the fittest path
				
				//tmp_score *= paired_signatures_score_paths(left, right,rpd);
				
			}
			//printf("best: %f tmp: %f\n", max_score, tmp_score);
			if (tmp_score > max_score) {//TODO: do an analysis of thresholds, to decide how to set it... and define a field in the rpda to hold it
				index_max = i;
				//printf("___________________Best path: %d (score:%f)\n", i, tmp_score);
				max_score = tmp_score;
			}
		}
		
		
	}
	
	if (index_max >= 0) {
		if(orientation == reverse){
			left = tmp;
			right = path_array_get(index_max, possible_paths);
		}else {
			right = tmp;
			left = path_array_get(index_max, possible_paths);
		}
		path_append(left,right);
		
		if(orientation == forward){
			path_reset(left);
			path_copy(left, right);
		}
	}else {
		if(orientation == reverse){
			path_copy(tmp, current);
			path_reset(current);
			path_reverse(tmp,current);
		}
	}
    
	if (path_get_length(left) > path_get_length(current)) {
		path_reset(current);
		path_copy(current, left);//TODO: Reduce the number of path inversions required...
        
	}
	if (path_get_length(right) > path_get_length(current)) {
		path_reset(current);
		path_copy(current, right);//TODO: Reduce the number of path inversions required...
		
	}
	
	//printf("final\n");
	
    
	//path_to_fasta(current, stdout);
	path_free_buffer_path(tmp);
	
	return current;
}

// Keeps growing path in current orientation until it can grow no more
void read_pair_grow_path(Path* path, ReadPairDescriptorArray* rpda, dBGraph* db_graph)
{
    Path *mate = path_get_buffer_path(); 
    boolean continue_growing;
    
    do {
        // Keep track of if we extend the path or not
        continue_growing = false;
        
        read_pair_find_mate(path, mate, rpda, db_graph);   
		
        if (mate->length == 0) {
#ifdef DEBUG_READ_PAIR
            printf("  No mate found\n");
#endif
        } else {
            // DEBUG
#ifdef DEBUG_READ_PAIR
            printf("  Adding mate of %i nodes\n", mate->length);
            printf("Got mate: ");
            path_to_fasta_debug(mate, stdout);
#endif
            
            // Add mate to path
            // We may need to validate VISITED_ flags
            assert(path->nodes[path->length-1] == mate->nodes[0]);
            assert(path->orientations[path->length-1] == mate->orientations[0]);
            continue_growing = path_append(path, mate);
            path_reset(mate);
            
            // DEBUG
#ifdef DEBUG_READ_PAIR
            printf("Merged path is: ");
            path_to_fasta_debug(path, stdout);
#endif            
        }
    } while (continue_growing);
    
	path_free_buffer_path(mate);    
}

void path_set_visited_forward_or_reverse(Path *path)
{
    void step_action_set(pathStep * ps) {
        db_node_action_set_flag(ps->node, ps->orientation == forward ? VISITED_FORWARD:VISITED_REVERSE);
    }
    path_iterator(&step_action_set, path);
}

void path_remove_visited_forward_or_reverse(Path *path)
{
    void step_action_unset(pathStep * ps) {
        db_node_action_unset_flag(ps->node, VISITED_FORWARD | VISITED_REVERSE);
	}
	path_iterator(&step_action_unset, path);    
}




//return true if path found

boolean read_pair_get_path(dBNode * node, void (*node_action) (dBNode * node), ReadPairDescriptorArray * rpda, dBGraph * db_graph, Path * path)
{	
    Path *reversed = path_get_buffer_path();
    char seq[1024];
    boolean ret = false;
	
    binary_kmer_to_seq(&node->kmer, db_graph->kmer_size, seq);
    
#ifdef DEBUG_READ_PAIR
    printf("\nTrying to get path from node %s\n", seq);
#endif
    
    // Count visited nodes
	int count_new=0;
	void count_new_nodes(dBNode* n) {
		if (db_node_check_flag_visited(n) == false) {
			count_new++;
		}
	}
    
    // Do a Y walk from this node - that gives us a path with a Y node at one or both ends
    // First thing y_walk does is reset path
	path = y_walk_get_path(node, reverse, &count_new_nodes, db_graph, true, path); //Getting a starting point    
    
	//we only take the path if it isn't repetitive 
	double avg_coverage;
	int min_coverage;
	int max_coverage;
	path_get_statistics(&avg_coverage, &min_coverage, &max_coverage, path);
    
#ifdef DEBUG_READ_PAIR
	printf(" New nodes %i total %i avg coverage %f length %i min length %i\n", count_new, path->length, avg_coverage, path->length,rpda->minimum_start_length);	
#endif  
    
    // Mark path as visited
    if (path) {
        path_mark_as_visited(path);
    }
    
    // Only bother to build path if...
    if ((path) && (path->length > 1) &&
        (count_new > rpda->minimum_kmers) &&
        (avg_coverage < rpda->maximum_coverage) &&
        (path->length > rpda->minimum_start_length))
    {			        
		ret = true;
        
        // Need to mark the nodes in this starting path as VISITED_FORWARD or VISITED_REVERSE
        // TODO: This must be done in the pre_step_action and cleaned in the post_step_action
        path_set_visited_forward_or_reverse(path);
        
        // DEBUG
#ifdef DEBUG_READ_PAIR
        printf("Got path: ");
        path_to_fasta_debug(path, stdout);
#endif
        
#ifdef DEBUG_PRINT_LABELS
#ifdef ENABLE_READ_PAIR_OLD
        //check consistency of labels within path
        printf("Entering print_labels_within_insert_size\n");
        print_labels_within_insert_size(path, rpda, rpda->pair[0]);
#endif
#endif
        
#ifdef READ_PAIR_DEBUG_GRAPH
        output_bitfield_data(path, rpda, rpda->pair[0]);
#endif
        
        // Grow in one direction
        read_pair_grow_path(path, rpda, db_graph);
        
        // Reverse and grow in the other
#ifdef DEBUG_READ_PAIR
        printf("  Trying other direction\n");
#endif
        path_reverse(path, reversed);
        read_pair_grow_path(reversed, rpda, db_graph);
        
        // Copy back into path for returning from function
        path_copy(path, reversed);
    }
    
    // Go through all of path, mark as visited and remove VISITED_FORWARD and VISITED_REVERSE
    void step_action_unset(pathStep * ps) {
        db_node_action_set_flag_visited(ps->node);
        db_node_action_unset_flag(ps->node, VISITED_FORWARD | VISITED_REVERSE);
	}
	path_iterator(&step_action_unset, path);
    
#ifdef DEBUG_READ_PAIR
    printf("  Marked as visited: %'lli of %'lli\n", get_visited_count(), db_graph->unique_kmers);
#endif
    
    // Free up buffers so they can be reused
    path_free_buffer_path(reversed);
    
	return ret; 
    
    
    
}





void read_pair_print_paths(char *filename, int max_length, boolean with_coverages, boolean with_viz, ReadPairDescriptorArray * rpda,dBGraph * db_graph)
{	
	Path *path = path_get_buffer_path();	//We will try to use only this buffer path.
	Path *next_path = path_get_buffer_path();
    
	path_reset(path);
	
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
	
	//perfect_path_get_funtions(&wf);
	//wf.continue_traversing = &continue_traversing;
	//wf.find_first_node = &always_false;
	
	void print_path(dBNode * node)
    {
		count_kmers++;
		if (db_node_check_flag_visited(node) == false) {
			db_node_action_set_flag_visited(node);
			
			boolean found_path = read_pair_get_path(node, &db_node_action_set_flag_visited, rpda, db_graph, path);
            
			if (found_path) {
                //pathStep last_step;
                //boolean found_next_path = false;
                
				if (path_get_nodes_count(path) > 1) {
					if (with_coverages) {
						path_to_coverage(path, fout_cov);
					}
                    
					if (with_viz) {
						path_graphviz_line(fout_viz, path);
					}
					
                    
                    path_to_fasta(path, fout);
                    
                    
                    if (path->length == max_length) {
                        printf("contig length equals max length [%i] for node_%i\n",
                               max_length, count_nodes);
                    }
                    
                    count_nodes++;
                    
					path_increase_id(path);
				} else {
					count_sing++;
				}
			}
		}
	}
	
	//find and mark double Ys 
	Path * buff = path_get_buffer_path();
	//db_graph_reset_flags(db_graph);
	mark_double_y(db_graph, buff);	
	
	log_and_screen_printf("Printing paths\n");
    
	hash_table_traverse(&print_path, db_graph);
	log_and_screen_printf("%'d nodes visited [%'qd singletons]\n", count_nodes, count_sing);
    
	//path_destroy(path);
	fclose(fout);
	if (with_coverages) {
		fclose(fout_cov);
	}
	if (with_viz) {
		path_graphviz_close_header(fout_viz);
		fclose(fout_viz);
	}
    
    path_free_buffer_path(path);
    path_free_buffer_path(buff);
    path_free_buffer_path(next_path);
    
	//	y_node_inited = false;	//Since we can clean the flags before running the node, we need to make sure next time we will find them	
}

void read_pair_enrich_graph(FILE * f1, FILE * f2, ReadPairDescriptor * rpd, int fastq_ascii_offset, dBGraph * db_graph)
{
	int max_read_length = rpd->max_read_length + 1; 
	Sequence * seq1 = sequence_new(max_read_length, LINE_MAX, 0);
	Sequence * seq2 = sequence_new(max_read_length, LINE_MAX, 0);
	int length_1, length_2, nkmers1, nkmers2;
    int missed_kmers_1, missed_kmers_2;
	short kmer_size = db_graph->kmer_size;
    
	boolean search_jumps = rpd->jumps == NULL? false: true; //Validate if the descriptor has space for the index of jumps. 
	
	//max_read_length/(kmer_size+1) is the worst case for the number of sliding windows, ie a kmer follow by a low-quality/bad base
	int max_windows = max_read_length / (kmer_size + 1);
	
	//number of possible kmers in a 'perfect' read
	int max_kmers = max_read_length - kmer_size + 1;
	
	KmerSlidingWindowSet * windows = binary_kmer_sliding_window_set_new(max_windows, max_kmers);
	windows->kmer_size = kmer_size;
	
	BinaryKmer tmp_key;
	int quality_cut_off = 0; //The cutoff comes from the previous CTXs, or the cleaning algorithm. 
	Path * p1 = path_new(max_read_length, kmer_size);
	Path * p2 = path_new(max_read_length, kmer_size);
	dBNode * current_node;
	pathStep ps;
	ps.label = Undefined;
	ps.orientation = forward;
	long long seq_count = 0;
	int ignore_count = 0;
	int used_count = 0;
	
	pathStep first_read;
	pathStep second_read;
	boolean first_blunt;
	boolean second_blunt;
	int first_length;
	int second_length;
	
	while((length_1 = read_sequence_from_fastq(f1, seq1, rpd->max_read_length) )&& (length_2 = read_sequence_from_fastq(f2, seq2, rpd->max_read_length))){
		seq_count++;
        
        if ((seq_count % 1000000) == 0) {
            printf("Enriched %lli sequences\n", seq_count);
        }
        
		path_reset(p1);
        missed_kmers_1 = 0;
        nkmers1 = get_sliding_windows_from_sequence(seq1->seq, seq1->qual, seq1->length, quality_cut_off, kmer_size, windows, max_windows, max_kmers,false,0);					
        
		while (binary_kmer_sliding_window_set_get_next(windows)) {
			element_get_key(&windows->current,kmer_size, &tmp_key); 		
			
			current_node = hash_table_find(&tmp_key, db_graph);
			if (current_node != NULL) {
				ps.node = current_node;
				
				if (binary_kmer_comparison_operator(tmp_key,windows->current)) {
					ps.orientation = forward;
				} else {
					ps.orientation = reverse;
				}
				path_add_node(&ps, p1);
			} else {
                missed_kmers_1++;
            }
		}
		
		path_reset(p2);
        missed_kmers_2 = 0;
        nkmers2 = get_sliding_windows_from_sequence(seq2->seq, seq2->qual, seq2->length, quality_cut_off, kmer_size, windows, max_windows, max_kmers, false, 0);
        
		while (binary_kmer_sliding_window_set_get_next(windows)) {
			element_get_key(&windows->current,kmer_size, &tmp_key); 
			current_node = hash_table_find(&tmp_key, db_graph);
			if(current_node != NULL) {
				ps.node = current_node;
				
				if (binary_kmer_comparison_operator(tmp_key, windows->current)) {
					ps.orientation = forward;
				} else {
					ps.orientation = reverse;
				}
				
				path_add_node(&ps, p2);
			} else {
                missed_kmers_2++;
            }
		}
        
        if ((missed_kmers_1 <= rpd->allowed_misses) && (missed_kmers_2 <= rpd->allowed_misses)) {
            used_count++;
        } else {
            ignore_count++;
        }
        
		if ((path_get_length(p1) > 0) && (path_get_length(p2) > 0) && (missed_kmers_1 <= rpd->allowed_misses) && (missed_kmers_2 <= rpd->allowed_misses)) { 
			// Validating that both reads are still in the graph. If they are not, it might mean that it makes no sense to enrich the pair. 
			
			long long label = seq_count + rand();
            
            
            
            short pair = 0; // RML: Update with pair number...
            
            void enrich_first_path(pathStep *ps) {
                if (ps->orientation == forward) {
                    db_node_action_add_read_pair(label, pair, Bitfield_First, ps->node);
                } else {
                    db_node_action_add_read_pair(label, pair, Bitfield_First_Inverse, ps->node);                    
                }
            }
            
            void enrich_second_path(pathStep *ps) {
                if (ps->orientation == forward) {
                    db_node_action_add_read_pair(label, pair, Bitfield_Second, ps->node);
                } else {
                    db_node_action_add_read_pair(label, pair, Bitfield_Second_Inverse, ps->node);                    
                }
            }
            
            path_iterator(&enrich_first_path, p1);
            
			path_iterator(&enrich_second_path, p2);
            
			if (search_jumps) {
				//TODO: This only works with illumina, we need to fix the mate pairs to work with other platforms. 
				// <---A----     ----B--->	(illumina mates). 
				// ---A---->     ----B--->	(illumina pairs). 
				
				
				//Leave the orientation of the second one,as it is, we will always look backards from the first
				path_get_step_reverse(&second_read, p2, 0);
				
				if (db_node_has_precisely_one_edge
					(second_read.node, second_read.orientation, &second_read.label) && db_node_check_flag_visited(second_read.node)) {
                    
					void is_second_blunt(Path * tmp_p){
						second_length = tmp_p->length;
						second_blunt = path_is_blunt(forward, tmp_p);
						//path_step_assign(&second_read,
                        path_get_last_step_reverse(&second_read, p2);
					}
                    
					perfect_path_get_path_from_step_with_callback(second_read, &db_node_action_set_flag_visited,&is_second_blunt ,db_graph);
					
					
				}
				
				
				//int perfect_path_get_path_from_step_with_callback(pathStep  first,
				//	void (*node_action) (dBNode * node),
				//	void (*path_action) (Path * path),
				//	dBGraph * db_graph)
				
				if(second_blunt){//We only want to do this if we already know that one half is blunt...
					if(rpd->mate){
						path_get_step_reverse(&first_read, p1, 0); 
					}else{
						path_get_last_step(&first_read, p1);
					}
					if (db_node_has_precisely_one_edge
                        (first_read.node, first_read.orientation, &first_read.label) && db_node_check_flag_visited(first_read.node)) {
						void is_first_blunt(Path * tmp_p){
							first_length = tmp_p->length;
							first_blunt = path_is_blunt(forward, tmp_p);
							//path_step_assign(&second_read,
							path_get_last_step_reverse(&first_read, tmp_p);
						}
						
						perfect_path_get_path_from_step_with_callback(first_read, &db_node_action_set_flag_visited,&is_first_blunt ,db_graph);
                        
					}
				}
				
				if(second_blunt && first_blunt ){
					read_pair_jump_add_to_descriptor(&first_read, &second_read, rpd->distance - first_length - second_length,  rpd);
				}		
				//read_pair_jump_add_to_descriptor(pathStep * a, pathStep * b, int distance, ReadPairDescriptor * rpd
			}
		} else {
            //RHRG: RML can bring it back when he fancies to debug       printf("p1: %d p2: %d mk1: %d mk2: %d\n", path_get_length(p1), path_get_length(p2), missed_kmers_1, missed_kmers_2);    
        }
	}
    
	path_destroy(p1);
	path_destroy(p2);
	binary_kmer_free_kmers_set(&windows);
	free_sequence(&seq1);
	free_sequence(&seq2);
    
    printf("Read %lli sequences\n", seq_count);
    printf("Ignored reads: %i Used reads: %i\n", ignore_count, used_count);
}

// ----------------------------------------------------------------------------------------------------
// To do:
// - Check all path_appends that we check for success.
// ----------------------------------------------------------------------------------------------------


// ----------------------------------------------------------------------------------------------------
// read_pair_get_next_step
//
// Given a 'stack' path, attempts to find a new edge and, hence, new node, to extend the path. The new
// edge must have a label that is alphabetically greater than that of the last step currently in the
// path. If a new edge is found, then the function will adjust the edge of the last step in the 'stack'
// and will return the next step. It is up to the caller to add the new step to the stack. Returns a
// boolean to indicate if the new step was found. This routine can be used in a depth-first search to walk 
// graph ensuring that every possible path is visited as the main depth-first routine backtracks. 
// ----------------------------------------------------------------------------------------------------
boolean read_pair_get_next_step(Path* stack, pathStep* next_step, dBGraph* db_graph)
{
    pathStep last_step;
    Nucleotide next_nucleotide;
    boolean found_candidate = false;
    
    // Stack must always have at least step before this function can be called
    assert(stack->length > 0);
    
    // Get last step
    path_get_last_step(&last_step, stack);
    next_nucleotide = last_step.label;
    
    // Keep trying nucleotides until we find one where an edge exists
    // AND that edge does not lead to a cycle. When we get back round
    // to Undefined, we stop.
    do {
        switch(next_nucleotide) { //What about using the nucleotide iterator? would that make it cleaner?
            case Undefined: next_nucleotide = Adenine; break;
            case Adenine: next_nucleotide = Cytosine; break;
            case Cytosine: next_nucleotide = Guanine; break;
            case Guanine: next_nucleotide = Thymine; break;
            case Thymine: next_nucleotide = Undefined; break;
        }
        
        // Check if an edge exists with this nucleotide
        if ((next_nucleotide != Undefined) &&
            (db_node_edge_exist_any_colour(last_step.node, next_nucleotide, last_step.orientation))) {
            pathStep step; 
            pathStep reverse_step;
            
            // Find the node this edge leads to - we make a copy of last_step, with the new label
            // because we might not end up keeping this label.
            step.node = last_step.node;
            step.label = next_nucleotide;
            step.orientation = last_step.orientation;
            db_graph_get_next_step(&step, next_step, &reverse_step, db_graph);
            
            // If the node hasn't already been visited, then it's a suitable candidate. TODO Richard we might consider to have this check outside this routine
            if (!db_node_check_for_any_flag(next_step->node, next_step->orientation == forward? VISITED_FORWARD:VISITED_REVERSE)) {
                found_candidate = true;
                
                // Now we need to adjust the label of the last step  -- TODO Richard great, this is potentially confusing but you did the right thing              
                path_remove_last(stack);
                last_step.label = next_nucleotide;
                path_add_node(&last_step, stack);
                
                // DEBUG
#ifdef DEBUG_READ_PAIR
                printf("Last step adjusted to %c\n", binary_nucleotide_to_char(last_step.label));
#endif
            }
        }
    } while ((next_nucleotide != Undefined) && (!found_candidate));
    
    return found_candidate;
}

// ----------------------------------------------------------------------------------------------------
// read_pair_score_pair
// 
// Score a pair of read pair signatures
// ----------------------------------------------------------------------------------------------------
int read_pair_score_pair(ReadPairSignature control, ReadPairSignature test, ReadPairDescriptorArray *rpda, int * score)
{
    int control_bits_set = read_pair_count_bits(control);
    int test_bits_set = read_pair_count_bits(test);
    int bits_in_common = read_pair_count_common_bits(control, test);
    int ret = 0;
    int lr_squared;
    int likely_random_bits;
    
    // Calculate likely random bits
    if (test_bits_set > control_bits_set) {
        lr_squared = (test_bits_set*test_bits_set);
    } else {
        lr_squared = (control_bits_set*control_bits_set);
    }
    
    // Round up likely random bits
    likely_random_bits = lr_squared / READ_PAIR_LENGTH;
    if ((lr_squared % READ_PAIR_LENGTH) != 0) {
        likely_random_bits++;
    }
    
    *score = bits_in_common;
    
#ifdef SIMPLE_SCORING
    // For VERY simple examples, we need to be stricter
    if ((bits_in_common > 0) && (control == test)) {
        ret = 1;
    }
#else
    if ((bits_in_common >= rpda->minimum_bits) &&
        (bits_in_common >= likely_random_bits) &&
        (bits_in_common >= ((control_bits_set / 2) - 1))) {
        ret = 1;
    }
#endif   
    return ret;
}

#ifndef READ_PAIR_OLD
int read_pair_score_nodes(pathStep* stack_step, pathStep* test_step, ReadPairDescriptorArray* rpda, int* common_bits_left, int* common_bits_right, int* left_index, int* right_index)
{
    int pair = 0;
    boolean score_left;
    boolean score_right;
    int ret = 0;
    int left_stack_bitfield = -1;
    int left_test_bitfield = - 1;
    int right_stack_bitfield = -1;
    int right_test_bitfield = -1;
    
    /* For READ PAIRS... */
    if ((stack_step->orientation == forward) && (test_step->orientation == forward)) {
        left_stack_bitfield = Bitfield_First;
        left_test_bitfield = Bitfield_Second_Inverse;
        right_stack_bitfield = Bitfield_Second_Inverse;
        right_test_bitfield = Bitfield_First;     
    } else if ((stack_step->orientation == forward) && (test_step->orientation == reverse)) {
        left_stack_bitfield = Bitfield_First;
        left_test_bitfield = Bitfield_Second;
        right_stack_bitfield = Bitfield_Second_Inverse;
        right_test_bitfield = Bitfield_First_Inverse;     
    } else if ((stack_step->orientation == reverse) && (test_step->orientation == forward)) {
        left_stack_bitfield = Bitfield_First_Inverse;
        left_test_bitfield = Bitfield_Second_Inverse;
        right_stack_bitfield = Bitfield_Second;
        right_test_bitfield = Bitfield_First;     
    } else if ((stack_step->orientation == reverse) && (test_step->orientation == reverse)) {
        left_stack_bitfield = Bitfield_First_Inverse;
        left_test_bitfield = Bitfield_Second;
        right_stack_bitfield = Bitfield_Second;
        right_test_bitfield = Bitfield_First_Inverse;     
    }
    
    if ((left_stack_bitfield == -1) || (left_test_bitfield == -1) || (right_stack_bitfield == -1) || (right_test_bitfield == -1)) {
        printf("Error: read pair bitfield not set. Can't continue.\n");
        exit(1);
    }
    
    /* Assume we're the left read */
    score_left = read_pair_score_pair(db_node_get_signature(pair, left_stack_bitfield, stack_step->node),
                                      db_node_get_signature(pair, left_test_bitfield, test_step->node),
                                      rpda,
                                      common_bits_left);
    /* Assume we're the right read */
    score_right = read_pair_score_pair(db_node_get_signature(pair, right_stack_bitfield, stack_step->node),
                                       db_node_get_signature(pair, right_test_bitfield, test_step->node),
                                       rpda,
                                       common_bits_right);
    
    /* Index - for debugging */
    *left_index = (left_stack_bitfield * 4) + left_test_bitfield;
    *right_index = (right_stack_bitfield * 4) + right_test_bitfield;
    
    /* Overall score... this is a good match if either score_left is very good, or
     score_right is very good, or both are ok, but added together they gain
     extra weight. */
    if ((score_left > 0) || (score_right > 0)) {
        ret = 1;
    }
    
    return ret;
}
#endif

// ----------------------------------------------------------------------------------------------------
// read_pair_check_for_match
// 
// See if we have found a node which is a mate to a node in our start path.
// ----------------------------------------------------------------------------------------------------
boolean read_pair_check_for_match(Path* start_path, Path* stack, ReadPairDescriptorArray* rpda, int * score, dBGraph* db_graph)
{
    int rp_in_range = 0;
    int rp_matches = 0;
    int i;
    int j;
    
	*score = 0;
    
    for (i=0; i<rpda->number_of_pairs; i++) {
        ReadPairDescriptor* rpd = rpda->pair[i];
        int walk_size = rpd->tolerance;
        int pair_distance = rpd->distance;
        int index_into_start_path = start_path->length - 1 - (pair_distance - stack->length + 1);
        pathStep stack_step;
        pathStep start_step;
        
        // Check if this pair is within range
        if ((index_into_start_path >= 0) && (index_into_start_path < start_path->length)) {
            int best_score = -1;
            
            rp_in_range++;
            
            for (j=index_into_start_path-walk_size; j<=index_into_start_path+walk_size; j++) {
                if (j >= 0 && j < start_path->length) {
                    path_get_last_step(&stack_step, stack);
                    path_get_step_at_index(j, &start_step, start_path);
                    
                    // DEBUG
                    //#ifdef DEBUG_READ_PAIR
                    //char seq1[1024];
                    //char seq2[1024];
                    //binary_kmer_to_seq(&(stack_step.node->kmer), db_graph->kmer_size, seq1);
                    //binary_kmer_to_seq(&(start_step.node->kmer), db_graph->kmer_size, seq2);
                    //printf("Comparing %s with %s (distance %d) index %i\n", seq1, seq2, pair_distance, j);
                    //#endif
                    
                    int this_score;
                    // Assuming the step we're comparing with isn't marked as an 'N' (ie. a messy bit!),
                    // then compare signatures...
#ifdef ENABLE_READ_PAIR_OLD
                    if ((!is_step_marked_as_uncertain(j, start_path)) &&
                        (read_pair_score_pair(stack_step.node->signature[0][0], start_step.node->signature[0][0], rpda, &this_score) == 1)) {
#else
                        int common_left, common_right, left_index, right_index;
                        if ((!is_step_marked_as_uncertain(j, start_path)) &&
                            (read_pair_score_nodes(&stack_step, &start_step, rpda, &common_left, &common_right, &left_index, &right_index) == 1)) {
                            this_score = common_left + common_right;
#endif
                            // Only increment rp_matches the first time we find a scoring match
                            if (best_score == -1) {
                                rp_matches++;
                            }
                            
                            // Keep track of best score
                            if (this_score > best_score) {
                                best_score = this_score;
                            }
                            
#ifdef DEBUG_READ_PAIR
                            printf("Found match %i on pair %i!!!\n", rp_matches, i);
#endif
                        }
#ifdef ENABLE_READ_PAIR_OLD  
                    }
#else
                }
#endif
            }
            
            if (best_score != -1) {
                *score += best_score;
            }            
        }
    }
    
    return ((rp_matches == rp_in_range) && (rp_matches > 0)) ? true:false;
}

// ----------------------------------------------------------------------------------------------------
// read_pair_find_mate
//
// Given a path, do a depth-first search to try to locate a read pair mate path which can be used to
// extend the path. In general we assume the last node in the path is a Y node and we need to choose the best path
// out of here (though this is not a necessary condition)
// ----------------------------------------------------------------------------------------------------
int read_pair_find_mate(Path* start_path, Path* mate_path, ReadPairDescriptorArray* rpda, dBGraph* db_graph)
{
    Path* stack = path_get_buffer_path();      //keeps track of the depth-first stack
    Path* best_path = path_get_buffer_path();
    pathStep last_step;
    pathStep next_step;
    boolean extended_path = false;
    boolean found_match = false;
    int max_match_paths = rpda->maximum_paths;
    int n_paths = 0;
    int best_score = 0;
    boolean rc;
    
    // Get last step from start path and push onto stack to start the search. 
    // Last step label gets set, so we have to reset to undefined
    // Richard, I think this is right, we are sometimes a bit confused about a path.
    // A path is a sequence of steps, but n steps result in a path of length of n+1 (right?) TODO	
    path_get_last_step(&last_step, start_path);
    last_step.label = Undefined;
    rc = path_add_node(&last_step, stack);    
    assert(rc == true);
    
    
    // Mark as visited
    db_node_action_set_flag(last_step.node, last_step.orientation == forward ? VISITED_FORWARD:VISITED_REVERSE); //these flags are used to detect cycles
    
    // DEBUG
#ifdef DEBUG_READ_PAIR
    char seq[1024];
    binary_kmer_to_seq(&(last_step.node->kmer), db_graph->kmer_size, seq);
    printf("Finding mate from last step kmer: %s label %c Orientation %s\n", seq, binary_nucleotide_to_char(last_step.label), last_step.orientation == forward ? "Fwd":"Rev");    
#endif
    
    int score = 0;
    do {
        do {
            // Reset flags and get last step from stack
            extended_path = false;
            found_match = false;
            path_get_last_step(&last_step, stack);
            
            // Try and find next label
            if (read_pair_get_next_step(stack, &next_step, db_graph)) {
                // Found a new node - push it onto the stack
                path_add_node(&next_step, stack); //TODO path_add_node should be called path_add_step
                extended_path = true;
                
                // DEBUG
#ifdef DEBUG_READ_PAIR
                printf("Push: ");
                path_to_fasta_debug(stack, stdout);
#endif
                
                // Mark new node as visited - so we can tell if we're in a cycle               
                db_node_action_set_flag(next_step.node, next_step.orientation == forward ? VISITED_FORWARD:VISITED_REVERSE);
                
                // Check if we've found what we're looking for				
                found_match = read_pair_check_for_match(start_path, stack, rpda, &score, db_graph);
            }
        } while ((extended_path) && (!found_match) && (stack->length <= rpda->walk_distance));
        
        // If we found a match, store it
        if (found_match) {
            n_paths++;
            if (score > best_score) {
                path_reset(best_path);
                path_copy(best_path, stack);
                best_score = score;
            }
        }
        
        // Pop top and unmark visited flag
        path_remove_last(stack);
        db_node_action_unset_flag(last_step.node, last_step.orientation == forward ? VISITED_FORWARD : VISITED_REVERSE);
        
        // DEBUG
#ifdef DEBUG_READ_PAIR
        if (stack->length > 0) {
            printf("Pop: ");
            path_to_fasta_debug(stack, stdout);
        } else {
            printf("Pop: end of path\n");
        }
#endif
    } while ((stack->length > 0) && (n_paths < max_match_paths));
    
    // There might be something on the stack still - if so, remove VISITED_FORWARD and VISITED_REVERSE flags 
    path_remove_visited_forward_or_reverse(stack);
    
    // Clear stack and point at best path
    path_free_buffer_path(stack);
    
    //printf("  Search revealed %i paths with best score of %i\n", n_paths, best_score);
    
    // We need to mark our best path as visited
    path_set_visited_forward_or_reverse(best_path);
    
    // Now find the mate path from where we are (next_step still contains the matching step)
    if ((n_paths > 0) && (best_score > 0)) {        
        int i=0;
        
        // Need to change next_step to be last step in best_path
        path_get_last_step(&next_step, best_path);
        
        // Before finding the mate path, we're going to mark the nodes between the start of the best stack path and the
        // mate point, so that labels get printed as Ns when output. The mate point is actually just the last node.        
        for (i=0; i<best_path->length-1; i++) {
            path_step_mark_as_uncertain(i, best_path, rpda->print_stack_as_n);
        }
        
        // Get mate path - this goes only forward now (see the new boolean variable)
        //mate_path = y_walk_get_path(next_step.node, next_step.orientation, &db_node_action_do_nothing, db_graph, false, mate_path);
        perfect_path_get_path(next_step.node, next_step.orientation, &db_node_action_do_nothing, db_graph, mate_path);        
        
        if (mate_path->length == 0) {
#ifdef DEBUG_READ_PAIR
            printf("Unable to get mate path\n");
#endif
        } else {        
            // DEBUG
#ifdef DEBUG_READ_PAIR
            printf("Start mate path: ");
            path_to_fasta_debug(mate_path, stdout);
#endif
            
            // Set VISITED_FORWARD or VISITED_REVERSE
            path_set_visited_forward_or_reverse(mate_path);
            
            // Join to best_path
            assert(best_path->nodes[best_path->length-1] == mate_path->nodes[0]);
            assert(best_path->orientations[best_path->length-1] == mate_path->orientations[0]);
            path_append(best_path, mate_path);
            
            // DEBUG
#ifdef DEBUG_READ_PAIR
            printf("best_path joined path: ");
            path_to_fasta_debug(best_path, stdout);
#endif        
            
            // And copy back to mate_path
            path_copy(mate_path, best_path);
        }
    }
    
    path_free_buffer_path(best_path);
    
    return best_score;
}

char* print_binary(ReadPairSignature n, char* s)
{
    int number_of_bits = sizeof(ReadPairSignature)*8;
    int j = number_of_bits;
    int i;
    
    s[j--] = 0;
    
    for (i=0; i<number_of_bits; i++) 
    {
        if (n & 1)
            s[j--] = '1';
        else
            s[j--] = '0';
        n >>= 1;    
    }
    
    return s;
}

#ifdef ENABLE_READ_PAIR_OLD
void print_labels_within_insert_size(Path* path, ReadPairDescriptorArray *rpda, ReadPairDescriptor* rpd)
{
    void step_action(int index,pathStep* ps)
    {
        pathStep paired_step;
        pathStep off_step;
        int paired_index = index+rpd->distance;
        int off_index    = index+rpd->distance+4*rpd->max_read_length;	
        char b1[256], b2[256], b3[256];
        
        if (index>10 && paired_index + 10 < path->length && off_index < path->length) {
            
            path_get_step_at_index(off_index,&off_step,path);
            
            printf("\nTEST (sizeof %li)\n", sizeof(READ_PAIR_DATATYPE));
            
            printf(" Index:     %8i %s %3i %16llX\n",
                   index,
                   print_binary(ps->node->signature[0][0], b1),
                   ps->node->coverage[0],
                   (long long unsigned int)ps->node->signature[0][0]);
            
            int j;
            for (j=-10; j<=10; j++) {
                int     npaired_index = paired_index+j;
                path_get_step_at_index(npaired_index,&paired_step,path);
                int score;
                
                printf("PIndex: %3i %8i %s %3i %16llX %3i %1i\n",
                       j,
                       npaired_index,
                       print_binary(paired_step.node->signature[0][0], b2),
                       paired_step.node->coverage[0],
                       (long long unsigned int)paired_step.node->signature[0][0],
                       read_pair_count_common_bits(ps->node->signature[0][0], paired_step.node->signature[0][0]),
                       read_pair_score_pair(ps->node->signature[0][0], paired_step.node->signature[0][0], rpda, &score));
            }                        
            
            path_get_step_at_index(paired_index,&paired_step,path);
            int score;
            printf("Paired:     %8i %s %3i %16llX %3i %1i\n",
                   paired_index,
                   print_binary(paired_step.node->signature[0][0], b2),
                   paired_step.node->coverage[0],
                   (long long unsigned int)paired_step.node->signature[0][0],
                   read_pair_count_common_bits(ps->node->signature[0][0], paired_step.node->signature[0][0]),
                   read_pair_score_pair(ps->node->signature[0][0], paired_step.node->signature[0][0], rpda, &score));
            
            printf("Offset:     %8i %s %3i %16llX %3i %1i\n",
                   off_index,
                   print_binary(off_step.node->signature[0][0], b3),
                   off_step.node->coverage[0],
                   (long long unsigned int)off_step.node->signature[0][0],
                   read_pair_count_common_bits(ps->node->signature[0][0],off_step.node->signature[0][0]),				
                   read_pair_score_pair(ps->node->signature[0][0], off_step.node->signature[0][0], rpda, &score));
        }
    }
    
    printf("Path length: %i\n",path->length);
    path_iterator_with_index(&step_action,path);
}
#endif

char* get_bitfield_name(ReadPairBitfield b, char* s)
{
    switch(b) {
        case Bitfield_First: strcpy(s, "First"); break;
        case Bitfield_Second: strcpy(s, "Second"); break;
        case Bitfield_First_Inverse: strcpy(s, "FirstInverse"); break;
        case Bitfield_Second_Inverse: strcpy(s, "SecondInverse"); break;
        default: strcpy(s, "Unknown"); break;
    }
    
    return s;
}

#ifdef READ_PAIR_DEBUG_GRAPH
#ifdef ENABLE_READ_PAIR_OLD
void output_bitfield_data(Path* path, ReadPairDescriptorArray *rpda, ReadPairDescriptor* rpd)
{
    static int contig_count = 0;
    pathStep middleStep;
    int middle_of_path = path->length / 2;
    int start = middle_of_path - 800;
    int end = middle_of_path + 800;
    FILE* bitfield_fp = 0;
    char filename[1024];
    
    if (start < 0) start = 0;
    if (end > path->length-1) end = path->length-1;    
    
    sprintf(filename, "contig%05d.txt", contig_count++);
    bitfield_fp = fopen(filename, "w");
    if (!bitfield_fp) {
        printf("Can't open similarity file.\n");
        exit(1);
    }     
    
    fputs("StepIndex\tScore\tCoverage\tBitsInCommon\n", bitfield_fp);
    
    path_get_step_at_index(path->length / 2, &middleStep, path);    
    
    void step_action(int index, pathStep* ps)
    {
        int cov = element_get_coverage_all_colours(ps->node);
        int match = 0;
        int common = 0; 
        int c;
        match = read_pair_score_pair(middleStep.node->signature[0][0], ps->node->signature[0][0], rpda, &c);
        common = read_pair_count_common_bits(ps->node->signature[0][0], middleStep.node->signature[0][0]);
        
        fprintf(bitfield_fp, "%i\t%i\t%i\t%i\n", index, match, cov, common);
    }
    
    path_iterator_with_index(&step_action, path);
    
    fclose(bitfield_fp);
}

#else    

void output_bitfield_data(Path* path, ReadPairDescriptorArray *rpda, ReadPairDescriptor* rpd)
{
    static int contig_count = 0;
    pathStep middleStep;
    int middle_of_path = path->length / 2;
    int start = middle_of_path - 800;
    int end = middle_of_path + 800;
    char filename[1024];
    FILE* bitfield_fp = 0;
    int i, j;
    char temp[1024];
    
    if (start < 0) start = 0;
    if (end > path->length-1) end = path->length-1;
    
    sprintf(filename, "contig%05d.txt", contig_count++);    
    bitfield_fp = fopen(filename, "w");
    if (!bitfield_fp) {
        printf("Can't open similarity file.\n");
        exit(1);
    } 
    
    fputs("StepIndex\tScore\tCoverage\tBitsInCommonLeft\tBitsInCommonRight\tBitsSetFirst\tBitsSetSecond\tBitsSetFirstInverse\tBitsSetSecondInverse", bitfield_fp);
    for (i=0; i<4; i++) {
        for (j=0; j<4; j++) {
            fputs("\tCommonBits", bitfield_fp);
            fputs(get_bitfield_name(i, temp), bitfield_fp);
            fputs("And", bitfield_fp);
            fputs(get_bitfield_name(j, temp), bitfield_fp);
        }
    }
    fputs("\tLeftIndex\tRightIndex\n", bitfield_fp);
    
    path_get_step_at_index(middle_of_path, &middleStep, path);    
    
    void step_action(int index, pathStep* ps)
    {
        if ((index >= start) && (index <= end)) {
            int cov = element_get_coverage_all_colours(ps->node);
            int match = 0;        
            int common_left, common_right, left_index, right_index;
            
            match = read_pair_score_nodes(&middleStep, ps, rpda, &common_left, &common_right, &left_index, &right_index);
            
            fprintf(bitfield_fp, "%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i", index, match, cov, common_left, common_right,
                    read_pair_count_bits(ps->node->signatures[Bitfield_First]),
                    read_pair_count_bits(ps->node->signatures[Bitfield_Second]),
                    read_pair_count_bits(ps->node->signatures[Bitfield_First_Inverse]),
                    read_pair_count_bits(ps->node->signatures[Bitfield_Second_Inverse]));
            
            for (i=0; i<4; i++) {
                for (j=0; j<4; j++) {
                    fprintf(bitfield_fp, "\t%i", read_pair_count_common_bits(middleStep.node->signatures[i], ps->node->signatures[j])); 
                }
            }
            fprintf(bitfield_fp, "\t%i\t%i\n", left_index, right_index);
        }
    }
    
    path_iterator_with_index(&step_action, path);    
    
    fclose(bitfield_fp);
}
#endif
#endif

ReadPairJump * read_pair_jump_new()
{
    ReadPairJump * rpj = malloc(sizeof(ReadPairJump));
    if(rpj == NULL) {
        fprintf(stderr, "Not enough memory to allocate new read_pair_jump");
        exit(-1);
    }
    return rpj;
}

void read_pair_jump_destroy(ReadPairJump ** rpj)
{
    free(*rpj);
    *rpj = NULL;
}

/**
 * Auxiliary function that invalidates a jump in both sides if it has been found. 
 */
#if 0
static void invalidate_pair_jump()
{
    
}
#endif

/**
 * Adds a jump to the read pair. It adds in both directions, so it can search from any pair. 
 */
void read_pair_jump_add_to_descriptor(pathStep * a, pathStep * b, int distance, ReadPairDescriptor * rpd)
{	
    assert(rpd != NULL);
    assert(a != NULL);
    assert(b != NULL);
    assert(rpd->jumps != NULL);
    
    ReadPairJump * rpj1 = read_pair_jump_new();
    ReadPairJump * rpj2 = read_pair_jump_new();
    
    path_step_assign(&rpj1->first, a);
    path_step_assign(&rpj1->second, b);
    
    path_step_assign(&rpj2->first, b);
    path_step_assign(&rpj2->second, a);
    
    BinaryTree * bt = rpd->jumps;
    
    ReadPairJump * tmp1 = (ReadPairJump * )binary_tree_find(rpj1, bt);
    ReadPairJump * tmp2 = (ReadPairJump * )binary_tree_find(rpj2, bt);
    
    
    if (tmp1 == NULL && tmp2 == NULL ) { //The pair hasn´t been found, so we don´t know
        binary_tree_add_element(rpj1, bt);
        binary_tree_add_element(rpj2, bt);
        // We have to validate if the found pair is the same, if not, invalidate all the corresponding jumps
    }else if (tmp1 == NULL) { //tmp1 is null, but tmp2 is not, we have to invalidate tmp2 and add the invalid pair to mark the ambiguity
        //TODO: Fill the complex cases
    }else if (tmp2 == NULL) { //tmp2 is null, but tmp1 is not, we have to invalidate tmp2 and add the invalid pair to mark the ambiguity
        //TODO: Fill the complex cases
    }else { //both exist, if we find that both pairs exist, but they are not linked, invalidate all the combinations. 
        //TODO: Fill the complex cases
    }
}

/**
 * Compares two jumps, to keep them in a set. It only compares the first step. In that way we can search for both sides
 * independently on the same sorted set. 
 */
int read_pair_jump_compare(void * a, void * b)
{	
    ReadPairJump * tmp1 = (ReadPairJump * ) a;
    ReadPairJump * tmp2 = (ReadPairJump * ) b;
    
    return path_step_compare(&tmp1->first, &tmp2->first);
    
}

ReadPairJump* read_pair_find_jump(ReadPairDescriptorArray* rpda, dBNode* node)
{
    int i;
    ReadPairJump* mate = NULL; 
    
    for (i=0; i<rpda->number_of_pairs; i++) {
        ReadPairDescriptor* rpd = rpda->pair[i];
        //printf("Pair: %i Jumps: %x\n", i, rpd->jumps);
        if (rpd->jumps != NULL) {
            TreeElement* match;
            match = binary_tree_find(node, rpd->jumps);
            if (match) {
                mate = match->element;
                printf("AMAZING: Mate found!\n");
                break;
            }
        }
    }
    
    return mate;
}



//Here comes the read pair walk and walking functions related to it!. 1


static boolean continue_traversing(pathStep * current_step,
                                   pathStep * next_step,
                                   pathStep * reverse_step, Path * temp_path,
                                   dBGraph * db_graph)
{
    //pathStep first;
    
    boolean cont = true;
    
    boolean next_visited = 
    (db_node_check_for_flag(next_step->node, VISITED_FORWARD) == true && next_step->orientation==forward )||
    (db_node_check_for_flag(next_step->node, VISITED_REVERSE) == true && next_step->orientation==reverse); //TODO: make this walk a bit forwad to see if the read pair is able to solve the little loop
    
    if (next_visited && next_step->label != Undefined) {
        path_action_set_flag(temp_path, IS_CYCLE);
        path_add_stop_reason(LAST, PATH_FLAG_IS_CYCLE, temp_path);
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

static long long count_solved_x_nodes = 0;
static long long count_unsolved_x_nodes = 0;
static long long count_solved_short_x_path = 0;
static long long count_unsolved_short_x_path = 0;
static long long count_solved_long_x_path = 0;
static long long count_unsolved_long_x_path = 0;




//TODO make a smaller, more readable signature. 
static Nucleotide find_best_nucleatide_in_double_y(pathStep * current_step, pathStep * next_step, pathStep * reverse_step, dBGraph * db_graph){
    Nucleotide best = Undefined;
    int i, j;
   // ReadPairSignature current_rps = current_step->node
    ReadPairDescriptorArray * rpda = db_graph->rpda;
   //int last_in_index = path_get_index_of_last_in_node(current_step->path);
    int pos = current_step->path->length-rpda->pair[0]->distance + rpda->pair[0]->tolerance;
    if (pos <= 0) {
        return Undefined;
    }
    int last_in_index =  path_get_first_in_node_after(pos, current_step->path);
    if (last_in_index != -1) {
        if(last_in_index > 0){
            last_in_index--;
        }
        dBNode * neighbours[4];
        Nucleotide labels[4];
        int found = db_graph_get_neighbouring_nodes_all_colours(next_step->node, next_step->orientation, neighbours, labels, db_graph);
        pathStep compare_to;
       // path_get_step_at_index(last_in_index, &compare_to, current_step->path);
        
        ReadPairSignature rps = ~0;
        for (i = last_in_index; i < path_get_length(current_step->path); i++) {
            path_get_step_at_index(i, &compare_to, current_step->path);
            rps &= db_node_get_signature(0, Bitfield_First, compare_to.node); //TODO: make this for multiple read pairs
        }
        if(read_pair_count_bits(rps) > 2 *( READ_PAIR_LENGTH / 3) ){
            return Undefined;
        }
        boolean invalid = true;
        int  best_matches = 0, temp_matches, second_best_matches = 0;
        for (i = 0; i < found; i++) {
            temp_matches = 0;
            for(j = 0; j < NUMBER_OF_SIGNATURES; j++){
                temp_matches += read_pair_count_common_bits(rps, neighbours[i]->signatures[j]);//TODO: make this line cleaner and compatible with multiple read pairs
                
            }
            if(temp_matches > best_matches && temp_matches > db_graph->rpda->pair[0]->min_pair_coverage){//TODO: make this a variable. 
                second_best_matches = best_matches;
                best_matches = temp_matches;
                best = labels[i];
                invalid = false;
            }else if (temp_matches == best_matches) {
                invalid = true;
                best = Undefined;
            }
        }
        
        if(best_matches - second_best_matches < db_graph->rpda->pair[0]->min_pair_coverage){
            best = Undefined;
        }
        
    }
    
    return  best;
}

static pathStep *get_next_step(pathStep * current_step, pathStep * next_step,
                               pathStep * reverse_step, dBGraph * db_graph)
{
	pathStep *step = db_graph_get_next_step(current_step, next_step, reverse_step,
                                            db_graph);
	assert(step != NULL);
    Nucleotide n;
	int count_fwd, count_rev;
    if (step->node != NULL) {
        
		step->label = Undefined;
        count_fwd = db_node_edges_count(next_step->node, next_step->orientation);
        count_rev = db_node_edges_count(next_step->node, opposite_orientation(next_step->orientation));
		if(count_fwd == 0){
            step->label = Undefined;//We cant do anything. 
        }else if (db_node_has_precisely_one_edge(next_step->node, next_step->orientation, &n)) {//Simple case, no ambiguity, we don't need to do anything else. 
			next_step->label = n;
		}else if(count_fwd > 1 && count_rev > 1 ){ //This is an X node, try to resolve with the own reads signature. If it is not greater than 1, we need to get the previous
            //   path_get_step_at_index(path_get_length(current_path) - 2, &prevStep, current_path);
            next_step->label = db_graph_get_best_next_step_nucleotide(next_step->node, current_step->node, next_step->orientation, db_graph);
            if(next_step->label != Undefined){
                count_solved_x_nodes++;
            }else{
                count_unsolved_x_nodes++;
            }
        }
        if( step->label == Undefined && count_fwd > 1){//Here we try to do the read pair magic. 
            next_step->label = find_best_nucleatide_in_double_y(current_step, next_step, reverse_step,db_graph);
            if(next_step->label != Undefined){
                count_solved_short_x_path++;
            }else{
                count_unsolved_short_x_path++;
            }
        }//This are all the cases. 
        
	}
    
	return next_step;
}

Path * read_pair_get_path_single_walk(dBNode * node, Orientation orientation,
                                      void (*node_action) (dBNode * node),
                                      dBGraph * db_graph, boolean both_directions,Path * path)
{
    
    Path *buff1 = path_get_buffer_path();
    assert(buff1 != path);
    path_reset(buff1);
    path_reset(path);	//Make sure that the paths are clean
    
    pathStep first;
    boolean only_one_edge;
    first.node = node;
    first.orientation =  orientation;
    
    WalkingFunctions wf;
    perfect_path_get_funtions(&wf);
    assert(buff1 != path);
    void (* old_step_action)(pathStep * st) = wf.step_action;
    void step_action (pathStep * step){
        old_step_action(step);
        node_action(step->node);
    }
    
    
    wf.continue_traversing = &continue_traversing;
    
    wf.get_next_step =  &get_next_step;
    wf.step_action = &step_action;
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
    void callback(Path * p) {
        // We have to copy the path, either in reverse or as it is, otherwise we lose it.
        if (both_directions == true) {
            path_reverse(p, path);
        } else {
            path_copy(path, p);
        }
    }
    assert(buff1 != path);
    wf.output_callback = &callback;
    assert(buff1 != path);
    path_reset(buff1);
    if (db_node_edges_count(first.node, first.orientation) > 0) {	
        assert(buff1 != path);
        db_graph_generic_walk(&first, buff1, &wf, db_graph);
        assert(buff1 != path);
    } else {
        if (both_directions == true) {
            first.label = undefined;
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


static void path_step_mark_as_visited(pathStep * ps){
    db_node_action_set_flag_visited(ps->node);
}

void read_pair_print_paths_single_walk(char *filename, int max_length,boolean with_coverages, boolean with_viz, ReadPairDescriptorArray * rpda,dBGraph * db_graph){
    
    int singleton_length = db_graph->kmer_size *2; //TODO: make this a parameter
    Path *path = path_get_buffer_path();	
    path_reset(path);
    
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
    long long count_repeated = 0;
    double graph_cov = db_graph_get_average_coverage(db_graph);
    log_and_screen_printf("Average coverage: %5.2f \n", graph_cov);
    
    PathCounts counts;
    path_counts_reset(&counts); 
    
    
    void print_supernode(dBNode * node) {
        
		count_kmers++;
		if (db_node_check_flag_visited(node) == false && db_node_check_flag_not_pruned(node)) {
			read_pair_get_path_single_walk(node, undefined, &db_node_action_do_nothing,db_graph, true, path);
            
			if (path_is_singleton(singleton_length, path)) {
                count_sing++;
            }else  if(path_is_repetitive(graph_cov, path)){
                count_rep++;
            }else if(path_percentage_new_nodes(path) < 50){
                count_repeated++;
                path_iterator(path_step_mark_as_visited, path);
            }else{
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
                path_iterator(path_step_mark_as_visited, path);
				path_counts_add(path, &counts);
				count_nodes++;
				path_increase_id(path);
			}
		}
	}
    
	hash_table_traverse(&print_supernode, db_graph);
	log_and_screen_printf("%'d nodes visited [%'qd singletons, %'qd repetitive, %'qd less than 50%% new]\n", count_nodes, count_sing, count_repeated);
    log_and_screen_printf("solved x-nodes %'lld\n unsolved x-nodes: %'lld\n", count_solved_x_nodes, count_unsolved_x_nodes);
    log_and_screen_printf("solved x-paths %'lld\n unsolved x-paths: %'lld\n", count_solved_short_x_path, count_unsolved_short_x_path);
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
}



static void or_sign(pathStep * ps,  void * argsa)
{
    sign_args * args = (sign_args *) argsa;
    args->rps |= db_node_get_signature(args->pair, args->index, ps->node);
}

static void and_sign(pathStep * ps,  void * argsa)
{
    sign_args * args = (sign_args *) argsa;
    args->rps &= db_node_get_signature(args->pair, args->index, ps->node);
}

static ReadPairSignature get_or_signature_across_path(int pair, int sign_index,Path * p){
    sign_args sign;
    sign.index = sign_index;
    sign.rps = 0; // all 0
    sign.pair = pair;
    path_iterator_with_args(&or_sign, &sign, p);
    return sign.rps;
}
static ReadPairSignature get_and_signature_across_path(int pair, int sign_index,Path * p){
    sign_args sign;
    sign.index = sign_index;
    sign.rps = ~0; // all 1
    sign.pair = 0;
    path_iterator_with_args(&and_sign, &sign, p);
    return sign.rps;
}

static ReadPairSignature get_or_signature_across_all_paths(int pair, int sign_index,PathArray * pa){
    ReadPairSignature sign = 0;
    int i;
    for (i = 0; i < path_array_get_number_of_paths(pa); i++){
        sign |= get_or_signature_across_path(pair,sign_index, path_array_get(i, pa));
    }
    
    return sign;

}

ReadPairSignature get_free_signature_from_bitfield(ReadPairSignature rps){
    int i;
    int shift;
    boolean found = false;
    size_t sign_size = READ_PAIR_LENGTH;
    long long offset = rand() % sign_size;
    ReadPairSignature sign = rps;
    ReadPairSignature offseted_signature = sign >> offset;
    for (i = 0; i < sign_size - offset && !found; i++) {
        if((1 & offseted_signature) == 0){
            found = true;
            shift = i + offset;
        }
        offseted_signature >>= 1;
    }
    for (i = 0; !found && i < offset; i++) {
        if ((1 & offseted_signature) == 0 ) {
            found = true;
            shift = i;
        }
        sign >>= 1;
    }
  /*  if(found == false){
        printf("not found for %x", (unsigned int)rps);
    }*/
    return found?1<<shift:0;
    
}

static ReadPairSignature get_free_signature(Path * p,int index){
    ReadPairSignature sign = get_or_signature_across_path(0, index, p );
    return get_free_signature_from_bitfield(sign);
}

int read_pair_get_maximum_insert_size(ReadPairDescriptorArray * rpda){
    ReadPairDescriptor * rpd; 
    int pairs = rpda->number_of_pairs;
    assert(pairs > 0);
    rpd = rpda->pair[0];
    assert(rpd != NULL);
    int max = rpd->insert_size;
    int i;
    for (i = 1; i < pairs; i++) {
        rpd = rpda->pair[i];
        if(max < rpd->insert_size){
            max =rpd->insert_size;
        }
    }
    
    return max;
}

int read_pair_get_minimum_insert_size(ReadPairDescriptorArray * rpda){

    ReadPairDescriptor * rpd; 
    int pairs = rpda->number_of_pairs;
    assert(pairs > 0);
    rpd = rpda->pair[0];
    assert(rpd != NULL);
    int min = rpd->insert_size;
    int i;
    for (i = 1; i < pairs; i++) {
        rpd = rpda->pair[i];
        if(min > rpd->insert_size){
            min =rpd->insert_size;
        }
    }
    
    return min;
}

static ReadPairSignature get_or_signature_of_flanking_paths(int pair, int sign_index, int max_distance, Path * path, dBGraph * db_graph){
    ReadPairSignature rps = 0;
    pathStep first;
    pathStep last;
    path_get_last_step(&last, path);
    path_get_step_at_index(0, &first, path);
    first.orientation = opposite_orientation(first.orientation);
    first.label = Undefined;
    last.label = Undefined;
    
    void or_local_path(Path * p){
        rps |= get_or_signature_across_path(pair, sign_index, p);
    }
    branches_get_all_paths_from_with_callback(&first, max_distance, &or_local_path, db_graph);
    branches_get_all_paths_from_with_callback(&last, max_distance, &or_local_path, db_graph);
    
    return rps;
}


int read_pair_mark_adjacency(Path * path, dBGraph * db_graph, int id){
    pathStep first;
    pathStep last;
    int marked = 0;
    long long path_collisions = 0;
    path_get_last_step(&last, path);
    path_get_step_at_index(0, &first, path);
    first.orientation = opposite_orientation(first.orientation);
    first.label = Undefined;
    last.label = Undefined;
    int i, j;
    //PathArray * right_paths = branches_get_all_paths_from(&last, max_distance, db_graph);
    //PathArray * left_paths = branches_get_all_paths_from(&first, max_distance, db_graph);
    //path_array_merge(&left_paths, right_paths);
   // ReadPairSignature free_signs = get_or_signature_across_all_paths(0, Bitfield_PerfectPath, right_paths);
   // ReadPairSignature free_signs = get_or_signature_of_flanking_paths(0, Bitfield_PerfectPath, max_distance, path, db_graph);
   // ReadPairSignature next_sign = get_free_signature_from_bitfield(free_signs);
   // db_node_action_add_read_pair(<#long long count#>, <#short pair#>, <#ReadPairBitfield bitfield#>, <#dBNode *node#>)
    //TODO: make this in such a way you make sure that both sides are marked, will it be more efficient?
    ReadPairSignature next_sign = id;
    void enrich_path(pathStep *ps) {
        // db_node_action_add_read_pair(label, pair, Bitfield_First, ps->node);
       db_node_set_signature(next_sign, 0, Bitfield_PerfectPath, ps->node);
    }
    
    //  printf("full Linking path\n");
    // path_to_fasta(full_linking_path, stdout);
    
    
    if (next_sign) {
        marked++;
        path_iterator(&enrich_path, path);
        /*for (i = 0; i < path_array_get_number_of_paths(right_paths); i++) {
            left_paths = path_split_in_perfect_paths(path_array_get(i, right_paths));
            for (j = 0; next_sign && j < path_array_get_number_of_paths(left_paths);j++) {
                free_signs |= next_sign;
                next_sign = get_free_signature_from_bitfield(free_signs); 
                if (next_sign) {
                    marked++;
                    path_iterator(&enrich_path,  path_array_get(j, left_paths) );
                }
                
                
            }   
            path_array_free_from_buffer(left_paths);
            left_paths = NULL;
        }*/
    }else{
        
        path_collisions++; 
        
    }
    path_mark_as_visited(path);
    //log_and_screen_printf("Perfect path saturated paths: %lld\n", path_collisions);
     //path_array_free_from_buffer(right_paths);
   
    return marked; 
}


static void mark_adjacency_from_node(dBNode * node, void * argsa){
    mark_args * args = (mark_args *)argsa;
    Path * path;
    int len;
    int mark;
    if (!db_node_check_flag_visited(node)) {
        
        path = path_get_buffer_path();
        len = perfect_path_get_path(node, undefined, &db_node_action_do_nothing, args->graph, path );
        if(len > 0){
            mark = read_pair_mark_adjacency(path, args->graph, args->count +1);
            if(mark)
                args->count += mark;    
            else
                args->fail_count++;
        }
        path_free_buffer_path(path);
        path = NULL;
    }    
    
}

void read_pair_mark_contiguous_perfect_paths(int max_distance, dBGraph * db_graph){
    
    int threads =  1, i;
#ifdef THREADS
    threads = db_graph->number_of_threads;
#endif
    long long total_count = 0;
    long long total_fail = 0;
    mark_args ** args = calloc(threads, sizeof(mark_args  *));
    for (i=0; i<threads; i++) {
        args[i] = calloc(1, sizeof(mark_args ));
        args[i]->graph = db_graph;
        args[i]->max_distance = max_distance;
        args[i]->count = 0;
        args[i]->fail_count = 0;
    }
    log_write_timestamp(true);
    log_and_screen_printf("Marking perfect paths\n");
    hash_table_traverse_with_args(&mark_adjacency_from_node, (void **)args, db_graph);
    
    for (i=0; i < threads; i++) {
        total_count +=  args[i]->count;
        total_fail  +=  args[i]->fail_count;
        free(args[i]);
    }
    
    log_and_screen_printf("Marked %lld perfect paths, not marked %lld perfect paths\n", total_count, total_fail);
    free(args);

}

void read_pair_search_enrich_graph(FILE * f1, FILE * f2, ReadPairDescriptor * rpd, int fastq_ascii_offset, dBGraph * db_graph)
{
	int max_read_length = rpd->max_read_length + 1; 
	Sequence * seq1 = sequence_new(max_read_length, LINE_MAX, 0);
	Sequence * seq2 = sequence_new(max_read_length, LINE_MAX, 0);
	int length_1, length_2, nkmers1, nkmers2;
    int missed_kmers_1, missed_kmers_2;
	short kmer_size = db_graph->kmer_size;

    
    
	//boolean search_jumps = rpd->jumps == NULL? false: true; //Validate if the descriptor has space for the index of jumps. 
	
	//max_read_length/(kmer_size+1) is the worst case for the number of sliding windows, ie a kmer follow by a low-quality/bad base
	int max_windows = max_read_length / (kmer_size + 1);
	
	//number of possible kmers in a 'perfect' read
	int max_kmers = max_read_length - kmer_size + 1;
	
	KmerSlidingWindowSet * windows = binary_kmer_sliding_window_set_new(max_windows, max_kmers);
	windows->kmer_size = kmer_size;
	
	BinaryKmer tmp_key;
	int quality_cut_off = 0; //The cutoff comes from the previous CTXs, or the cleaning algorithm. 
	Path * path_read_one = path_get_buffer_path();
	Path * path_read_two = path_get_buffer_path();
    Path * linking_path = path_get_buffer_path();
    Path * full_linking_path = path_get_buffer_path();
    Path * linking_path_long = path_get_buffer_path();
    
	dBNode * current_node;
	pathStep ps;
	ps.label = Undefined;
	ps.orientation = forward;
    ps.flags = 0;
    long long seq_count = 0;
	long long ignore_count = 0;
	long long used_count = 0;
    long long on_perfect_path = 0;
    long long label_added = 0;
    long long not_found_label = 0;
    long long no_path_found = 0;
    long long already_connected = 0;
	long long not_perfect_path = 0;
    long long too_short = 0;
    //double coverage = db_graph_get_average_coverage(db_graph);
    //int cov = (int) coverage;
    int min_common_signatures = rpd->max_pair_coverage;//Make this dynamic? use the coverage? Seems the less the better! 
//int min_common_signatures = cov;//Make this dynamic? use the coverage?
	while((length_1 = read_sequence_from_fastq(f1, seq1, rpd->max_read_length) )&& (length_2 = read_sequence_from_fastq(f2, seq2, rpd->max_read_length))){
		seq_count++;
        path_reset(path_read_one);
        path_reset(path_read_two);
        path_reset(linking_path);
        path_reset(linking_path_long);
        path_reset(full_linking_path);
        if ((seq_count % 1000000) == 0) {
            log_write_timestamp(1);
            log_and_screen_printf("Enriched %lli sequences\n", seq_count);
            log_and_screen_printf("Ignored reads: %lli Used reads: %lli, in perfect path:: %lli, labels added: %lli, not_found_label: %lli, path_too_short %lli, not_path_found %lli, already_connected: %lli , not_perfect_path: %lli, used_average_distance: %lli, used_average_depth: %lli\n", ignore_count, used_count, on_perfect_path, label_added, not_found_label, too_short, no_path_found, already_connected, not_perfect_path, rpd->insert_distance_sum/rpd->used_pairs, rpd->total_depth/rpd->used_pairs );  
        }
        
		path_reset(path_read_one);
        missed_kmers_1 = 0;
        nkmers1 = get_sliding_windows_from_sequence(seq1->seq, seq1->qual, seq1->length, quality_cut_off, kmer_size, windows, max_windows, max_kmers,false,0);					
        
		while (binary_kmer_sliding_window_set_get_next(windows)) {
			element_get_key(&windows->current,kmer_size, &tmp_key); 		
			current_node = hash_table_find(&tmp_key, db_graph);
			if (current_node != NULL) {
				ps.node = current_node;
				if (binary_kmer_comparison_operator(tmp_key,windows->current)) {
					ps.orientation = forward;
				} else {
					ps.orientation = reverse;
				}
				path_add_node(&ps, path_read_one);
			} else {
                missed_kmers_1++;
            }
		}
		
		path_reset(path_read_two);
        missed_kmers_2 = 0;
        nkmers2 = get_sliding_windows_from_sequence(seq2->seq, seq2->qual, seq2->length, quality_cut_off, kmer_size, windows, max_windows, max_kmers, false, 0);
        
		while (binary_kmer_sliding_window_set_get_next(windows)) {
			element_get_key(&windows->current,kmer_size, &tmp_key); 
			current_node = hash_table_find(&tmp_key, db_graph);
			if(current_node != NULL) {
				ps.node = current_node;
				
				if (binary_kmer_comparison_operator(tmp_key, windows->current)) {
					ps.orientation = forward;
				} else {
					ps.orientation = reverse;
				}
				
				path_add_node(&ps, path_read_two);
			} else {
                missed_kmers_2++;
            }
		}
        
        if ((path_get_length(path_read_one) > 0) && (path_get_length(path_read_two) > 0) && (missed_kmers_1 <= rpd->allowed_misses) && (missed_kmers_2 <= rpd->allowed_misses)) {
            used_count++;
        } else {
            ignore_count++;
        }
        
		if ((path_get_length(path_read_one) > 0) && (path_get_length(path_read_two) > 0) && (missed_kmers_1 <= rpd->allowed_misses) && (missed_kmers_2 <= rpd->allowed_misses)) { 
			// Validating that both reads are still in the graph. If they are not, it might mean that it makes no sense to enrich the pair. 
			
            
            pathStep first;
            pathStep last;
            path_get_step_at_index(0, &first, path_read_one);
            path_get_step_at_index(0, &last, path_read_two);
            //path_get_last_step(&last, path_read_two);
            short pair = 0; //Do we really want many graphs for different insert sizes.
          //  ReadPairSignature rps = db_node_get_signature(pair, Bitfield_First, first.node) & db_node_get_signature(pair, Bitfield_First, first.node);
            if (db_node_get_signature(pair, Bitfield_PerfectPath, first.node) == db_node_get_signature(pair, Bitfield_PerfectPath, last.node)) {
               on_perfect_path++;  
            }else if (read_pair_count_common_bits(db_node_get_signature(pair, Bitfield_First, first.node), db_node_get_signature(pair, Bitfield_First, last.node)) < min_common_signatures) {//We only search them if they don't seem connected. We add up to 3 signatures to search. 
                
                
                 Path * the_path = branches_get_path_between_paths(path_read_one, path_read_two, rpd->distance + rpd->tolerance, rpd->depth,linking_path, db_graph );
               
                //TODO: get the perfect paths from the first node, last node and concatenate them for the signature. 
                //  if(linking_path->in_nodes_count > 0 && linking_path->out_nodes_count > 0){
               // if(the_path != NULL){
                if(the_path != NULL && the_path->in_nodes_count > 0 && the_path->out_nodes_count > 0){
                    //only enrich if this solves a double Y
                    
                    int distance  = path_get_length(the_path);
                    if (distance < rpd->distance - rpd->tolerance) {
                        too_short++;
                    }else{
                    rpd->total_depth += the_path->out_nodes_count;
                    rpd->insert_distance_sum += distance;
                    rpd->used_pairs++;
                    path_append(linking_path_long, linking_path);
                                    
                   void append_forward(Path * tmp_path){                        
                        path_append(linking_path_long, tmp_path);
                    }
                    
                    void append_backward(Path * tmp_path2){
                        path_reverse(tmp_path2, full_linking_path);
                        path_append(full_linking_path, linking_path);
                    }
                    
                    
                    pathStep lasto;
                    path_get_last_step( &lasto,linking_path );
                    if(db_node_has_precisely_one_edge(lasto.node, lasto.orientation, &lasto.label)){
                        perfect_path_get_path_from_step_with_callback(lasto, &db_node_action_do_nothing, &append_forward, db_graph);
                    }
                    
                    pathStep firsto;
                    path_get_step_at_index(0, &firsto,linking_path );
                    firsto.orientation = opposite_orientation(firsto.orientation);
                    if(db_node_has_precisely_one_edge(firsto.node, firsto.orientation, &firsto.label)){
                        perfect_path_get_path_from_step_with_callback(firsto, &db_node_action_do_nothing, &append_backward, db_graph);
                    }else{
                        path_copy(full_linking_path, linking_path);
                    }
                    ReadPairSignature label = get_free_signature( full_linking_path, Bitfield_First);//to avoid overflowing the signature.    
                    
                    void enrich_path(pathStep *ps) {
                        db_node_set_signature(label, pair, Bitfield_First, ps->node);
                    }
                    

                    if (label) {
                        path_iterator(&enrich_path, full_linking_path);
                        label_added++;
                    }else{
                        not_found_label++;
                    }
                    
                    not_perfect_path++;
                    }
                }else{
                    no_path_found++;
                }
            }else{
                already_connected++;
            }
		}
	}
    path_free_buffer_path(linking_path_long);
	path_free_buffer_path(path_read_one);
    path_free_buffer_path(path_read_two);
    path_free_buffer_path(linking_path);
    path_free_buffer_path(full_linking_path);
    
	binary_kmer_free_kmers_set(&windows);
	free_sequence(&seq1);
	free_sequence(&seq2);
    
    log_and_screen_printf("Read %lli sequences\n", seq_count);
    log_write_timestamp(1);
    log_and_screen_printf("Enriched %lli sequences\n", seq_count);
    log_and_screen_printf("Ignored reads: %lli Used reads: %lli, in perfect path:: %lli, labels added: %lli, not_found_label: %lli, not_path_found %lli, already_connected: %lli , not_perfect_path: %lli\n", ignore_count, used_count, on_perfect_path, label_added, not_found_label, no_path_found, already_connected, not_perfect_path);  
}

void read_pair_enrich_graph_of_graph(FILE * f1, FILE * f2, ReadPairDescriptor * rpd, int fastq_ascii_offset, dBGraph * db_graph)
{
	int max_read_length = rpd->max_read_length + 1; 
	Sequence * seq1 = sequence_new(max_read_length, LINE_MAX, 0);
	Sequence * seq2 = sequence_new(max_read_length, LINE_MAX, 0);
	int length_1, length_2, nkmers1, nkmers2;
    int missed_kmers_1, missed_kmers_2;
	short kmer_size = db_graph->kmer_size;
//    int number_of_perfect_paths = db_g
    
    
	int max_windows = max_read_length / (kmer_size + 1);
	
	//number of possible kmers in a 'perfect' read
	int max_kmers = max_read_length - kmer_size + 1;
	
	KmerSlidingWindowSet * windows = binary_kmer_sliding_window_set_new(max_windows, max_kmers);
	windows->kmer_size = kmer_size;
	
	BinaryKmer tmp_key;
	int quality_cut_off = 0; //The cutoff comes from the previous CTXs, or the cleaning algorithm. 
	Path * path_read_one = path_get_buffer_path();
	Path * path_read_two = path_get_buffer_path();
    Path * linking_path = path_get_buffer_path();
    Path * full_linking_path = path_get_buffer_path();
    Path * linking_path_long = path_get_buffer_path();
    
	dBNode * current_node;
	pathStep ps;
	ps.label = Undefined;
	ps.orientation = forward;
    ps.flags = 0;
    long long seq_count = 0;
	long long ignore_count = 0;
	long long used_count = 0;
    long long on_perfect_path = 0;
    long long label_added = 0;
    long long not_found_label = 0;
    long long no_path_found = 0;
    long long already_connected = 0;
	long long not_perfect_path = 0;
    long long too_short = 0;
    //double coverage = db_graph_get_average_coverage(db_graph);
    //int cov = (int) coverage;
    int min_common_signatures = rpd->max_pair_coverage;//Make this dynamic? use the coverage? Seems the less the better! 
    //int min_common_signatures = cov;//Make this dynamic? use the coverage?
	while((length_1 = read_sequence_from_fastq(f1, seq1, rpd->max_read_length) )&& (length_2 = read_sequence_from_fastq(f2, seq2, rpd->max_read_length))){
		seq_count++;
        path_reset(path_read_one);
        path_reset(path_read_two);
        path_reset(linking_path);
        path_reset(linking_path_long);
        path_reset(full_linking_path);
        if ((seq_count % 1000000) == 0) {
            log_write_timestamp(1);
            log_and_screen_printf("Enriched %lli sequences\n", seq_count);
            log_and_screen_printf("Ignored reads: %lli Used reads: %lli, in perfect path:: %lli, labels added: %lli, not_found_label: %lli, path_too_short %lli, not_path_found %lli, already_connected: %lli , not_perfect_path: %lli, used_average_distance: %lli, used_average_depth: %lli\n", ignore_count, used_count, on_perfect_path, label_added, not_found_label, too_short, no_path_found, already_connected, not_perfect_path, rpd->insert_distance_sum/rpd->used_pairs, rpd->total_depth/rpd->used_pairs );  
        }
        
		path_reset(path_read_one);
        missed_kmers_1 = 0;
        nkmers1 = get_sliding_windows_from_sequence(seq1->seq, seq1->qual, seq1->length, quality_cut_off, kmer_size, windows, max_windows, max_kmers,false,0);					
        
		while (binary_kmer_sliding_window_set_get_next(windows)) {
			element_get_key(&windows->current,kmer_size, &tmp_key); 		
			current_node = hash_table_find(&tmp_key, db_graph);
			if (current_node != NULL) {
				ps.node = current_node;
				if (binary_kmer_comparison_operator(tmp_key,windows->current)) {
					ps.orientation = forward;
				} else {
					ps.orientation = reverse;
				}
				path_add_node(&ps, path_read_one);
			} else {
                missed_kmers_1++;
            }
		}
		
		path_reset(path_read_two);
        missed_kmers_2 = 0;
        nkmers2 = get_sliding_windows_from_sequence(seq2->seq, seq2->qual, seq2->length, quality_cut_off, kmer_size, windows, max_windows, max_kmers, false, 0);
        
		while (binary_kmer_sliding_window_set_get_next(windows)) {
			element_get_key(&windows->current,kmer_size, &tmp_key); 
			current_node = hash_table_find(&tmp_key, db_graph);
			if(current_node != NULL) {
				ps.node = current_node;
				
				if (binary_kmer_comparison_operator(tmp_key, windows->current)) {
					ps.orientation = forward;
				} else {
					ps.orientation = reverse;
				}
				
				path_add_node(&ps, path_read_two);
			} else {
                missed_kmers_2++;
            }
		}
        
        if ((path_get_length(path_read_one) > 0) && (path_get_length(path_read_two) > 0) && (missed_kmers_1 <= rpd->allowed_misses) && (missed_kmers_2 <= rpd->allowed_misses)) {
            used_count++;
        } else {
            ignore_count++;
        }
        
		if ((path_get_length(path_read_one) > 0) && (path_get_length(path_read_two) > 0) && (missed_kmers_1 <= rpd->allowed_misses) && (missed_kmers_2 <= rpd->allowed_misses)) { 
			// Validating that both reads are still in the graph. If they are not, it might mean that it makes no sense to enrich the pair. 
			
            
            pathStep first;
            pathStep last;
            path_get_step_at_index(0, &first, path_read_one);
            path_get_step_at_index(0, &last, path_read_two);
            //path_get_last_step(&last, path_read_two);
            short pair = 0; //Do we really want many graphs for different insert sizes.
            //  ReadPairSignature rps = db_node_get_signature(pair, Bitfield_First, first.node) & db_node_get_signature(pair, Bitfield_First, first.node);
            if (db_node_get_signature(pair, Bitfield_PerfectPath, first.node) == db_node_get_signature(pair, Bitfield_PerfectPath, last.node)) {
                on_perfect_path++;  
            }else if (read_pair_count_common_bits(db_node_get_signature(pair, Bitfield_First, first.node), db_node_get_signature(pair, Bitfield_First, last.node)) < min_common_signatures) {//We only search them if they don't seem connected. We add up to 3 signatures to search. 
                
                
                Path * the_path = branches_get_path_between_paths(path_read_one, path_read_two, rpd->distance + rpd->tolerance, rpd->depth,linking_path, db_graph );
                
                //TODO: get the perfect paths from the first node, last node and concatenate them for the signature. 
                //  if(linking_path->in_nodes_count > 0 && linking_path->out_nodes_count > 0){
                // if(the_path != NULL){
                if(the_path != NULL && the_path->in_nodes_count > 0 && the_path->out_nodes_count > 0){
                    //only enrich if this solves a double Y
                    
                    int distance  = path_get_length(the_path);
                    if (distance < rpd->distance - rpd->tolerance) {
                        too_short++;
                    }else{
                        rpd->total_depth += the_path->out_nodes_count;
                        rpd->insert_distance_sum += distance;
                        rpd->used_pairs++;
                        path_append(linking_path_long, linking_path);
                        
                        void append_forward(Path * tmp_path){                        
                            path_append(linking_path_long, tmp_path);
                        }
                        
                        void append_backward(Path * tmp_path2){
                            path_reverse(tmp_path2, full_linking_path);
                            path_append(full_linking_path, linking_path);
                        }
                        
                        
                        pathStep lasto;
                        path_get_last_step( &lasto,linking_path );
                        if(db_node_has_precisely_one_edge(lasto.node, lasto.orientation, &lasto.label)){
                            perfect_path_get_path_from_step_with_callback(lasto, &db_node_action_do_nothing, &append_forward, db_graph);
                        }
                        
                        pathStep firsto;
                        path_get_step_at_index(0, &firsto,linking_path );
                        firsto.orientation = opposite_orientation(firsto.orientation);
                        if(db_node_has_precisely_one_edge(firsto.node, firsto.orientation, &firsto.label)){
                            perfect_path_get_path_from_step_with_callback(firsto, &db_node_action_do_nothing, &append_backward, db_graph);
                        }else{
                            path_copy(full_linking_path, linking_path);
                        }
                        ReadPairSignature label = get_free_signature( full_linking_path, Bitfield_First);//to avoid overflowing the signature.    
                        
                        void enrich_path(pathStep *ps) {
                            db_node_set_signature(label, pair, Bitfield_First, ps->node);
                        }
                        
                        
                        if (label) {
                            path_iterator(&enrich_path, full_linking_path);
                            label_added++;
                        }else{
                            not_found_label++;
                        }
                        
                        not_perfect_path++;
                    }
                }else{
                    no_path_found++;
                }
            }else{
                already_connected++;
            }
		}
	}
    path_free_buffer_path(linking_path_long);
	path_free_buffer_path(path_read_one);
    path_free_buffer_path(path_read_two);
    path_free_buffer_path(linking_path);
    path_free_buffer_path(full_linking_path);
    
	binary_kmer_free_kmers_set(&windows);
	free_sequence(&seq1);
	free_sequence(&seq2);
    
    log_and_screen_printf("Read %lli sequences\n", seq_count);
    log_write_timestamp(1);
    log_and_screen_printf("Enriched %lli sequences\n", seq_count);
    log_and_screen_printf("Ignored reads: %lli Used reads: %lli, in perfect path:: %lli, labels added: %lli, not_found_label: %lli, not_path_found %lli, already_connected: %lli , not_perfect_path: %lli\n", ignore_count, used_count, on_perfect_path, label_added, not_found_label, no_path_found, already_connected, not_perfect_path);  
}






#endif
