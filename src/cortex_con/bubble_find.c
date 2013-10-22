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
 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <binary_kmer.h>
#include <element.h>
#include <dB_graph.h>
#include <bubble_find.h>
#include <logger.h>

#define MAX_PATHS_IN_ARRAY 4096

#ifdef DEBUG_BUBBLES
int BUBBLEDEBUG = 1;
#else
int BUBBLEDEBUG = 0;
#endif

// ----------------------------------------------------------------------
// Identify branch points of potential SNPs
// ----------------------------------------------------------------------
void db_graph_identify_branches(int max_length, dBGraph * db_graph)
{
	int branchNodes = 0;
	
	db_graph_reset_flags(db_graph);

	// Hash table iterator to label nodes
	void find_end_nodes(dBNode * node) {
		if (db_node_check_flag_not_pruned(node)) {
			boolean branch = false;

			// Look for Y shape branch forward orientation
			// The nodes at the top of the Y should contain different colours
			if (db_node_edges_count_all_colours(node, forward) > 1
			    && db_node_edges_count_all_colours(node, reverse) == 1) {
				db_node_action_set_flag(node, BRANCH_NODE_FORWARD);
				branchNodes++;
				branch = true;
			}
			// Look for Y shape branch reverse orientation                     
			if (db_node_edges_count_all_colours(node, reverse) > 1
			    && db_node_edges_count_all_colours(node, forward) == 1) {
				db_node_action_set_flag(node, BRANCH_NODE_REVERSE);
				branchNodes++;
				branch = true;
			}
			// Look for X-shaped branch
			if (db_node_edges_count_all_colours(node, reverse) > 1
			    && db_node_edges_count_all_colours(node, forward) > 1) {
				db_node_action_set_flag(node, X_NODE);
				branchNodes++;
				branch = true;
			}

			if (DEBUG) {
				if (branch == true) {
					char node_kmer[db_graph->kmer_size + 1];
					binary_kmer_to_seq(element_get_kmer(node), db_graph->kmer_size, node_kmer);
					printf("Marked kmer %s\n", node_kmer);
				}
			}
		}
	}

	if (DEBUG) {
		printf("[db_graph_identify_branches] Finding branches...\n");
	}
	hash_table_traverse(&find_end_nodes, db_graph);

	if (DEBUG) {
		printf("[db_graph_identify_branches] Branch point identification complete: %d nodes labelled.\n", branchNodes);
	}
}

// ----------------------------------------------------------------------
// Add a path to the path array, after checking it's not circular
// ----------------------------------------------------------------------
void db_graph_check_and_add_path(Path * path, PathArray * patharray)
{
	if ((path->length > 1) && (path->nodes[0] == path->nodes[path->length - 1])) {
		printf("Circular path - not adding\n");
	} else {
		if (patharray->number_of_paths == patharray->capacity) {
			printf("Warning: Trying to add a path to an array which has reached capacity (%d). Path ignored, but some bubbles may be missed.\n", patharray->capacity);
		} else {
			// Copy the path - because this one will get destroyed when we go back down the recursion stack
			Path *p = path_new(path->length + 1, path->kmer_size);
			path_copy(p, path);
			path_array_add_path(p, patharray);
		}
	}
}



// ----------------------------------------------------------------------
// Called after a looping path from a branch point has been detected to
// see if, without the path, we would still count the node as a branch
// point.
// ----------------------------------------------------------------------
void db_graph_check_still_branch_after_loop_detected(pathStep * step, dBGraph * db_graph)
{
	dBNode *node = step->node;
	char tmp_seq[db_graph->kmer_size + 1];
	int forward_n = db_node_edges_count_all_colours(node, forward);
	int reverse_n = db_node_edges_count_all_colours(node, reverse);

	if (DEBUG) {
		binary_kmer_to_seq(element_get_kmer(node), db_graph->kmer_size, tmp_seq);
		printf("First step: label=%c, orientation=%s, kmer=%s\n",
		       binary_nucleotide_to_char(step->label),
		       step->orientation == forward ? "forward" : "reverse",
		       tmp_seq);
	}
	// Subtract one from the count of paths in the orientation of the loop
	if (step->orientation == forward) {
		forward_n--;
	} else {
		reverse_n--;
	}

	// If now there is just one forward and one reverse, then remove the branch flag
	// and set IGNORE_START_NODE, so we know to ignore the paths being built
	if ((forward_n == 1) && (reverse_n == 1)) {
		if (DEBUG) {
			printf("[db_graph_check_still_branch_after_loop_detected] Unsetting branch flags\n");
		}
		db_node_action_unset_flag(node, BRANCH_NODE_FORWARD | BRANCH_NODE_REVERSE);
		db_node_action_set_flag(node, IGNORE_START_NODE);
	}
}

// ----------------------------------------------------------------------
// Walk path from branch node
// ----------------------------------------------------------------------
void db_graph_walk_from_node(dBNode * node, Path * current_path, int orientation, int depth, int max_depth,
                             int max_length, PathArray * patharray, dBGraph * db_graph)
{
	pathStep first_step;
	char node_kmer[db_graph->kmer_size + 1];
	char tmp_seq[db_graph->kmer_size + 1];

	if (BUBBLEDEBUG) {
		int i;
		for (i = 0; i < depth; i++) {
			printf("    ");
        }
		binary_kmer_to_seq(element_get_kmer(node), db_graph->kmer_size, node_kmer);
		printf("Path so far %s Start node %s orientation %s\n", current_path->seq, node_kmer, orientation == forward ? "forward" : "reverse");
	}
    
	// Function to go through each nucleotide, walk the path if an edge exists and store path
	void walk_if_exists(Nucleotide n) {
		if (db_node_edge_exist_any_colour(node, n, orientation)) {
			Path *new_path;
			Path *merged_path;
			dBNode *end_node;
			int end_orientation;
			first_step.node = node;
			first_step.orientation = orientation;
			first_step.label = n;

			new_path = path_new(max_length, db_graph->kmer_size);
			merged_path = path_new(max_length, db_graph->kmer_size);
			if ((!new_path) || (!merged_path)) {
				fprintf(stderr,	"\n[db_graph_walk_from_node] Couldn't get memory new_path.\n\n");
				exit(-1);
			}
			// Get the next bit of path
			db_graph_get_perfect_path_with_first_edge_all_colours(&first_step, &db_node_action_do_nothing, new_path, db_graph);

			// Make a new path
			if (current_path->length == 0) {
				path_copy(merged_path, new_path);
			} else {
				path_copy(merged_path, current_path);
				path_append(merged_path, new_path);
			}
			end_node = new_path->nodes[new_path->length - 1];
			end_orientation = new_path->orientations[new_path->length - 1];

			if (BUBBLEDEBUG) {
				int i;
				for (i = 0; i < depth; i++) {
					printf("    ");
                }
			}
            
			// What to do with the new path?
			// - If path ends at a blunt end, we can go down no deeper
			// - If first node in path is equal to the end node, it's a loop, so we don't add it to the
			//   patharray and we check if we can ignore this node completely.
			// - If we've reached the maximum length, then we look no deeper and add the path to the 
			//   path to the patharray.
			// - If we've reached the maxium depth, then we look no deeper and add the path to the
			//   patharray.
			// - Otherwise, we search deeper...
			if (db_node_is_blunt_end_all_colours(end_node, end_orientation)) {
				if (BUBBLEDEBUG) {
					binary_kmer_to_seq(element_get_kmer(end_node), db_graph->kmer_size, tmp_seq);
					printf("  Path ends at %s (orientation %s) - adding...\n", tmp_seq, end_orientation == forward ? "forward" : "reverse");
				}
				db_graph_check_and_add_path(merged_path, patharray);
			} else if (merged_path->nodes[0] == merged_path->nodes[merged_path->length - 1]) {
				pathStep step;
				if (BUBBLEDEBUG) {
					printf("  Path is complete circle - ignoring...\n");
				}
				db_graph_check_still_branch_after_loop_detected(path_get_step_at_index(0, &step, merged_path), db_graph);
			} else if (path_is_cycle(merged_path)) {
				if (BUBBLEDEBUG) {
					printf("  Path contains cycle.. Quitting and adding...\n");
				}
				db_graph_check_and_add_path(current_path, patharray);
			} else if (current_path->length + new_path->length >= max_length) {
				if (BUBBLEDEBUG) {
					printf("  Maximum path length reached... Quitting and adding...\n");
				}
				db_graph_check_and_add_path(merged_path, patharray);
			} else if ((depth + 1) > max_depth) {
				if (BUBBLEDEBUG) {
					printf("  Max depth reached. Adding path.\n");
				}
				db_graph_check_and_add_path(merged_path, patharray);
			} else {
				if (BUBBLEDEBUG) {
					printf("  New path is %s\n", new_path->seq);
				}
				db_graph_walk_from_node(end_node, merged_path, end_orientation, depth + 1, max_depth, max_length, patharray, db_graph);
			}

			// Destroy paths not needed
			path_destroy(merged_path);
			path_destroy(new_path);
		}
	}

	nucleotide_iterator(&walk_if_exists);

	if (BUBBLEDEBUG) {
		int i;
		for (i = 0; i < depth; i++) {
			printf("    ");
        }
		printf("End node %s\n", node_kmer);
	}
}

// ----------------------------------------------------------------------
// Identify branch points of potential SNPs and indels and write
// sequences to file.
// ----------------------------------------------------------------------
void db_graph_walk_branches(char *filename, int total_max_length, int bubble_max_length, int bubble_max_depth, dBGraph * db_graph)
{
	pathStep end_step;
	PathArray *patharray;
	int orientation;
	int number_walked = 0;
	int contigs_output = 0;
	int id = 0;
	int max_path_array_size = MAX_PATHS_IN_ARRAY;

	//log_and_screen_printf("Looking for branches to walk.\n");

	// Create output files
	db_graph_prepare_output_files(filename);

	// Function to go through hash table and find points marked as branches
	void find_branch_points(dBNode * node) {
		// If node is marked as forward branch, reverse branch or X
		// node, then try and walk and compare...
		if (db_node_check_for_any_flag(node, BRANCH_NODE_FORWARD | BRANCH_NODE_REVERSE | X_NODE))
		{
			Path *initial_path = NULL;
			number_walked++;

			patharray = path_array_new(max_path_array_size);

			// Slightly different behaviour for Y nodes and X nodes
			if (db_node_check_for_any_flag(node, BRANCH_NODE_FORWARD | BRANCH_NODE_REVERSE)) {
				if (DEBUG) {
					printf("[find_branch_points] Got Y node\n");
				}
				// Set orientation (used by walk_if_exists)
				orientation = db_node_check_for_any_flag(node, BRANCH_NODE_REVERSE) ? reverse : forward;

				// Build array of paths
				if (BUBBLEDEBUG) {
					char tmp_seq[db_graph->kmer_size + 1];
					binary_kmer_to_seq(element_get_kmer(node), db_graph->kmer_size, tmp_seq);
					printf("\nWalking from node %s orientation %s\n", tmp_seq, orientation == forward ? "forward" : "reverse");
				}
                
				initial_path = path_new(bubble_max_length, db_graph->kmer_size);
   				db_graph_walk_from_node(node, initial_path,	orientation, 0,	bubble_max_depth, bubble_max_length, patharray, db_graph);

				// Check if the IGNORE_START_NODE flag was set for this node, and if so do nothing with the paths
				if (db_node_check_for_any_flag(node, IGNORE_START_NODE)) {
					if (DEBUG) {
						printf("IGNORE_START_NODE set\n");
					}
					db_node_action_unset_flag(node, IGNORE_START_NODE);
				} else {
					db_graph_walk_display_paths(patharray);
					if (db_graph_compare_paths(patharray, &end_step, db_graph->kmer_size)) {
						id = db_graph_found_matched_paths(patharray, total_max_length, orientation, node, &end_step, filename, db_graph);
						if (id >= 0) {
							contigs_output = id;
                        }
					}
				}
			} else if (db_node_check_for_any_flag(node, X_NODE)) {
				if (BUBBLEDEBUG) {
					char tmp_seq[db_graph->kmer_size + 1];
					binary_kmer_to_seq(element_get_kmer(node), db_graph->kmer_size, tmp_seq);
					printf("\nWalking from X-node %s orientation forward\n", tmp_seq);
				}

				if (DEBUG) {
					printf("[find_branch_points] Got X node\n");
				}
                
				// Look in both directions.. forward first
				if (DEBUG) {
					printf("[find_branch_points] Orientation forward\n");
				}
                
				// Build array of paths forward
				orientation = forward;
				initial_path = path_new(bubble_max_length, db_graph->kmer_size);
				db_graph_walk_from_node(node, initial_path,	orientation, 0,	bubble_max_depth, bubble_max_length, patharray, db_graph);
				if (db_node_check_for_any_flag(node, IGNORE_START_NODE)) {
					db_node_action_unset_flag(node, IGNORE_START_NODE);
				} else {
					db_graph_walk_display_paths(patharray);
					if (db_graph_compare_paths(patharray, &end_step, db_graph->kmer_size)) {
						id = db_graph_found_matched_paths(patharray, total_max_length, orientation, node, &end_step, filename, db_graph);
						if (id >= 0) {
							contigs_output = id;
                        }
					}
				}

				// Destroy paths
				path_destroy(initial_path);
				path_array_destroy(patharray);
				patharray = path_array_new(max_path_array_size);

				if (BUBBLEDEBUG) {
					char tmp_seq[db_graph->kmer_size + 1];
					binary_kmer_to_seq(element_get_kmer(node), db_graph->kmer_size, tmp_seq);
					printf("\nWalking from X-node %s orientation reverse\n", tmp_seq);
				}
                
				// then reverse
				if (DEBUG) {
					printf("[find_branch_points] Orientation reverse\n");
				}
                
				// Build array of paths reverse
				orientation = reverse;
				initial_path = path_new(bubble_max_length, db_graph->kmer_size);
				db_graph_walk_from_node(node, initial_path,	orientation, 0,	bubble_max_depth, bubble_max_length, patharray, db_graph);
				if (db_node_check_for_any_flag(node, IGNORE_START_NODE)) {
					db_node_action_unset_flag(node, IGNORE_START_NODE);
				} else {
					db_graph_walk_display_paths(patharray);
					if (db_graph_compare_paths(patharray, &end_step, db_graph->kmer_size)) {
						id = db_graph_found_matched_paths(patharray, total_max_length, orientation, node, &end_step, filename, db_graph);
						if (id >= 0) {
							contigs_output = id;
                        }
					}
				}
			}
			// Destroy paths
			path_destroy(initial_path);
			path_array_destroy(patharray);
		}
	}

	hash_table_traverse(&find_branch_points, db_graph);

	log_and_screen_printf("Finished walking - %d contigs output.\n", contigs_output);
}

/*
// OLD VERSION
// Identify branch points of potential SNPs and indels and write sequences to file
void db_graph_walk_branches(char * filename, int max_length, dBGraph * db_graph) {
	pathStep first_step;
	pathStep end_step;
	PathArray * patharray;
	int orientation;
	int number_walked = 0;
	int id = 0;

	printf("Looking for branches to walk.\n");
	
	// Create output files
	db_graph_prepare_output_files(filename);
	
	// Function to go through hash table and find points marked as branches
	void find_branch_points(dBNode * node) {
		// Function to go through each nucleotide, walk the path if an edge exists and store path
		void walk_if_exists(Nucleotide n) {
			if (db_node_edge_exist_any_colour(node,n,orientation)) {
				Path * new_path;
				first_step.node = node;
				first_step.orientation = orientation;
				first_step.label = n;
				
				if (DEBUG)
					printf("[find_branch_points] Walking path for %c\n", binary_nucleotide_to_char(n));
				
				if (patharray->number_of_paths > patharray->capacity) {
					fprintf(stderr, "\n[find_branch_points] Too many paths.\n\n");
					exit(-1);
				}
				
				new_path = path_new(max_length, db_graph->kmer_size);
				if (!new_path) {
					fprintf(stderr, "\n[find_branch_points] Couldn't get memory new_path.\n\n");
					exit(-1);						
				}

				path_array_add_path(new_path, patharray);
				path_reset(patharray->paths[patharray->number_of_paths-1]);
				db_graph_get_perfect_path_with_first_edge_all_colours(
					&first_step,
					&db_node_action_do_nothing,
					patharray->paths[patharray->number_of_paths-1],
					db_graph);
			}
		}

		// If node is marked as forward branch, reverse branch or X
		// node, then try and walk and compare...
		if (db_node_check_for_any_flag(node, BRANCH_NODE_FORWARD | BRANCH_NODE_REVERSE | X_NODE)) {
			number_walked++;

			patharray = path_array_new(4);
		
			// Slightly different behaviour for Y nodes and X nodes
			if (db_node_check_for_any_flag(node, BRANCH_NODE_FORWARD | BRANCH_NODE_REVERSE)) {
				if (DEBUG)
					printf("[find_branch_points] Got Y node\n");

				// Set orientation (used by walk_if_exists)
				orientation = db_node_check_for_any_flag(node, BRANCH_NODE_REVERSE) ? reverse:forward;
							
				// Walk the 4 potential nucleotide paths, building an array
				// of paths. Call function to compare these paths and then
				// another separate function to deal with matches...
				nucleotide_iterator(&walk_if_exists);
				db_graph_walk_display_paths(patharray);
				if (db_graph_compare_paths(patharray, &end_step, db_graph->kmer_size)) {					
					id = db_graph_found_matched_paths(patharray, max_length, orientation, node, &end_step, filename, db_graph);
				}
			}
			else if (db_node_check_for_any_flag(node, X_NODE)) {
				if (DEBUG)
					printf("[find_branch_points] Got X node\n");

				// Look in both directions.. forward first
				if (DEBUG)
					printf("[find_branch_points] Orientation forward\n");
				orientation = forward;
				nucleotide_iterator(&walk_if_exists);
				db_graph_walk_display_paths(patharray);
				if (db_graph_compare_paths(patharray, &end_step, db_graph->kmer_size)) {					
					id = db_graph_found_matched_paths(patharray, max_length, orientation, node, &end_step, filename, db_graph);
				}

				// Destroy paths
				path_array_destroy(patharray);
				patharray = path_array_new(4);

				// then reverse
				if (DEBUG)
					printf("[find_branch_points] Orientation reverse\n");
				orientation = reverse;
				nucleotide_iterator(&walk_if_exists);
				db_graph_walk_display_paths(patharray);
				if (db_graph_compare_paths(patharray, &end_step, db_graph->kmer_size)) {					
					id = db_graph_found_matched_paths(patharray, max_length, orientation, node, &end_step, filename, db_graph);
				}
			}

			// Destroy paths
			path_array_destroy(patharray);
		}
	}
	
	hash_table_traverse(&find_branch_points, db_graph);

	printf("Finished walking - %d contigs output.\n", id);
}
*/

// ----------------------------------------------------------------------
// Return point at which two paths diverge. If paths are same up until
// the end of one of them, divergence point will equal shortest length.
// ----------------------------------------------------------------------
void db_graph_find_divergence(Path * path_a, Path * path_b, int *a_ctr, int *b_ctr)
{
	pathStep step_a, step_b;

	while ((*a_ctr < path_a->length) && (*b_ctr < path_b->length)) {
		path_get_step_at_index(*a_ctr, &step_a, path_a);
		path_get_step_at_index(*b_ctr, &step_b, path_b);
		if (!path_step_equals(&step_a, &step_b)) {
			break;
		}
		(*a_ctr)++;
		(*b_ctr)++;
	}
}

// ----------------------------------------------------------------------
// Return point at which two paths converge. If they don't, return -1.
// ----------------------------------------------------------------------
int db_graph_check_for_convergence(Path * path_a, Path * path_b, int *a_ctr, int *b_ctr)
{
	pathStep step_a, step_b;
	int index = -1;
	int k;

	while ((index == -1) && (*a_ctr < path_a->length)) {
		// Get next step from path A
		path_get_step_at_index(*a_ctr, &step_a, path_a);

		// See if it appears in path B
		for (k = *b_ctr; k < path_b->length; k++) {
			path_get_step_at_index(k, &step_b, path_b);
			if (path_step_equals_without_label(&step_a, &step_b)) {
				// Update b_ctr
				*b_ctr = k;

				// The return value (index) is the lowest position this step occurs in the two paths
				index = *a_ctr < k ? *a_ctr : k;

				// Go no further round the loop
				break;
			}
		}

		// Only increment a_ctr if we've not found a convergence point
		if (index == -1) {
			(*a_ctr)++;
		}
	}

	return index;
}

// ----------------------------------------------------------------------
// Count how many paths in the array contain the given step
// ----------------------------------------------------------------------
int db_graph_count_paths_with_step(pathStep * step, PathArray * patharray)
{
	int i;
	int count = 0;

	for (i = 0; i < patharray->number_of_paths; i++) {
		if (path_contains_step(step, patharray->paths[i])) {
			count++;
		}
	}

	return count;
}

// ----------------------------------------------------------------------
// Compare paths to see if they end at the same node
// ----------------------------------------------------------------------
boolean db_graph_compare_paths(PathArray * patharray, pathStep * end, int kmer_size)
{
	int i, j;
	int end_paths = 0;
	int end_index = 0;
	int index = 0;
	int n_paths = 0;
	boolean match_found = false;
	pathStep step;
	int converge;
	int i_ctr;
	int j_ctr;

	// No point comparing if only 1 path
	if (patharray->number_of_paths < 2) {
		return false;
	}
	// Reset match flags
	for (i = 0; i < patharray->number_of_paths; i++) {
		path_action_unset_flag(patharray->paths[i], MATCH_FOUND);
	}

	if (DEBUG) {
		printf("[db_graph_compare_paths] Starting to compare paths...\n");
	}
	// Compare each path against each other 
	for (i = 0; i < patharray->number_of_paths; i++) {
		for (j = i + 1; j < patharray->number_of_paths; j++) {
			if (BUBBLEDEBUG) {
				printf("[db_graph_compare_paths] Comparing path %d with path %d...\n", i, j);
			}
            
			// ALGORITHM
			//
			// Given a branch point and a set of paths leading out of that branch point
			//
			// To find the end point for a given pair of paths:                     
			// 1. Find path step that these paths diverge on.
			// 2. Now find the step that they converge - store the node as the end node for these two paths.
			// 3. There may be other paths that converge here also, so count total number of paths that contain this step. 
			// 4. Keep looking further down the path until either:
			//    - we reach a node that splits in the orientation of travel, or
			//    - we reach the end of the path
			// 5. If, whilst continuing to walk down this path, we find a node where another path joins the paths we already have (2 or more of them), then this is now the end node.
			//              
			// To decide which end point represents the one end point for this branch point:
			// 1. Choose the one which has most paths meeting at it.
			// 2. If two end points have the same number of paths, then choose the closest one to the branch point.

			// Find the first point these paths diverge
			i_ctr = 0;
			j_ctr = 0;
			db_graph_find_divergence(patharray->paths[i], patharray->paths[j], &i_ctr, &j_ctr);

			if (BUBBLEDEBUG) {
				printf("[db_graph_compare_paths] Divergence counters %d and %d\n", i_ctr, j_ctr);
			}
			// Move counters past the divergence point, otherwise this will be picked up as a convergence point too!
			i_ctr++;
			j_ctr++;

			// Now find if they converge
			converge = db_graph_check_for_convergence(patharray->paths[i], patharray->paths[j], &i_ctr, &j_ctr);

			if (BUBBLEDEBUG) {
				printf("[db_graph_compare_paths] Convergence point %d\n", converge);
			}
			// Convergence index will be -1 if they don't converge, otherwise it will be the index where they do
			if (converge != -1) {
				// Count paths at this convergence point
				path_get_step_at_index(i_ctr, &step, patharray->paths[i]);
				n_paths = db_graph_count_paths_with_step(&step, patharray);
				if ((n_paths > end_paths) || ((n_paths == end_paths)  && (converge < end_index))) {
					path_get_step_at_index(i_ctr, end, patharray->paths[i]);
					end_index = converge;
				}
				// Now follow rest of path and see if we get to a point where there are more paths converging
				while ((i_ctr < patharray->paths[i]->length) && (j_ctr < patharray->paths[j]->length)) {
					path_get_step_at_index(i_ctr, &step, patharray->paths[i]);
					n_paths = db_graph_count_paths_with_step(&step, patharray);
					index = i_ctr < j_ctr ? i_ctr : j_ctr;

					if ((n_paths > end_paths) || ((n_paths == end_paths) && (index < end_index))) {
						path_step_assign(end, &step);
						end_paths = n_paths;
						end_index = index;
					}
                    
					// If there is more than one path out of this node, we don't go any further
					if (db_node_edges_count_all_colours (step.node, step.orientation) > 1) {
						break;
					}

					i_ctr++;
					j_ctr++;
				}
			}
		}
	}

	if (end_paths > 0) {
		match_found = true;
	}

	if (BUBBLEDEBUG) {
		if (!match_found) {
			printf("[db_graph_compare_paths] No matches.\n");
		}
	}

	return match_found;
}

/*
// Compare paths to see if they end at the same node
boolean db_graph_compare_paths(PathArray * patharray, pathStep * end, int kmer_size)
{
	int i,j,k;
	int diverge_index;
	int length;
	char tmp_seq[kmer_size+1];
	boolean match_found = false;
	
	if (patharray->number_of_paths < 2) {
		return false;
	}
	
	for (i=0; i<patharray->number_of_paths;  i++) {
		path_action_unset_flag(patharray->paths[i], MATCH_FOUND);
	}
	
	if (DEBUG) {
		printf("[db_graph_compare_paths] Starting to compare paths...\n");
	}
	
	void check_each_node(pathStep * step) {
		if (!match_found) {
			if (path_contains_step(step, patharray->paths[j])) {
				// Store this node as the end node
				path_step_assign(end, step);
				match_found = true;
				//path_action_set_flag(patharray->paths[i], MATCH_FOUND);
				//path_action_set_flag(patharray->paths[j], MATCH_FOUND);
				if (DEBUG) {
					printf("[db_graph_compare_paths] Match found.\n");
					printf("[db_graph_compare_paths] Matching paths %d and %d.\n", i, j);
					printf("[db_graph_compare_paths] Matching kmer is %s\n",
						   binary_kmer_to_seq(element_get_kmer(step->node), patharray->paths[j]->kmer_size, tmp_seq));
				}
			}
		}
	}
	
	// Compare each path against each other to see if it ends at a common point
	for (i=0; i<patharray->number_of_paths; i++) {
		for (j=i+1; j<patharray->number_of_paths; j++) {
			if (DEBUG) {
				printf("[db_graph_compare_paths] Comparing path %d with path %d...\n", i, j);
			}
			
			// Find the first point these paths diverge
			diverge_index = 0;
			length = patharray->paths[i]->length;
			if (patharray->paths[j]->length < length) {
				length = patharray->paths[j]->length;
			}
			
			for (k=0; k<length; k++) {
				pathStep step_a, step_b;
				path_get_step_at_index(k, &step_a, patharray->paths[i]);
				path_get_step_at_index(k, &step_b, patharray->paths[j]);
				if (!path_step_equals(&step_a, &step_b)) {
					diverge_index = k;
					break;
				}
			}
			
			// Go through each step of path i and see if node is in path j                  
			if (!match_found) {
				path_iterator_from_index(diverge_index, &check_each_node, patharray->paths[i]);
			}
		}
	}
	
	if (DEBUG) {
		if (!match_found) {
			printf("[db_graph_compare_paths] No matches.\n");
		}
	}
	
	return match_found;
}
*/

// ----------------------------------------------------------------------
// Having found matched pairs, assemble a complete sequence with prefix
// and suffix and print it to screen.
// ----------------------------------------------------------------------
int db_graph_found_matched_paths(PathArray * patharray, int max_length, Orientation orientation, dBNode * start_node,
                                 pathStep * end_step, char *filename, dBGraph * db_graph)
{
	Path *prefix_path = 0;
	Path *suffix_path = 0;
	Path *merged_path;
	dBNode *end_node = end_step->node;
	pathStep step;
	char *stats_string;
	char type_string[3];
	boolean success;
	boolean unique_path;
	int last_unique_path = -1;
	int i, j;
	int id;
	int steps_same;

	// Remove all path after the end node
	for (i = 0; i < patharray->number_of_paths; i++) {
		if (path_contains_step_without_label
		    (end_step, patharray->paths[i])) {
			if (BUBBLEDEBUG) {
				printf("[db_graph_found_matched_paths] Path %d before=%s ", i, patharray->paths[i]->seq);
			}

			do {
				path_get_last_step(&step, patharray->paths[i]);
				if (path_step_equals_without_label(&step, end_step)) {
					patharray->paths[i]->labels[patharray->paths[i]->length-1] = Undefined;
				} else {
					path_remove_last(patharray->paths[i]);
				}
			}
			while ((patharray->paths[i]->length > 0) && (!path_step_equals_without_label(&step, end_step)));

			// Check the path is still unique
			unique_path = true;
			for (j = 0; j < i; j++) {
				if (paths_equal(patharray->paths[i], patharray->paths[j])) {
					unique_path = false;
					break;
				}
			}

			if (BUBBLEDEBUG) {
				printf(" after=%s\n", patharray->paths[i]->seq);
			}

			if (unique_path) {
				path_action_set_flag(patharray->paths[i], MATCH_FOUND);
				last_unique_path = i;
			} else {
				if (BUBBLEDEBUG) {
					printf("Throwing away non-unique path\n");
				}
			}
		}
	}

	if (last_unique_path == -1) {
		log_and_screen_printf("WARNING: No unique paths for match. Something went wrong...\n");
		return -1;
	}
	// Now, for those paths that we've found share a common end point, check where the common start point is
	// Put the index of divergence into diverge_index.
	pathStep step_a, step_b;
	int diverge_index = 0;
	for (j = 0; j < patharray->paths[last_unique_path]->length; j++) {
		path_get_step_at_index(j, &step_a, patharray->paths[last_unique_path]);
		steps_same = 1;
		for (i = 0; i < patharray->number_of_paths; i++) {
			if (path_check_for_flag(patharray->paths[i], MATCH_FOUND)) {
				if (j < patharray->paths[i]->length) {
					path_get_step_at_index(j, &step_b, patharray->paths[i]);
					if (!path_step_equals(&step_a, &step_b)) {
						steps_same = 0;
					} else {
						diverge_index = j;
					}
				} else {
					steps_same = 0;
				}
			}

			if (!steps_same) {
				break;
			}
		}

		if (!steps_same) {
			break;
		}
	}

	// If the paths don't diverge at the first node, then see if the point where they do diverge is also a branch point.
	// If it is, we can ignore these paths...
	if (diverge_index != 0) {
		if (db_node_edges_count_all_colours(step_a.node, reverse) > 1 || db_node_edges_count_all_colours(step_a.node, forward) > 1) {
			if (BUBBLEDEBUG) {
				log_and_screen_printf("WARNING: Common prefix, so ignoring...\n");
			}
			return -1;
		}
	}
	// Output debugging
	if (DEBUG) {
		char seq[db_graph->kmer_size + 1];
		printf("[db_graph_found_matched_paths] Orientation is %d\n", orientation);
		binary_kmer_to_seq(element_get_kmer(start_node), db_graph->kmer_size, seq);
		printf("[db_graph_found_matched_paths] Start node is: %s\n", seq);
		binary_kmer_to_seq(element_get_kmer(end_node), db_graph->kmer_size, seq);
		printf("[db_graph_found_matched_paths] End node is: %s\n", seq);
	}
	// Get prefix path, if we can, by reversing the path from the start node
	int forward_n = db_node_edges_count_all_colours(start_node, forward);
	int reverse_n = db_node_edges_count_all_colours(start_node, reverse);
	if ((forward_n > 1) && (reverse_n == 1)) {
		// It's a BRANCH_NODE_FORWARD
		prefix_path = db_graph_get_surrounding_path(start_node, reverse, true, max_length, db_graph);
		if (db_node_count_number_of_colours_out(start_node, forward) > 1)
			type_string[0] = 'F';
		else
			type_string[0] = 'f';
	} else if ((forward_n == 1) && (reverse_n > 1)) {
		// It's a BRANCH_NODE_REVERSE
		prefix_path = db_graph_get_surrounding_path(start_node, forward, true, max_length, db_graph);
		if (db_node_count_number_of_colours_out(start_node, reverse) > 1)
			type_string[0] = 'R';
		else
			type_string[0] = 'r';
	} else {
		type_string[0] = 'X';
	}

	// Get suffix path
	forward_n = db_node_edges_count_all_colours(end_node, forward);
	reverse_n = db_node_edges_count_all_colours(end_node, reverse);
	if ((forward_n > 1) && (reverse_n == 1)) {
		// It's a BRANCH_NODE_FORWARD
		suffix_path = db_graph_get_surrounding_path(end_node, reverse, false, max_length, db_graph);
		if (db_node_count_number_of_colours_out(end_node, forward) > 1)
			type_string[1] = 'F';
		else
			type_string[1] = 'f';
	} else if ((forward_n == 1) && (reverse_n > 1)) {
		// It's a BRANCH_NODE_REVERSE
		suffix_path = db_graph_get_surrounding_path(end_node, forward, false, max_length, db_graph);
		if (db_node_count_number_of_colours_out(end_node, reverse) > 1)
			type_string[1] = 'R';
		else
			type_string[1] = 'r';
	} else {
		type_string[1] = 'X';
	}

	type_string[2] = 0;

	// Check first step of prefix path and suffix paths isn't in main path. If they are, get rid of them - it's a loop
	for (i = 0; i < patharray->number_of_paths; i++) {
		if (prefix_path) {
			if (path_contains_step(path_get_step_at_index(0, &step, prefix_path), patharray->paths[i])) {
				path_destroy(prefix_path);
				prefix_path = 0;
			}
		}
		if (suffix_path) {
			if (path_contains_step(path_get_step_at_index(0, &step, suffix_path), patharray->paths[i])) {
				path_destroy(suffix_path);
				suffix_path = 0;
			}
		}
	}

	/*
	   if (!db_node_check_for_any_flag(end_node, X_NODE)) {
	   Orientation o_end = db_node_check_for_any_flag(end_node, BRANCH_NODE_REVERSE)?reverse:forward;
	   suffix_path = db_graph_get_surrounding_path(end_node, opposite_orientation(o_end), false, max_length, db_graph);
	   }
	 */

	if (DEBUG) {
		printf("[db_graph_found_matched_paths] Found %d sequences:\n", patharray->number_of_paths);
	}
	// Build output paths, containing prefix, difference, suffix
	for (i = 0; i < patharray->number_of_paths; i++) {
		if (path_check_for_flag(patharray->paths[i], MATCH_FOUND)) {
			if (DEBUG) {
				printf("[db_graph_found_matched_pairs] MATCH_FOUND set on path %d\n", i);
			}
			// Generate colour coverage stats
			stats_string = db_graph_generate_colour_stats_string(patharray->paths[i], prefix_path, suffix_path, type_string);

			// Merge prefix, path and suffix
			merged_path = path_new(max_length * 3, db_graph->kmer_size);
			if (!merged_path) {
				fprintf(stderr, "\n[db_graph_found_matched_pairs] Couldn't get memory for merged path.\n\n");
				exit(-1);
			}

			merged_path->flags = patharray->paths[i]->flags;

			success = true;

			if (prefix_path) {
				success = path_append(merged_path, prefix_path);
            }

			if (success) {
				success = path_append(merged_path, patharray->paths[i]);
            }
            
			if ((success) && (suffix_path)) {
				success = path_append(merged_path, suffix_path);
            }
            
			if (!success) {
				char *message = "warning:path_append_failed ";
				int j;
				for (j = strlen(stats_string); j >= 0; j--) {
					stats_string[j + strlen(message)] = stats_string[j];
                }
				strncpy(stats_string, message, strlen(message));
			}

			merged_path->header = stats_string;

			if (DEBUG) {
				printf("\n");
				printf("Path %d:      ", i);
				if (prefix_path != 0)
					printf("%s-", prefix_path->seq);
				else
					printf("XXXXX-");

				printf("%s", patharray->paths[i]->seq);

				if (suffix_path != 0)
					printf("-%s\n", suffix_path->seq);
				else
					printf("-XXXXX\n");
				printf("Merged path: %s\n", merged_path->seq);
				printf("\n");
			}
			// Replaces path in paths array with merged path
			if (DEBUG) {
				printf("[db_graph_found_matched_pairs] Destroying old path.\n");
            }

			path_destroy(patharray->paths[i]);
			patharray->paths[i] = merged_path;
		}
	}

	// Clear prefix and suffix paths
	if (suffix_path)
		path_destroy(suffix_path);
	if (prefix_path)
		path_destroy(prefix_path);

	// Output paths to file
	id = db_graph_output_search_paths(filename, patharray);

	// Unset flags... but X_NODE stays...
	db_node_action_unset_flag(end_node, BRANCH_NODE_FORWARD | BRANCH_NODE_REVERSE);

	return id;
}

// ----------------------------------------------------------------------
// Get a prefix or suffix path (flanking).
// ----------------------------------------------------------------------
Path* db_graph_get_surrounding_path(dBNode * node, Orientation orientation, boolean reverse_path, int max_length, dBGraph * db_graph)
{
	Nucleotide n;
	Path *path = path_new(max_length, db_graph->kmer_size);

	if (!path) {
		fprintf(stderr, "\n[db_graph_get_surrounding_path] Can't get memory for new path.\n\n");
		exit(-1);
	}

	if (db_node_has_precisely_one_edge_all_colours(node, orientation, &n)) {
		pathStep first_step;
		first_step.node = node;
		first_step.orientation = orientation;
		first_step.label = n;
		db_graph_get_perfect_path_with_first_edge_all_colours(&first_step, &db_node_action_do_nothing, path, db_graph);

		if (DEBUG)
			printf("[db_graph_get_surrounding_path] Path: %s\n", path->seq);

		if (reverse_path) {
			Path *reversed_path = path_new(max_length, db_graph->kmer_size);
			if (reversed_path) {
				path_reverse(path, reversed_path);
				path_destroy(path);
				path = reversed_path;
			} else {
				fprintf(stderr, "\n[db_graph_get_surrounding_path] Couldn't get memory for reversed path.\n\n");
				exit(-1);
			}

			if (DEBUG)
				printf("[db_graph_get_surrounding_path] Path reversed: %s\n", path->seq);
		}
	} else {
		fprintf(stderr, "[db_graph_get_surrounding_path] Something went wrong - more than one edge in opposite orientation.\n");
		return 0;
	}

	return path;
}

// ----------------------------------------------------------------------
// Display paths to screen
// ----------------------------------------------------------------------
void db_graph_walk_display_paths(PathArray * patharray)
{
	if (BUBBLEDEBUG) {
		printf("Displaying %d paths...\n", patharray->number_of_paths);
	}
#ifndef SHORT_FLAGS
	int i;
	for (i = 0; i < patharray->number_of_paths; i++) {
		db_node_action_unset_flag(patharray->paths[i]->nodes[0], PRINT_FORWARD | PRINT_REVERSE);
		db_node_action_set_flag(patharray->paths[i]->nodes[0], patharray->paths[i]->orientations[0] == forward ? PRINT_FORWARD:PRINT_REVERSE);

		if (BUBBLEDEBUG) {
			printf("  Path: %s\n", patharray->paths[i]->seq);
		}
	}
#endif
	
}

// ----------------------------------------------------------------------
// Create blank files, which we'll append data to later
// ----------------------------------------------------------------------
void db_graph_prepare_output_files(char *filename)
{
	char pathname[4096];
	FILE *fp;

	sprintf(pathname, "%s.fasta", filename);
	if (DEBUG) {
		printf("[db_graph_prepare_output_files] Creating output file: %s\n", pathname);
    }
    
	fp = fopen(pathname, "w");
	fclose(fp);

	sprintf(pathname, "%s.coverage", filename);
	if (DEBUG) {
		printf("[db_graph_prepare_output_files] Creating output file: %s\n", pathname);
    }
    
	fp = fopen(pathname, "w");
	fclose(fp);
}

void db_graph_change_case_at(PathArray * patharray, int index)
{
	char case_change = 'a' - 'A';
	int p;

	// Change case of character at given index for all paths where MATCH_FOUND is set.
	for (p = 0; p < patharray->number_of_paths; p++) {
		if (path_check_for_flag(patharray->paths[p], MATCH_FOUND)) {
			if (index < strlen(patharray->paths[p]->seq)) {
				patharray->paths[p]->seq[index] += case_change;
			}
		}
	}
}

// ----------------------------------------------------------------------
// Make differences between two paths lower case.
// ----------------------------------------------------------------------
void db_graph_make_differences_lower_case(PathArray * patharray)
{
	int p;

	for (p = 0; p < patharray->number_of_paths; p++) {
		if (path_check_for_flag(patharray->paths[p], MATCH_FOUND)) {
			char *length_string = strstr(patharray->paths[p]->header, "pre_length");
			int pre = 0, mid = 0, post = 0;

			if (length_string) {
				sscanf(length_string, "pre_length:%i mid_length:%i post_length:%i", &pre, &mid, &post);
				if (pre > 0) {
					pre = pre - patharray->paths[p]->kmer_size;
				} else {
					mid = mid - patharray->paths[p]->kmer_size;
				}

				if (pre + mid >
				    strlen(patharray->paths[p]->seq)) {
					pre = 0;
					mid = 0;
					post = 0;
				}
			}

			if (mid > 0) {
				int i;
				for (i = pre; i < pre + mid; i++) {
					patharray->paths[p]->seq[i] = tolower(patharray->paths[p]->seq[i]);
				}
			}
		}
	}
}

// ----------------------------------------------------------------------
// Write paths to file.
// ----------------------------------------------------------------------
int db_graph_output_search_paths(char *filename, PathArray * patharray)
{
	static int id_counter = 0;
	char id[1024];
	char pathname[4096];
	FILE *fp_seq;
	FILE *fp_cov;
	int pathcount = 0;
	int i, c;

	// Highlight differences in lower case
	db_graph_make_differences_lower_case(patharray);

	// Open files to append...
	sprintf(pathname, "%s.fasta", filename);
	fp_seq = fopen(pathname, "a");
	if (!fp_seq) {
		fprintf(stderr, "\n[db_graph_output_search_paths] Couldn't open fasta file.\n\n");
		exit(-1);
	}

	sprintf(pathname, "%s.coverage", filename);
	fp_cov = fopen(pathname, "a");
	if (!fp_cov) {
		fprintf(stderr, "\n[db_graph_output_search_paths] Couldn't open coverage file.\n\n");
		fclose(fp_seq);
		exit(-1);
	}
	// Output entries
	for (i = 0; i < patharray->number_of_paths; i++) {
		if (path_check_for_flag(patharray->paths[i], MATCH_FOUND)) {
			// Make id
			sprintf(id, "match_%d_path_%d", id_counter, pathcount++);

			// Output path sequence
			path_to_fasta_colour(patharray->paths[i], fp_seq, id);

			if (BUBBLEDEBUG) {
				printf("Outputting %s\n", patharray->paths[i]->seq);
			}
			// Output coverage
			for (c = 0; c < NUMBER_OF_COLOURS; c++) {
				path_to_coverage_colour(patharray->paths[i], fp_cov, id, c);
			}
		}
	}
	fprintf(fp_seq, "\n");
	fclose(fp_seq);
	fprintf(fp_cov, "\n");
	fclose(fp_cov);

	// Increment id counter
	id_counter++;

	return id_counter;
}

// ----------------------------------------------------------------------
// Generate stats string for FASTA file output.
// ----------------------------------------------------------------------
char* db_graph_generate_colour_stats_string(Path * p, Path * prefix, Path * suffix, char *type)
{
	int coverage[NUMBER_OF_COLOURS];
	char tmpstring[128];
	char *string = malloc(256);
	int i, c;
	int n = 0;
	int prefix_length = 0;
	int mid_length = 0;
	int suffix_length = 0;

	// Set all coverages to zero
	for (c = 0; c < NUMBER_OF_COLOURS; c++) {
		coverage[c] = 0;
    }

	// Go through all nodes, adding up colour coverages.
	// Note: ignore first and last nodes, as they are common.
	for (i = 1; i < (p->length - 1); i++) {
		for (c = 0; c < NUMBER_OF_COLOURS; c++) {
			coverage[c] += p->nodes[i]->coverage[c];
		}
		n++;
	}

	// Prefix, mid and suffix lengths
	if (prefix) {
		prefix_length = prefix->length - 1;
	}
	prefix_length += p->kmer_size;
	mid_length = n;
	if (suffix) {
		suffix_length = suffix->length;
	} else {
		suffix_length = 1;
	}

	// Output the type of node.
	sprintf(string, "type:%s", type);

	// Output lengths of prefix, bubble, suffix
	sprintf(tmpstring, " pre_length:%i mid_length:%i post_length:%i", prefix_length, mid_length, suffix_length);
	strcat(string, tmpstring);

	// Now work out averages
	for (c = 0; c < NUMBER_OF_COLOURS; c++) {
		double av = (double)coverage[c] / (double)n;
		sprintf(tmpstring, " c%d_average_coverage:%5.2f", c, av);
		strcat(string, tmpstring);
	}

	return string;
}
