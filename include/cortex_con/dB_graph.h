/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo 
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
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
 
#ifndef DB_GRAPH_H_
#define DB_GRAPH_H_

#include <open_hash/hash_table.h>
#include <stdio.h>
#include <structs.h>
#include <path.h>

typedef HashTable dBGraph;

#ifndef MAX_STACKED_FUNCTIONS
#define MAX_STACKED_FUNCTIONS 5
#endif

typedef struct {
	dBNode * node;
	Path * path_so_far;
	Orientation orientation;
	int depth;
} NodesToWalkItem;

typedef struct {
	int max_size;
	int number_of_items;
	int head;
	int tail;
	NodesToWalkItem** items;
} NodesToWalkQueue;

typedef struct{
    void (*callback[MAX_STACKED_FUNCTIONS])(Path * path, void * args);
    void *  args[MAX_STACKED_FUNCTIONS];
    int used;
} PathCallbackArray;

typedef struct{
    void (*callback[MAX_STACKED_FUNCTIONS])(pathStep * step, void * args);
    void *  args[MAX_STACKED_FUNCTIONS];
    int used;
} PathStepActionCallbackArray;

typedef struct{
    void (*callback[MAX_STACKED_FUNCTIONS])(dBNode * node, void * args);
    void *  args[MAX_STACKED_FUNCTIONS];
    int used;
} NodeActionCallbackArray;


typedef struct {
	pathStep *  (*get_starting_step)(pathStep * step, dBGraph * db_graph );//Used to find the first  step for a given random node.
	boolean (*continue_backwards)(Path * path, dBGraph * db_graph);
	void (* post_step_action)(pathStep * step);
	void (* pre_step_action)(pathStep * step);
	pathStep * (*get_next_step)(pathStep * current_step, pathStep * next_step, pathStep * reverse_step, dBGraph * db_graph);
	void (* step_action)(pathStep * step);
	boolean (*continue_traversing)(pathStep * current_step, pathStep * next_step, pathStep * reverse_step, Path * path,  dBGraph * db_graph);
	void (* output_callback)(Path * path);
	PathArray * buffers;
    
    PathStepActionCallbackArray step_actions;
    PathCallbackArray   path_callbacks;
    NodeActionCallbackArray node_callbacks;
    dBGraph * db_graph; 
} WalkingFunctions;

boolean db_graph_remove_path_callback(WalkingFunctions * wf, void * funct);
void db_graph_add_path_callback(WalkingFunctions * wf, void (*path_callback)(Path * path));
void db_graph_add_step_action(WalkingFunctions * wf, void (*step_action)(pathStep * ps));
void db_graph_add_node_action(WalkingFunctions * wf, void (*node_action)(dBNode * node));

void db_graph_add_path_callback_with_args(WalkingFunctions * wf, void (*path_callback)(Path * path, void * arg), void * args);
void db_graph_add_step_action_with_args(WalkingFunctions * wf, void (*step_action)(pathStep * ps, void * arg), void * args);
void db_graph_add_node_action_with_args(WalkingFunctions * wf, void (*node_action)(dBNode * node, void * arg), void * args);


int db_graph_generic_walk(pathStep * first_step, Path * path, WalkingFunctions * functions, dBGraph * db_graph);

pathStep * db_graph_get_next_step(pathStep * current_step,
		pathStep * next_step, pathStep * rev_step, dBGraph * db_graph);

int db_graph_get_perfect_path_with_first_edge_all_colours(pathStep * first_step,
														  void(*node_action)(dBNode * node),
														  Path * path,
														  dBGraph * db_graph);

int db_graph_get_perfect_path_with_first_edge_long(dBNode * node,
		Orientation orientation, int limit, Nucleotide fst_nucleotide,
		void(*node_action)(dBNode * node), dBNode * * path_nodes,
		Orientation * path_orientations, Nucleotide * path_labels, char * seq,
		double * avg_coverage, int * min_coverage, int * max_coverage,
		boolean * is_cycle, dBGraph * db_graph);
void db_graph_print_status(dBGraph * db_graph);
int db_graph_clip_tip(dBGraph * db_graph);

void db_graph_print_supernode_from_begining(char * filename, int max_length,
		boolean with_coverages, dBGraph * db_graph);
int db_graph_get_perfect_path_with_first_edge_multipath(dBNode * node,
		Orientation orientation, Nucleotide fst_nucleotide, void(*node_action)(
				dBNode * node), Path * path, dBGraph * db_graph);

void db_graph_path_calculate_coverage(dBNode * * path_nodes,
		double * avg_coverage, int * min_coverage, int * max_coverage,
		int length, boolean notHead);

int db_graph_db_node_prune_edges_with_single_coverage(dBNode * node,
		int coverage, void(*node_action)(dBNode * node), dBGraph * db_graph);

long long db_graph_prune_low_coverage_edges(int coverage_threshold,
		dBGraph * db_graph);

int db_graph_get_perfect_path(dBNode * node, Orientation orientation,
		void(*node_action)(dBNode * node), dBGraph * db_graph, Path * path);

boolean db_graph_path_splits(dBNode * node, Orientation o, dBGraph * db_graph,
		int limit);
		
double db_graph_get_average_coverage(dBGraph * db_graph);

Orientation db_graph_get_reverse_edge_between(dBNode * src, dBNode * tgt,
		dBGraph * db_graph);
/*boolean db_graph_detect_bubble(dBNode * node, Orientation orientation,
		int limit, void(*node_action)(dBNode * node), int * length1,
		dBNode ** path_nodes1, Orientation * path_orientations1,
		Nucleotide * path_labels1, char * seq1, double * avg_coverage1,
		int * min_coverage1, int * max_coverage1, int * length2,
		dBNode ** path_nodes2, Orientation * path_orientations2,
		Nucleotide * path_labels2, char * seq2, double * avg_coverage2,
		int * min_coverage2, int * max_coverage2, dBGraph * db_graph);
*/

boolean db_graph_detect_bubble(dBNode * node, Orientation orientation,
							   PathArray * pa, dBGraph * db_graph);
							   
boolean db_graph_db_node_smooth_bubble(dBNode * node, Orientation orientation,
		int limit, int coverage_limit, void(*node_action)(dBNode * node),
		dBGraph * db_graph);

boolean db_graph_db_node_prune_low_coverage(dBNode * node, int coverage,
		void(*node_action)(dBNode * node), dBGraph * db_graph);

// limit is the max length
// min_coverage, max_coverage and avg_coveragte refer to the internal nodes
int db_graph_supernode(dBNode * node, void(*node_action)(dBNode * node),
		Path * path, dBGraph * db_graph);

void printNode(dBNode * dbn, short int kmerSize);

int db_graph_db_node_clip_tip_with_orientation(dBNode * node,
		Orientation orientation, int limit, void(*node_action)(dBNode * node),
		dBGraph * db_graph);
		
int db_graph_db_node_clip_tip(dBNode * node, int limit, void(*node_action)(
		dBNode * node), dBGraph * db_graph);

boolean db_graph_is_condition_true_for_all_nodes_in_supernode(dBNode * node,
		int limit, boolean(*condition_for_all_nodes)(dBNode * node),
		void(*node_action)(dBNode * node), dBNode * * path_nodes,
		Orientation * path_orientations, Nucleotide * path_labels,
		int* path_length, char * string, double * avg_coverage, int * min,
		int * max, boolean * is_cycle, dBGraph * db_graph);

boolean db_graph_is_condition_true_for_at_least_one_node_in_supernode(
		dBNode * node, int limit, boolean(*condition_for_all_nodes)(
				dBNode * node), void(*node_action)(dBNode * node),
		dBNode * * path_nodes, Orientation * path_orientations,
		Nucleotide * path_labels, int* path_length, char * string,
		double * avg_coverage, int * min, int * max, boolean * is_cycle,
		dBGraph * db_graph);

boolean
db_graph_is_condition_true_for_start_and_end_but_not_all_nodes_in_supernode(
		dBNode * node, int limit, boolean(*condition_for_all_nodes)(
				dBNode * node), void(*node_action)(dBNode * node),
		int min_start, int min_end, int min_diff, dBNode * * path_nodes,
		Orientation * path_orientations, Nucleotide * path_labels,
		int * path_length, char * string, double * avg_coverage, int * min,
		int * max, boolean * is_cycle, dBGraph * db_graph);

void
db_graph_print_supernodes_where_condition_is_true_for_all_nodes_in_supernode(
		dBGraph * db_graph, boolean(*condition)(dBNode * node),
		int min_covg_required, FILE* fout, boolean is_for_testing,
		char** for_test_array_of_supernodes, int* for_test_index);

void
		db_graph_print_supernodes_where_condition_is_true_for_at_least_one_node_in_supernode(
				dBGraph * db_graph, boolean(*condition)(dBNode * node),
				int min_covg_required, FILE* fout, boolean is_for_testing,
				char** for_test_array_of_supernodes, int* for_test_index);

void
		db_graph_print_supernodes_where_condition_is_true_at_start_and_end_but_not_all_nodes_in_supernode(
				dBGraph * db_graph, boolean(*condition)(dBNode * node),
				int min_covg_required, int min_start, int min_end,
				int min_diff, FILE* fout, boolean is_for_testing,
				char** for_test_array_of_supernodes, int* for_test_index);

void db_graph_print_supernodes(char * filename, int max_length,
		boolean with_coverages, dBGraph * db_graph);

void db_graph_print_coverage(FILE * f, dBGraph * db_graph);

void db_graph_print_kmer_coverage(FILE * out, dBGraph * db_graph);

int db_graph_clip_tips(int threshold, dBGraph * db_graph); 

//int db_graph_remove_low_coverage_nodes(int coverage, dBGraph * db_graph);

void db_graph_dump_binary(char * filename, boolean(*condition)(dBNode * node),
		dBGraph * db_graph);

void db_graph_dump_binary_by_colour(char * filename, boolean(*condition)(dBNode * node),
									short colour, dBGraph * db_graph);

void db_graph_smooth_bubbles(int coverage, int limit, dBGraph * db_graph);

void db_graph_detect_vars(int max_length, dBGraph * db_graph);

void db_graph_traverse_with_array(void(*f)(HashTable*, Element *, int**, int),
		HashTable * hash_table, int** array, int length_of_array);
int db_graph_get_N50_of_supernodes(dBGraph* db_graph);

int int_cmp(const void *a, const void *b);

dBNode * db_graph_get_first_from_Y(dBNode * n, Orientation o, int limit,
		dBGraph * db_graph, boolean rev);

pathStep* db_graph_get_first_node_in_supernode_containing_given_node(
		pathStep* step, dBGraph* db_graph);

dBNode* db_graph_get_next_node_in_supernode(dBNode* node,
		Orientation orientation, Orientation* next_orientation,
		dBGraph* db_graph);

void db_graph_get_supernode_length_marking_it_as_visited(dBGraph* db_graph,
		Element* node, int** array_of_supernode_lengths, int length_of_array);

int db_graph_db_node_clip_tip(dBNode * node, int limit,
			      void (*node_action)(dBNode * node),
			      dBGraph * db_graph);
			      
dBNode * db_graph_get_next_node(dBNode * current_node, Orientation current_orientation, 
			       Orientation * next_orientation,
			       Nucleotide edge, Nucleotide * reverse_edge,dBGraph * db_graph);
			       
void db_graph_write_graphviz_file(char * filename, dBGraph * db_graph);

int db_graph_get_neighbouring_nodes_all_colours(dBNode* start, Orientation orientation, dBNode* neighbours[], Nucleotide labels[], dBGraph * db_graph);

void db_graph_cleanup_graph(dBGraph * db_graph);

void db_graph_reset_flags(dBGraph * db_graph);

Nucleotide db_graph_get_best_next_step_nucleotide(dBNode * from, dBNode * previous,  Orientation orientation, dBGraph * db_graph );
void db_graph_calculate_stats(dBGraph * db_graph);

double db_graph_get_average_coverage_by_colour(short colour, dBGraph * db_graph);

#endif /* DB_GRAPH_H_ */
