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
 
#ifndef PATH_H_
#define PATH_H_


#ifndef MAX_PATH_LENGTH
#define MAX_PATH_LENGTH 1000000
#endif

#ifndef MAX_PATH_BUFFERS
#define MAX_PATH_BUFFERS 2000
#endif

#ifndef STARTING_PATH_BUFFERS
#define STARTING_PATH_BUFFERS 20
#endif

#define PATH_MAX_DOUBLE_Y_COMPLEXITY 8;

#define PATH_FASTA_LINE 80
#define PATH_MAX_IN_NODES 100
#define PATH_IN_NODES_CAPACITY_INCREASE 10

#define PATH_DEFAULT_SUPERNODES 100
#define PATH_SUPERNODES_CAPACITY_INCREASE 10

#define PATH_FLAG_STOP_BLUNT_END        1
#define PATH_FLAG_CONVERGING_PATHS      2
#define PATH_FLAG_DIVERGING_PATHS       4
#define PATH_FLAG_IS_DOUBLE_Y           8
#define PATH_FLAG_IS_CYCLE              16
#define PATH_FLAG_LONGER_THAN_BUFFER    32
#define PATH_TOO_COMPEX                 64

//Block of flags for the path (may be useful)
#define PRINT_FIRST          	(1 <<  0) //x0000 0001
#define IS_CYCLE			 	(1 <<  1) //x0000 0002
#define NEW_PATH 			 	(1 <<  2) //x0000 0004
#define REVERSE_PATH 		 	(1 <<  3) //x0000 0008
#define EMPTY_PATH 			 	(1 <<  4) //x0000 0010
#define ERROR_PATH			 	(1 <<  5) //x0000 0020
#define FIND_AGAIN 			 	(1 <<  6) //x0000 0040
#define MATCH_FOUND			 	(1 <<  7) //x0000 0080
#define BRANCHED_REVERSE	 	(1 <<  9) //x0000 0100
#define GREEDY_IN_NODE_SEARCH	(1 << 10) //x0000 0200
//Functions

// Block of flags for path step
#define PATH_STEP_VISITED_A    (1 <<  0) //x0000 0001
#define PATH_STEP_VISITED_C    (1 <<  1) //x0000 0002
#define PATH_STEP_VISITED_T    (1 <<  2) //x0000 0004
#define PATH_STEP_VISITED_G    (1 <<  3) //x0000 0008
#define PATH_STEP_MASK_VISITED (~15)     //x1111 1110
#define PATH_STEP_ALL_VISITED  ( 15)     //x0000 000F
#define PRINT_LABEL_AS_N       (1 <<  4) //x0000 0010
#define PRINT_LABEL_LOWERCASE  (1 <<  5) //x0000 0020


typedef enum  {NONE = 0,  FIRST = 1, LAST = 2 }PathEnd ;

typedef struct {
	long long id;
	dBNode * * nodes;
	Orientation * orientations;
	Nucleotide * labels;
    Flags * step_flags;
	int * in_nodes;
	char * seq;
	char * header;
	
	int length;
	int max_length;
    int max_virtual_length; //A soft limit to be used when you want a limit smaller than the buffer size. 
    int new_nodes;
	Flags flags;
	int in_nodes_count;
	int in_nodes_capacity;
    int out_nodes_count;
	short kmer_size;
	boolean used;
	
	Flags stop_reasons_first;
	Flags stop_reasons_last;
    
#ifdef ENABLE_MARK_PAIR
    uint32_t * supernodes;
    uint32_t * supernodes_count;
    uint32_t last_supernode;
    int supernodes_capacity;
    int supernodes_total_count;
#endif
    
} Path;

typedef struct {
	dBNode * node;
	Orientation orientation;
	Nucleotide label;
    Flags flags;//This should be used as read only, it reflects the status of the flag when the step was queried from the path. Whe a step path is added, The flags are added to whatever flags set internally would be used. 
    Path * path;//Pointer to the path to which the step belongs, if any. If NULL, it doesnt really matters. 
} pathStep;

typedef struct{
	pathStep step[4];
	short count;
} nextSteps;



typedef struct{
	Path  ** paths;
	int number_of_paths;
	int capacity;
    short kmer_size;
#ifdef  THREADS
    pthread_mutex_t mutex;
#endif
} PathArray;



typedef struct{
	uint64_t blunt_ends;
	uint64_t converging_paths;
	uint64_t diverging_paths;
	uint64_t is_double_y;//used to count how many paths are double y
	uint64_t is_cycle;
	uint64_t longer_than_buffer;
	 
	uint64_t minimum_double_y; //count how many double Y were marked
	uint64_t total_double_y_lenght; //counts the total lenght of the double Ys. We divide later by the minimum_double_y to get the average. However
}PathCounts;

void path_counts_reset(PathCounts * pc);

void path_counts_add(Path * p, PathCounts * pc);

void path_counts_print_and_log(PathCounts * pc);

void path_do_nothing(Path * p);

Path * path_new(int max_length, short kmer_size);

void path_destroy(Path * path);

void path_copy(Path * to, Path * from);

void path_increase_id(Path * path);

void path_to_fasta(Path * path, FILE * fout);

void path_to_fasta_debug(Path * path, FILE * fout);

void path_to_fasta_colour(Path * path, FILE * fout, char * id);

void path_get_statistics(double * avg_coverage, int * min_coverage,
		int * max_coverage, Path * path);

int path_get_nodes_count(Path * path);

int path_get_edges_count(Path * path);

void path_reset(Path * path);

boolean path_contains_step(pathStep * step, Path * path);

boolean path_contains_step_without_label(pathStep * step, Path * path);

boolean path_contains_step_from_index(int index, pathStep * step, Path * path);

boolean path_is_first_equals_to_last(Path * path);

boolean path_add_node(pathStep * step, Path * path);

boolean path_has_in_step(pathStep * step, Path * path);

void path_set_limit(int limit, Path * p);

int path_get_limit(Path *p);

dBNode * path_last_node(Path * path);

boolean path_is_singleton(int length, Path * path);

boolean path_has_space(Path * path);

pathStep * path_get_last_step_reverse(pathStep * step, Path * path);

pathStep * path_get_step_reverse(pathStep * step, Path * path, int index);

pathStep * path_get_last_step(pathStep * ps, Path * path);

int path_index_of_step(pathStep * step, Path * path);

int path_get_length(Path * path);


PathArray * path_split_in_perfect_paths(Path * p);

Orientation path_last_orientation(Path * path);
/**
 * This method is the one to be called to call the squence. 
 * The grace  of it is that it will remove the last char of 
 * seq whenever the sequence is cyclic when the first and 
 * last steps are the same. in this way we can have 
 * cycles resolved by, for example, read pairs in the path. 
 * 
 */
char * path_get_seq(char * tmp, Path * path);

void path_remove_last(Path * path);

void path_to_coverage(Path * path, FILE * fout);

void path_to_coverage_colour(Path * path, FILE * fout, char * id, short colour);

void path_modfy_last_label(Nucleotide n, Path * p);

/**
 * If the walking implementation supports it, it tells what is the percentage of new
 * nodes. Useful to avoid printing paths that had been printed already. 
 */
int path_percentage_new_nodes(Path * path);

void path_iterator_from_index(int index, void (*step_action)(pathStep * step),Path * path);

void path_iterator(void (*step_action)(pathStep * step),Path * path);

void path_iterator_with_args(void (*step_action) (pathStep * , void *),void * args,  Path * path);

void path_inner_iterator_with_args(void (*step_action) (pathStep * step, void*),void * args,  Path * path);

void path_iterator_reverse(void (*step_action)(pathStep * step),Path * path);

void path_iterator_with_index(void (*step_action) (int index, pathStep * step), Path * path);

void path_inner_iterator(void (*step_action) (pathStep * step), Path * path);

int path_step_compare(pathStep * a, pathStep * b);

void path_get_statistics(double *avg_coverage, int *min_coverage, int *max_coverage, Path * path);

boolean path_is_empty(Path * path);

boolean path_is_cycle(Path * path);

boolean path_to_retry(Path * path);

boolean path_is_blunt(Orientation o, Path * p);

int path_get_index_of_last_in_node(Path * p);

int path_get_first_in_node_after(int pos, Path * path);

boolean path_is_repetitive(double graph_cov, Path * p);

void path_add_stop_reason(PathEnd o, Flags f, Path * p);

boolean path_has_stop_reason(PathEnd o, Flags f, Path * p);

boolean path_has_any_stop_reason(PathEnd o, Flags f, Path * path);

void path_clean_stop_reason(Path * p);

Nucleotide path_last_nucleotide(Path * path);

void path_free_buffer_path(Path * path);

void path_graphviz_open_header(FILE * f);

void path_graphviz_line(FILE * f, Path * p);

void path_graphviz_close_header(FILE * f);

//FLAGS FUCNTIONS
void path_action_clear_flags(Path * node);

void path_action_set_flag(Path * node, Flags f);

void path_action_unset_flag(Path * node, Flags f);

Flags path_get_flags(Path * node, Flags f);

boolean path_check_for_flag(Path * node, Flags flag);

boolean path_check_for_any_flag(Path * node, Flags flag);

void path_reverse(Path * source, Path * destination);

boolean path_append(Path * destination, Path * source);

boolean paths_equal(Path * path_a, Path * path_b);

//STEP functions

void path_step_assign(pathStep*to, pathStep*from);

boolean path_step_equals(pathStep * step, pathStep * other);

boolean path_step_equals_without_label(pathStep * step, pathStep * other);

boolean path_step_has_unvisited_edge_all_colours(pathStep * step);

Nucleotide path_step_get_unvisited_edge_all_colours(pathStep * step);

boolean unlabelled_path_step_equals(pathStep * step, pathStep * other);

//Searches for a node in the path, regardles of the next node. 
boolean path_contains_step_from_index_without_label(int index, pathStep * step, Path * path);

//Adds the node and the index where it is, of a node in the path which contains more than one incoming edge
void path_add_in_step(pathStep * step, Path * path);


//returns how many free spaces the path has to allocate more nodes. 
int path_free_spaces(Path * path);
int path_index_of_last_in_node(Path * path);
pathStep * path_get_step_at_index(int index, pathStep * step, Path * path);

void path_step_print(pathStep * step, int kmer_size, FILE * f);

//PathArray functions
PathArray * path_array_new(short number_of_paths);

void path_array_destroy(PathArray * pa);


void path_array_destroy_buffers();

void path_array_destroy_struct(PathArray ** pa);

int path_array_get_number_of_paths(PathArray *);

boolean path_array_add_path(Path * p, PathArray *pa);

void path_array_merge(PathArray ** from, PathArray * to);

Path * path_array_get(int path, PathArray *pa);
//Buffer functions
void path_array_initialise_buffers(short kmer_size);

void path_array_to_fasta(FILE * f, PathArray * pa);

 /**
  * This assumes that the program is always using the same kmer size! 
  */ 
Path * path_get_buffer_path();

void path_free_buffer_path(Path * path);

Flags path_last_flags(Path * path);
//Buffered array functions
void path_array_free_from_buffer(PathArray * pa);

PathArray * path_array_get_from_buffer_with_size(short size);

void path_array_destroy_buffers();

void path_step_mark_as_uncertain(int i, Path * path, boolean as_n);

boolean is_step_marked_as_uncertain(int i, Path * path);

void path_mark_as_visited(Path* path);

void path_pairs_to_fasta(PathArray* pa, int distances[], FILE* fout);
Path *path_get_buffer_path();Path *path_get_buffer_path();

#ifdef SOLID
NucleotideBaseSpace path_get_first_base(Path * p);

void path_to_base_space_fasta(Path * path, FILE * fout);

#endif

#ifdef ENABLE_MARK_PAIR

void path_add_supernode(uint32_t id,  Path * path);

void path_remove_supernode(uint32_t id, Path * path);

#endif

#endif /* PATH_H_ */

