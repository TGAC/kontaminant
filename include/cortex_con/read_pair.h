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
 
 
#ifndef READ_PAIR_H_
#define READ_PAIR_H_

/**
 * 
 * ReadPairDescriptor contains the description of a read pair
 * distance	 .-	The distance between the two reads
 * tolerance .-	Half of the size of the window to search. The best
 * 			  	match for a kmer i will be in the range 
 * 			  	( i + distance - tolerance, i + distance + tolerance)
 * first	 .-	The first index in the element->signature
 * second	 .-	The second index in the element->signature
 * colour	 .-	Colour of the read pair. 
 * score_path.-	Function that gives a score on how "fit" is the extending 
 * 				path. the return value is between 0 and 1, to be able to
 * 				multiply results of different functions on different 
 * 				descriptors. The last parameter, is a recursive pointer, 
 * 				since the information in the descriptor may be usefull to 
 * 				the scoring function. 
 * 				
 */
typedef struct read_pair_descriptor{
    long long insert_distance_sum;
    long long total_depth;
    long long used_pairs;
	int distance;
	int tolerance;
    int allowed_misses;
	int first;
	int second;
	int colour;
    int insert_size;
	int max_read_length;
    int max_pair_coverage;
    int min_pair_coverage; 
    int depth;
	double (* score_paths)(Path * current, Path * extending, struct read_pair_descriptor *rpda );
    unsigned char * connections;
	BinaryTree * jumps;
	boolean mate;
} ReadPairDescriptor;


/**
 * Array with all the mappings for the read pairs. 
 */ 
typedef struct read_pair_descriptor_array{
	ReadPairDescriptor * * pair;
	
    int number_of_pairs;
	int capacity;
    int walk_distance;
    int maximum_coverage;
    int minimum_bits;
    int minimum_start_length;
    int maximum_paths;
    int minimum_kmers;
    int minimum_jump_insert_size;
    boolean print_stack_as_n;
}ReadPairDescriptorArray;

typedef struct read_pair_jump{
	pathStep first;
	pathStep second;
	int distance;
	boolean valid;
}ReadPairJump;

typedef struct {
    int index;
    int pair;
    ReadPairSignature rps;
}sign_args;

typedef struct {
    dBGraph * graph;
    long long count;
    long long fail_count;
    int max_distance;
}mark_args;

ReadPairJump * read_pair_jump_new();

void read_pair_jump_destroy(ReadPairJump ** rpj);

void read_pair_jump_add_to_descriptor(pathStep * a, pathStep * b, int distance, ReadPairDescriptor * rpd);

int read_pair_jump_compare(void * a, void * b);


ReadPairDescriptorArray * new_read_pair_descriptor_array(char capacity, int walk_distance, int max_coverage,int min_bits, int start_length, int max_paths, int min_kmers, boolean stack_as_n);

void read_pair_descriptor_init_connections();

void destroy_read_pair_descriptor_array(ReadPairDescriptorArray ** rpda);

void read_pair_descriptor_array_add_pair(char first, char second, int colour, int insert, int read_length, int tolerance, int misses, int  min_pair_coverage, int max_pair_coverage,  int depth,  double (* score_paths)(Path * current, Path * extending, ReadPairDescriptor *rpda ),
                                         ReadPairDescriptorArray * rpda);

void read_pair_print_paths(char *filename, int max_length, boolean with_coverages, boolean with_viz, ReadPairDescriptorArray * rpda,dBGraph * db_graph);

void read_pair_print_paths_single_walk(char *filename, int max_length, boolean with_coverages, boolean with_viz, ReadPairDescriptorArray * rpda,dBGraph * db_graph);

int read_pair_get_minimum_length(ReadPairDescriptorArray *);

void read_pair_mark_contiguous_perfect_paths(int max_distance, dBGraph * db_graph);

double paired_signatures_score_paths(Path * left, Path * right, ReadPairDescriptor * rpd);

double simple_score_paths(Path * left, Path * right, ReadPairDescriptor * rpd);

void read_pair_enrich_graph(FILE * f1, FILE * f2, ReadPairDescriptor * rpd,int fastq_ascii_offset, dBGraph * db_graph);

double simple_iterated_score_paths(Path * left, Path * right, ReadPairDescriptor * rpd);

int read_pair_count_bits(ReadPairSignature n);

int read_pair_get_maximum_insert_size(ReadPairDescriptorArray * rpda);
int read_pair_count_common_bits(ReadPairSignature left, ReadPairSignature right);

void read_pair_search_enrich_graph(FILE * f1, FILE * f2, ReadPairDescriptor * rpd,int fastq_ascii_offset, dBGraph * db_graph);


#endif 

