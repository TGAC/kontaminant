//
//  mark_pair.h
//  Untitled
//
//  Created by Ricardo Ramirez-Gonzalez (TGAC) on 03/11/2011.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef mark_pair_h
#define mark_pair_h

typedef struct {
    uint32_t from;
    uint32_t to;
    uint32_t coverage;
}PathConnection;

typedef struct {
    dBGraph * graph;
    long long count;
    long long fail_count;
}mark_args;

typedef unsigned short mp_cov_counter ;
#define MAX_MP_COV USHRT_MAX
#define MAX_READ_LIBRARIES 2
//Should we make this dynamic? if so we need to get some thread safe mechanism to pass buffers to avoid mallocs 
#define MAX_BRANCHES_TO_EXPLORE 4 

typedef struct {
    mp_cov_counter * links;
    long long number_of_paths;
    long long distance_sum;//To try to guess the distance in the library. 
    long long used_reads;
    long long unused_reads;
    long long connected_paths;//To keep a count and evaluate how "sparse" is the matrix. 
    long long overflows;
    int expected_distance;
    int max_distance;
    int min_distance;
    
} ReadConnectionGraph;


void mark_pair_load_pair_data(char* filename, char fastq_ascii_offset,dBGraph * db_graph);

void mark_pair_search_enrich_graph(int library, char * f1, char * f2, char fastq_ascii_offset, dBGraph * db_graph);

long long mark_pair_mark_contiguous_perfect_paths(dBGraph * db_graph);

//mp_cov_counter mark_pair_connections_between(int library, uint32_t p1, uint32_t p2, dBGraph * db_graph);
mp_cov_counter mark_pair_connections_between(uint32_t p1, uint32_t p2, ReadConnectionGraph * rcg);
void mark_pair_print_paths(char *filename, int max_length, boolean with_coverages, boolean with_viz, dBGraph * db_graph);

void mark_pair_validate_file(char * filename);

#endif
