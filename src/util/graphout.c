/*
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

/*----------------------------------------------------------------------*
 * File:    graphout.c                                                  *
 * Purpose: Output a graph file for a portion of a CTX file (or set of) *
 * Author:  Richard Leggett                                             *
 * History: 12-Aug-10: RML: Created                                     *
 *----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <execinfo.h>
#include <signal.h>
#include <ctype.h>
#include <limits.h>
#include <stdint.h>
#include "global.h"
#include "binary_kmer.h"
#include "flags.h"
#include "element.h"
#include "seq.h"
#include "open_hash/hash_table.h"
#include "file_reader.h"
#include "dB_graph.h"
#include "logger.h"
#include "graph_tools.h"
#include "graph_formats.h"
#include "node_queue.h"
#include "coverage_walk.h"
#include "perfect_path.h"

typedef struct {
    dBNode* seed_node;
    int graph_size;
} SubGraphInfo;

/*----------------------------------------------------------------------*
 * Constants                                                            *
 *----------------------------------------------------------------------*/
#define MAX_SEARCH_KMERS 1000
#define MAX_PATHS 100
#define ANALYSIS_QUEUE_SIZE 10000000
#define MAX_SEEDS 100000000
//#define MAX_SEEDS 1000000

/*----------------------------------------------------------------------*
 * Global variables                                                     *
 *----------------------------------------------------------------------*/
char* file_of_filenames = 0;
char* base_filename = 0;
char* debug_filename = 0;
char* analysis_filename = 0;
char* consensus_contigs_filename = 0;
char* search_kmers[MAX_SEARCH_KMERS];
boolean output_ctx = false;
boolean output_ubigraph = false;
boolean output_biolayout = false;
boolean output_gexf = false;
boolean output_gml = false;
boolean output_graphml = false;
boolean output_graphviz = false;
boolean analyse_graph = false;
int kmer_size = 5;
int hash_key_bits = 17;
int bucket_size = 60;
int number_of_search_kmers = 0;
int min_subgraph_kmers = 0;
int clean_graph = 0;
GraphToolsState gt_state;
GraphToolsOptions gt_options;

/*----------------------------------------------------------------------*
 * Function: clean_graph                                                *
 * Purpose:  Remove arcs to non-existent kmers.                         *
 * Params:   None                                                       *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void cleanup_graph(dBGraph* db_graph)
{
    Orientation orientation;
    
    void clean_node(dBNode * node) {
        void check_edge(Nucleotide nucleotide) {
            // Check if edge exists
            if (db_node_edge_exist_any_colour
                (node, nucleotide, orientation)) {
                // Now see if kmer it leads to exists and is flagged for output
                pathStep current_step, next_step, rev_step;
                current_step.node = node;
                current_step.label = nucleotide;
                current_step.orientation = orientation;
                db_graph_get_next_step(&current_step,
                                       &next_step, &rev_step,
                                       db_graph);
                if ((next_step.node == NULL) || 
                    (!db_node_check_flag_visited(next_step.node))) {
                    db_node_reset_edge_all_colours(node,
                                                   orientation,
                                                   nucleotide);
                }
            }
        }
        
        orientation = forward;
        nucleotide_iterator(&check_edge);
        
        orientation = reverse;
        nucleotide_iterator(&check_edge);
    }
    
    log_and_screen_printf("Cleaning graph\n");
    hash_table_traverse(&clean_node, db_graph);
}

/*----------------------------------------------------------------------*
 * Function: write_ctx_files                                            *
 * Purpose:  Write output CTX files of visited nodes                    *
 * Params:   filename -> base filename to use (no extension)            *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void write_ctx_files(char* f, dBGraph* db_graph)
{
    int c;
    char filename[1024];
    
    for (c=0; c<NUMBER_OF_COLOURS; c++) {
        sprintf(filename, "%s_c%d.ctx", f, c);
        log_and_screen_printf("Writing colour %d as file %s...\n", c, filename);
        db_graph_dump_binary_by_colour(filename, &db_node_check_flag_visited, c, db_graph);
    }
}

/*----------------------------------------------------------------------*
 * Function: parse_string                                               *
 * Purpose:  Return a string parameter (command line parsing)           *
 * Params:   argc = number of arguments                                 *
 *           argv -> argument array                                     *
 *           i -> pointer to argument number counter                    *
 * Returns:  a string and updates i                                     *
 *----------------------------------------------------------------------*/
char* parse_string(int argc, char* argv[], int* i)
{
    char* token = 0;
    
    if (strlen(argv[*i]) > 2)
        token = argv[*i] + 2;
    else if (*i < (argc-1))	{
        *i = *i + 1;
        token = argv[*i];
    }
    
    if ((token) && (token[0] == '-'))
        token = 0;
    
    return token;
}

/*----------------------------------------------------------------------*
 * Function: parse_int                                                  *
 * Purpose:  Return an integer parameter (command line parsing)         *
 * Params:   argc = number of arguments                                 *
 *           argv -> argument array                                     *
 *           i -> pointer to argument number counter                    *
 * Returns:  a number and updates i                                     *
 *----------------------------------------------------------------------*/
int parse_int(int argc, char* argv[], int* i)
{
    int v = -1;
    char* token = parse_string(argc, argv, i);
    
    if (token)
        v = atoi(token);
    
    return v;
}

/*----------------------------------------------------------------------*
 * Function: parse_command_line_args                                    *
 * Purpose:  Deal with command line arguments                           *
 * Params:   argc = number of arguments                                 *
 *           argv -> argument array                                     *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void parse_command_line_args(int argc, char* argv[])
{
    int i = 1;
    char* ks;
    
    if (argc < 5)
    {
        printf("Syntax: graphout [-i filename] [-k int] [options]\n");
        printf("where [-i filename] specifies the name of a file of files (CTX only)\n");
        printf("      [-k int] specifies kmer size\n");
        printf("      [-A filename] will analyse graph structure\n");
        printf("      [-f kmer] specifies a kmer to start from (1 or more may be specified)\n");
        printf("      [-l filename] speficies the name of a debugging log file\n");
        printf("      \nAnalysis optiokns:\n");
        printf("      [-S filename] outputs consenus contigs for each subgraph\n");
        printf("      [-T int] sets the minimum number of nodes in a subgraph to output contigs (default 0)\n");
        printf("      \nGraph output formats:\n");
        printf("      [-O basename] specifies the basename (prefix) of output files\n");
        printf("      [-B] will output a BioLayout file (.txt) for each search kmer\n");
        printf("      [-C] will output CTX files for each search kmer\n");
        printf("      [-D] will output a DOT format file (.gv) for each search kmer\n");
        printf("      [-G] will output a GraphML file (.graphml) for each search kmer\n");
        printf("      [-M] will output a GML file (.gml) for each search kmer\n");
        printf("      [-U] will output a Ubigraph file (.csv) for each search kmer\n");
        printf("      [-X] will output a GEXF file (.gexf.xml) for each search kmer\n");
        printf("      \nGraph output options:\n");        
        printf("      [-a int] specifies the maximum length of path to add each step (default 10)\n");
        printf("      [-d int] specifies the maximum branch level (default 10)\n");
        printf("      [-m int] specifies the maximum number of nodes to output (default 200)\n"); 
        printf("      [-n int] specifies number of buckets in hash table (default 17)\n");
        printf("      [-b int] specifies size of hash table buckets (default 60)\n");
        printf("      [-j] outputs only branching nodes\n");
        printf("      [-s] makes super nodes small (DOT format only)\n");
        printf("      [-c] makes branching nodes use circles instead of ellipses (DOT format only)\n");
        printf("      [-z] cleans output .ctx files to remove arcs to non-existent kmers\n");
        printf("\n");
        exit(1);
    }
    
    while(i < argc)
    {
        char* parameter = argv[i];
        if (parameter[0] == '-')
        {
            char option = parameter[1];
            
            switch (option)
            {
                case 'i':
                    file_of_filenames = parse_string(argc, argv, &i);
                    if (!file_of_filenames) {
                        printf("ERROR: Invalid filename for -i parameter.\n");
                        exit(1);
                    }
                    break;
                case 'A':
                    analysis_filename = parse_string(argc, argv, &i);
                    if (!analysis_filename) {
                        printf("ERROR: Invalid filename for -A parameter.\n");
                        exit(1);
                    }
                    analyse_graph = true;
                    break;
                case 'O':
                    base_filename = parse_string(argc, argv, &i);
                    if (!base_filename) {
                        printf("ERROR: Invalid filename for -O parameter.\n");
                        exit(1);
                    }
                    break;
                case 'S':
                    consensus_contigs_filename = parse_string(argc, argv, &i);
                    if (!consensus_contigs_filename) {
                        printf("ERROR: Invalid filename for -S parameter.\n");
                        exit(1);
                    }
                    break;
                case 'f':
                    ks = parse_string(argc, argv, &i);
                    if (!ks) {
                        printf("ERROR: Invalid filename for -f parameter.\n");
                        exit(1);
                    }										
                    if (number_of_search_kmers >= MAX_SEARCH_KMERS) {
                        printf("Error: Too many search kmers.\n");
                        exit(1);
                    }					
                    search_kmers[number_of_search_kmers++] = ks;
                    break;
                case 'l':
                    debug_filename = parse_string(argc, argv, &i);
                    if (!debug_filename) {
                        printf("ERROR: Invalid filename for -l parameter.\n");
                        exit(1);
                    }
                    break;
                case 'k':
                    kmer_size = parse_int(argc, argv, &i);
                    gt_options.kmer_size = kmer_size;
                    break;
                case 'B': output_biolayout = true; break;
                case 'D': output_graphviz = true; break;
                case 'G': output_graphml = true; break;
                case 'M': output_gml = true; break;
                case 'U': output_ubigraph = true; break;
                case 'X': output_gexf = true; break;
                case 'C': output_ctx = true; break;
                case 'a': gt_options.max_add_length = parse_int(argc, argv, &i); break;
                case 'd': gt_options.max_node_depth = parse_int(argc, argv, &i); break;
                case 'm': gt_options.max_nodes_to_output = parse_int(argc, argv, &i); break;
                case 'n': hash_key_bits = parse_int(argc, argv, &i); break;
                case 'b': bucket_size = parse_int(argc, argv, &i); break;
                case 's': gt_options.make_minor_nodes_small = 1; break;
                case 'c': gt_options.circles_for_major_nodes = 1; break;
                case 'j': gt_options.only_major_nodes = 1;	break;
                case 'z': clean_graph = 1; break;
                case 'T': min_subgraph_kmers = parse_int(argc, argv, &i); break;
                default:
                    printf("ERROR: Invalid option -%c\n", option);
                    exit(1);
            }
        }
        
        i++;
    }
    
    if ((kmer_size < 1) || (kmer_size > 255)) {
        printf("ERROR: Invalid kmer size.\n");
        exit(1);
    }
    
    if (!file_of_filenames) {
        printf("ERROR: You must specify an input file of files.\n");
        exit(1);
    }
    
    if ((!analyse_graph) && (!base_filename)) {
        printf("ERROR: You must specify an output base filename.\n");
        exit(1);
    }
    
    if ((!analyse_graph) && (number_of_search_kmers == 0)) {
        printf("ERROR: You must specify a search kmer.\n");
        exit(1);
    }
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
void load_ctx_files(char* file_of_filenames, dBGraph* graph)
{
    FILE* fp_fnames;
    char filename[1024];
    
    // Load CTX file
    log_and_screen_printf("Loading CTX files\n");
    log_and_screen_printf("Loading file of files %s\n", file_of_filenames);
    
    // Open file of file names
    fp_fnames= fopen(file_of_filenames, "r");
    if (!fp_fnames) {
        log_and_screen_printf("Can't open file of files\n");
        exit(1);
    }
    
    // For each file 
    while (!feof(fp_fnames)) {
        short colour = 0;
        fscanf(fp_fnames, "%s %hd\n", filename, &colour);	
        log_and_screen_printf("Loading %s\n", filename);
        load_binary_from_filename_into_graph(filename, graph, colour, 0);
        log_and_screen_printf("Loaded. Unique kmers: %i\n", graph->unique_kmers);
    }
        
    fclose(fp_fnames);
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
void find_kmers_and_output_graphs(dBGraph* graph)
{
    BinaryKmer b;
    char filename[1024];
    int i;
    int r;
    
    for (i=0; i<number_of_search_kmers; i++) {
        // Find kmer
        char* kmer_string = search_kmers[i];
        log_and_screen_printf("Finding kmer %s\n", kmer_string);
        seq_to_binary_kmer(kmer_string, kmer_size, &b);
        
        // Walk graph
        log_and_screen_printf("Walking\n");
        r = graph_tools_walk_subgraph_for_kmer(&b, &gt_options, &gt_state, graph);
        if (r > 0) {            
            if (output_graphviz) {
                sprintf(filename, "%s_%s.gv", base_filename, kmer_string);
                log_and_screen_printf("Writing file %s...\n", filename);
                graph_tools_output_subgraph_for_kmer(filename, GRAPH_FORMAT_GRAPHVIZ, &gt_options, &gt_state, graph);
            }
            
            if (output_biolayout) {
                sprintf(filename, "%s_%s.txt", base_filename, kmer_string);
                log_and_screen_printf("Writing file %s...\n", filename);
                graph_tools_output_subgraph_for_kmer(filename, GRAPH_FORMAT_BIOLAYOUT, &gt_options, &gt_state, graph);
            }
            
            if (output_ubigraph) {
                sprintf(filename, "%s_%s.csv", base_filename, kmer_string);
                log_and_screen_printf("Writing file %s...\n", filename);
                graph_tools_output_subgraph_for_kmer(filename, GRAPH_FORMAT_UBIGRAPH, &gt_options, &gt_state, graph);
            }
            
            if (output_graphml) {
                sprintf(filename, "%s_%s.graphml", base_filename, kmer_string);
                log_and_screen_printf("Writing file %s...\n", filename);
                graph_tools_output_subgraph_for_kmer(filename, GRAPH_FORMAT_GRAPHML, &gt_options, &gt_state, graph);
            }
            
            if (output_gexf) {
                sprintf(filename, "%s_%s.gexf", base_filename, kmer_string);
                log_and_screen_printf("Writing file %s...\n", filename);
                graph_tools_output_subgraph_for_kmer(filename, GRAPH_FORMAT_GEXF, &gt_options, &gt_state, graph);
            }
            
            if (output_gml) {
                sprintf(filename, "%s_%s.gml", base_filename, kmer_string);
                log_and_screen_printf("Writing file %s...\n", filename);
                graph_tools_output_subgraph_for_kmer(filename, GRAPH_FORMAT_GML, &gt_options, &gt_state, graph);
            }
            
            // Output CTX files
            if (output_ctx) {
                if (clean_graph) {
                    cleanup_graph(graph);
                }
                
                log_and_screen_printf("Writing CTX output\n");
                sprintf(filename, "%s_%s", base_filename, kmer_string);
                write_ctx_files(filename, graph);
            }
        }
    }    
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
int grow_graph_from_node(dBNode* start_node, dBNode** best_node, dBGraph* graph)
{                         
    Queue* nodes_to_walk;
    dBNode* node;
    int orientation;
    int depth;
    int current_graph_size = 0;
    int best_coverage = 0;
    int best_edges = 0;
    
    *best_node = 0;
    
    // Nucleotide iterator, used to walk all possible paths from a node
    void walk_if_exists(Nucleotide n) {
        //if (debug) printf("Trying nucleotide %i\n", n);
        
        // If there is an edge in any colour for this nucleotide...
        if (db_node_edge_exist_any_colour(node, n, orientation)) {
            
            //if (debug) printf("  Edge exists\n");
            
            // Get first node along this edge and check we've not already visited it...
            Orientation next_orientation;
            Nucleotide reverse_nucleotide;
            dBNode * next_node;
            next_node = db_graph_get_next_node(node, orientation, &next_orientation, n, &reverse_nucleotide, graph);
            if (!next_node) {
                log_and_screen_printf("Error: Something went wrong with db_graph_get_next_node\n");
                exit(1);
            }
            
            // If not already visited the first node, walk it...
            if (!db_node_check_flag_visited(next_node)) {
                pathStep first_step;
                Path * new_path;
                dBNode* end_node; 
                int i = 0;
                
                //if (debug) printf("  Not already visited\n");
                
                // Get path				
                first_step.node = node;
                first_step.orientation = orientation;
                first_step.label = n;
                new_path = path_new(gt_options.max_nodes_to_output, graph->kmer_size);
                if (!new_path) {
                    log_and_screen_printf("ERROR: Not enough memory to allocate new path.\n");
                    exit(1);
                }
                
                db_graph_get_perfect_path_with_first_edge_all_colours(&first_step, &db_node_action_do_nothing, new_path, graph);
                
                
                // Add end node to list of nodes to visit
                end_node = new_path->nodes[new_path->length-1];
                if (!db_node_check_flag_visited(end_node)) {
                    if (!db_node_is_blunt_end_all_colours(end_node, new_path->orientations[new_path->length-1])) {
                        if (queue_push_node(nodes_to_walk, end_node, depth+1) == NULL) {
                            log_and_screen_printf("Queue too large. Ending.\n");
                            exit(1);
                        }                        
                    }
                }
                
                // Now go through all nodes, look for best and mark all as visited
                for (i=0; i<new_path->length; i++) {
                    if (!db_node_check_flag_visited(new_path->nodes[i])) {
                        int this_coverage = element_get_coverage_all_colours(new_path->nodes[i]);
                        int this_edges = db_node_edges_count_all_colours(new_path->nodes[i], forward) + db_node_edges_count_all_colours(new_path->nodes[i], reverse);
                        
                        if ((best_node == 0) ||
                            (this_coverage > best_coverage) ||
                            ((this_coverage == best_coverage) && (this_edges < best_edges)))
                        {
                            best_coverage = this_coverage;
                            best_edges = this_edges;
                            *best_node = new_path->nodes[i];                            
                        }
                        
                        db_node_action_set_flag_visited(new_path->nodes[i]);
                        current_graph_size++;                        
                    }
                }
                
                // Clean up
                path_destroy(new_path);
            }
        }
    }
    
    // Start a queue of nodes to walk
    nodes_to_walk = queue_new(ANALYSIS_QUEUE_SIZE);
    if (!nodes_to_walk) {
        log_and_screen_printf("Couldn't get memory for node queue.\n");
        exit(1);
    }
    
    // Add start node to list of nodes to visit
    if (queue_push_node(nodes_to_walk, start_node, 0) == NULL) {
        log_and_screen_printf("Queue too large. Ending.\n");
        exit(1);        
    }
    
    if (!db_node_check_flag_visited(start_node)) {
        db_node_action_set_flag_visited(start_node);
        current_graph_size++;
    }
    
    // Now keep visiting nodes and walking paths
    while (nodes_to_walk->number_of_items > 0) {
        // Take top node from list
        node = queue_pop_node(nodes_to_walk, &depth);
        
        // Look at all paths out from here
        orientation = forward;
        nucleotide_iterator(&walk_if_exists);
        orientation = reverse;
        nucleotide_iterator(&walk_if_exists);				
    }
    
    queue_free(nodes_to_walk);
    
    return current_graph_size;
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
void do_graph_analysis(dBGraph* graph)
{
    FILE* fp;
    int total_nodes = 0;
    SubGraphInfo* sub_graphs;
    int max_length = 200000;
    Path *path_fwd = path_new(max_length, graph->kmer_size);
    Path *path_rev = path_new(max_length, graph->kmer_size);
    Path *final_path = path_new(max_length, graph->kmer_size);
    int n_seeds = 0;
    int i;
    char seq[256];
    
    path_array_initialise_buffers(graph->kmer_size);
      
    sub_graphs = calloc(MAX_SEEDS, sizeof(SubGraphInfo));
    if (!sub_graphs) {
        log_and_screen_printf("ERROR: Can't get memory for subgraphs\n");
        exit(1);
    }
    
    fp = fopen(analysis_filename, "w");
    if (!fp) {
        log_and_screen_printf("ERROR: Can't open analysis file.\n");
        exit(1);
    }
    
    log_and_screen_printf("Allocated space for %i seeds\n", MAX_SEEDS);
    
    // For each node, if it's not pruned or visited, try and grow a graph
    void explore_node(dBNode * node) {
        if (db_node_check_for_any_flag(node, PRUNED | VISITED) == false) {
            int nodes_in_graph;
            
            // Grow graph from this node
            nodes_in_graph = grow_graph_from_node(node, &(sub_graphs[n_seeds].seed_node), graph);
            total_nodes += nodes_in_graph;
            
            binary_kmer_to_seq(&(node->kmer), graph->kmer_size, seq);            
            fprintf(fp, "%i\t%i\t%i\t%s\t", n_seeds, nodes_in_graph, total_nodes, seq);
            sub_graphs[n_seeds].graph_size = nodes_in_graph;
            binary_kmer_to_seq(&(sub_graphs[n_seeds].seed_node->kmer), graph->kmer_size, seq);
            fprintf(fp, "%s\n", seq);
            n_seeds++;
            
            if (n_seeds == MAX_SEEDS) {
                log_and_screen_printf("Error: MAX_SEEDS exceeded. Quitting.\n");
                exit(1);
            }
        }
    }
    
    // Traverse each node...
    log_and_screen_printf("Doing analysis...\n");
    hash_table_traverse(&explore_node, graph);
    log_and_screen_printf("Finished. Total: %i\n", total_nodes);
    fclose(fp);
    
    
    if (consensus_contigs_filename) {
        fp = fopen(consensus_contigs_filename, "w");
        if (!fp) {
            log_and_screen_printf("ERROR: Can't open contig file.\n");
            exit(1);
        }

        db_graph_reset_flags(graph);    
        log_and_screen_printf("Outputting contigs...\n");
        for (i=0; i<n_seeds; i++) {
            log_printf("Graph %i\n", i);           
            if (sub_graphs[i].graph_size >= min_subgraph_kmers) {            
                binary_kmer_to_seq(&(sub_graphs[i].seed_node->kmer), graph->kmer_size, seq);
                //log_printf("  Seed %s Orientation fwd\n", seq);
                coverage_walk_get_path(sub_graphs[i].seed_node, forward, NULL, graph, path_fwd);
                //log_printf("  Seed %s Orientation rev\n", seq);
                coverage_walk_get_path(sub_graphs[i].seed_node, reverse, NULL, graph, path_rev);
                path_reverse(path_fwd, final_path);
                path_append(final_path, path_rev);
                final_path->id = i;
                path_to_fasta(final_path, fp);
                log_printf("  Seed %s\tFwd path length %i\tRev path length %i\tFinal path length %i\n", seq, path_fwd->length, path_rev->length, final_path->length);
                path_reset(path_fwd);
                perfect_path_get_path(sub_graphs[i].seed_node, forward, &db_node_action_do_nothing, graph, path_fwd);
                log_printf("\t\tPerfect path fwd length %i\n", path_fwd->length);
                path_reset(path_rev);
                path_reset(final_path);
            } else {
                log_printf("  Number of nodes (%i} too small. Not outputting contig.\n", sub_graphs[i].graph_size);
            }
            
        }
        log_and_screen_printf("Finished contig output.\n");    
        fclose(fp);
    }

    free(sub_graphs);
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
void display_arguments(void)
{
    int i;
    int total_mem = 0;
    
    log_printf("\ngraphout - generate visualisations of parts of cortex assemblies\n\n");    
    log_and_screen_printf(SVN_VERSION);
    log_and_screen_printf(SVN_COMMIT_DATE);
    log_and_screen_printf("Compiled on %s at %s \n\n", __DATE__, __TIME__);
    
    log_and_screen_printf("Memory requirements:\n");
    log_and_screen_printf("Seed information: %i Mb\n", (MAX_SEEDS * sizeof(dBNode*)) / 1024 / 1024);
    log_and_screen_printf("      Graph size: %i Mb\n", (MAX_SEEDS * sizeof(dBNode*)) / 1024 / 1024);
    total_mem += MAX_SEEDS * sizeof(dBNode*) * 2;
    log_and_screen_printf("      Node queue: %i Mb\n", (ANALYSIS_QUEUE_SIZE * sizeof(QueueItem*)) / 1024 / 1024);
    total_mem += ANALYSIS_QUEUE_SIZE * sizeof(QueueItem*);
    log_and_screen_printf("       All nodes: %i Mb\n", (ANALYSIS_QUEUE_SIZE * sizeof(QueueItem))/ 1024 / 1024);
    total_mem += ANALYSIS_QUEUE_SIZE * sizeof(QueueItem);
    log_and_screen_printf("           Total: %i Mb\n\n", total_mem / 1024 / 1024);
    
    // Print arguments to screen and log
    log_and_screen_printf("Input file of files: %s\n", file_of_filenames);
    
    if (consensus_contigs_filename) {
        log_and_screen_printf(" Consensus filename: %s\n", consensus_contigs_filename);        
        log_and_screen_printf(" Min subgraph kmers: %i\n", min_subgraph_kmers);
    }
    
    if (base_filename) {
        log_and_screen_printf("   Output base name: %s\n", base_filename);
    }
    
    for (i=0; i<number_of_search_kmers; i++) {
        log_and_screen_printf("        Search kmer: %s\n", search_kmers[i]);
    }
    log_and_screen_printf("          Kmer size: %d\n", kmer_size);
    log_and_screen_printf("     Max node depth: %d\n", gt_options.max_node_depth);
    log_and_screen_printf(" Max nodes in graph: %d\n", gt_options.max_nodes_to_output);
    log_and_screen_printf("      Hash key bits: %d\n", hash_key_bits);
    log_and_screen_printf("   Hash key buckets: %d\n", bucket_size);
    log_and_screen_printf("  Minor nodes small: %d\n\n", gt_options.make_minor_nodes_small);    
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
void handler(int sig)
{
    void *array[10];
    size_t size;
    
    // get void*'s for all entries on the stack
    size = backtrace(array, 10);
    
    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, 2);
    exit(1);
}

/*----------------------------------------------------------------------*
 * Function: main                                                       *
 * Purpose:  Entry point to program.                                    *
 * Params:   argc = number of arguments                                 *
 *           argv -> array of arguments                                 *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
int main (int argc, char * argv[])
{
    HashTable* db_graph = NULL;
    
    graph_tools_initialise_options(&gt_options);
    
    printf("\ngraphout - generate visualisations of parts of cortex assemblies\n\n");
    
    // Install signal handler for catching segmentation faults
    signal(SIGSEGV, handler);
    
    // Read command line arguments
    parse_command_line_args(argc, argv);
    
    // Open debug log
    log_start(debug_filename);
    
    // Dispplay program arguments
    display_arguments();
    
    // Create hash table
    log_and_screen_printf("Creating hash table (element size %i)\n", sizeof(Element));
    db_graph = hash_table_new(hash_key_bits, bucket_size, 20, kmer_size);
    if (!db_graph) {
        log_and_screen_printf("Error: Couldn't get memory for graph\n");
        exit(1);
    }
    
    // Load CTX files
    load_ctx_files(file_of_filenames, db_graph);
    
#ifdef DEBUG_CLEANUP
    db_graph_cleanup_graph(db_graph);	
#endif	
    
    // We either analyse the graph or find kmers and output graphs
    if (analyse_graph) {
        do_graph_analysis(db_graph);
    } else {
        find_kmers_and_output_graphs(db_graph);
    }
    
    log_and_screen_printf("Finished.\n");
    
    return 0;
}
