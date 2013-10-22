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
#include <string.h>
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
#include "metacortex.h"

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
                exit(-1);
            }
            
            // If not already visited the first node, walk it...
            if (!db_node_check_flag_visited(next_node)) {
                pathStep first_step;
                Path * new_path;
                dBNode* end_node; 
                int i = 0;
                                
                // Get path				
                first_step.node = node;
                first_step.orientation = orientation;
                first_step.label = n;
                new_path = path_new(MAX_EXPLORE_NODES, graph->kmer_size);
                if (!new_path) {
                    log_and_screen_printf("ERROR: Not enough memory to allocate new path.\n");
                    exit(-1);
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
    //log_and_screen_printf("Allocating %d Mb to store queue information (max %d nodes, when full each node could be %d)...\n", ((METACORTEX_QUEUE_SIZE * sizeof(QueueItem*)) / 1024) / 1024, METACORTEX_QUEUE_SIZE, sizeof(QueueItem));
    nodes_to_walk = queue_new(METACORTEX_QUEUE_SIZE);
    if (!nodes_to_walk) {
        log_and_screen_printf("Couldn't get memory for node queue.\n");
        exit(-1);
    }
    
    // Add start node to list of nodes to visit
    if (queue_push_node(nodes_to_walk, start_node, 0) == NULL) {
        log_and_screen_printf("Queue too large. Ending.\n");
        exit(-1);        
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
    
    // If we didn't find a start node, presumably this is a singleton?
    if (*best_node == 0) {
        printf("Note: didn't find a best node, setting to start node\n");
        *best_node = start_node;
    }
    
    return current_graph_size;
}

void metacortex_find_subgraphs(dBGraph* graph, char* consensus_contigs_filename, int min_subgraph_kmers)
{
    SubGraphInfo* sub_graphs;
    FILE* fp;
    Path *path_fwd = path_new(MAX_EXPLORE_PATH_LENGTH, graph->kmer_size);
    Path *path_rev = path_new(MAX_EXPLORE_PATH_LENGTH, graph->kmer_size);
    Path *final_path = path_new(MAX_EXPLORE_PATH_LENGTH, graph->kmer_size);
    char seq[256];
    char analysis_filename[strlen(consensus_contigs_filename) + 10];
    long int total_nodes = 0;
    int n_seeds = 0;
    int i;
    
    sprintf(analysis_filename, "%s.analysis", consensus_contigs_filename);
    log_and_screen_printf("Running metacortex subgraph analysis...\n");
    log_and_screen_printf("          Contig file: %s\n", consensus_contigs_filename);
    log_and_screen_printf("        Analysis file: %s\n", analysis_filename);
    log_and_screen_printf("Minimum subgraph size: %i\n", min_subgraph_kmers);
    
    /* Initialise temporaray path array buffers */
    path_array_initialise_buffers(graph->kmer_size);
    
    /* Create a list of subgraphs */
    log_and_screen_printf("Allocating %d Mb to store subgraph information (max %d seeds)...\n", ((MAX_SEEDS * sizeof(SubGraphInfo)) / 1024) / 1024, MAX_SEEDS);
    sub_graphs = calloc(MAX_SEEDS, sizeof(SubGraphInfo));
    if (!sub_graphs) {
        log_and_screen_printf("ERROR: Can't get memory for subgraphs\n");
        exit(-1);
    }

    /* Open the analysis file */
    fp = fopen(analysis_filename, "w");
    if (!fp) {
        log_and_screen_printf("ERROR: Can't open analysis file.\n");
        exit(-1);
    }
        
    /* For each node, if it's not pruned or visited, try and grow a graph */
    void explore_node(dBNode * node) {
        if (node == NULL) {
            log_and_screen_printf("Error: NULL node passed to explore_node.\n");
            exit(-1);
        }
        
        if (db_node_check_for_any_flag(node, PRUNED | VISITED) == false) {
            int nodes_in_graph;
            
            /* Grow graph from this node, returning the 'best' (highest coverage) node to store as seed point */
            nodes_in_graph = grow_graph_from_node(node, &(sub_graphs[n_seeds].seed_node), graph);
            total_nodes += nodes_in_graph;
            
            if (sub_graphs[n_seeds].seed_node == NULL) {
                printf("ERROR: Seed node is NULL, nodes in graph is %d\n", nodes_in_graph);
            } else {
                /* Write data to analysis file */
                binary_kmer_to_seq(&(node->kmer), graph->kmer_size, seq);            
                fprintf(fp, "%i\t%i\t%ld\t%s\t", n_seeds, nodes_in_graph, total_nodes, seq);
                binary_kmer_to_seq(&(sub_graphs[n_seeds].seed_node->kmer), graph->kmer_size, seq);
                fprintf(fp, "%s\n", seq);

                /* Store nodes in this subgraph */
                sub_graphs[n_seeds].graph_size = nodes_in_graph;
                n_seeds++;
                
                /* Check we've not run out of seed storage - in future, this should dynamically allocate */
                if (n_seeds == MAX_SEEDS) {
                    log_and_screen_printf("Error: MAX_SEEDS exceeded. Quitting.\n");
                    exit(-1);
                }
            }
        }
    }
    
    /* Traverse each node... */
    log_and_screen_printf("Finding subgraphs...\n");
    hash_table_traverse(&explore_node, graph);
    log_and_screen_printf("Finished. Total: %ld\n", total_nodes);
    fclose(fp);    
    
    /* Open consensus contigs file */
    fp = fopen(consensus_contigs_filename, "w");
    if (!fp) {
        log_and_screen_printf("ERROR: Can't open contig file.\n");
        exit(-1);
    }
    
    /* Now go through all the seed points and generate the consensus contigs by walking forward and backward from the seed */
    db_graph_reset_flags(graph);    
    log_and_screen_printf("Outputting contigs...\n");
	log_progress_bar(0);
	long long one_percent = n_seeds/100;
    int percent;
    
    if (one_percent < 1) {
        one_percent = 1;
    }
    
    for (i=0; i<n_seeds; i++) {
        if (i % one_percent == 0) {
            percent = (100 * i) / n_seeds;
            log_progress_bar(percent);
        } 
        
        //log_printf("Graph %i\n", i);           
        if (sub_graphs[i].graph_size >= min_subgraph_kmers) {            
            binary_kmer_to_seq(&(sub_graphs[i].seed_node->kmer), graph->kmer_size, seq);
            coverage_walk_get_path(sub_graphs[i].seed_node, forward, NULL, graph, path_fwd);
            coverage_walk_get_path(sub_graphs[i].seed_node, reverse, NULL, graph, path_rev);
            path_reverse(path_fwd, final_path);
            path_append(final_path, path_rev);
            final_path->id = i;
            path_to_fasta(final_path, fp);
            //log_printf("  Seed %s\tFwd path length %i\tRev path length %i\tFinal path length %i\n", seq, path_fwd->length, path_rev->length, final_path->length);
            path_reset(path_fwd);
            perfect_path_get_path(sub_graphs[i].seed_node, forward, &db_node_action_do_nothing, graph, path_fwd);
            //log_printf("\t\tPerfect path fwd length %i\n", path_fwd->length);
            path_reset(path_rev);
            path_reset(final_path);
        } else {
            log_printf("  Number of nodes (%i} too small. Not outputting contig.\n", sub_graphs[i].graph_size);
        }
        
    }
	log_progress_bar(100);
	printf("\n");
    log_and_screen_printf("Finished contig output.\n");    
    fclose(fp);
    
    free(sub_graphs);
}
