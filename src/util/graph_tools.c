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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <execinfo.h>
#include <signal.h>
#include <ctype.h>
#include <limits.h>
#include <stdint.h>
#include <assert.h>
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

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
void graph_tools_initialise_options(GraphToolsOptions* options)
{
    options->kmer_size = 0;
    options->only_major_nodes = false;
    options->circles_for_major_nodes = false;
    options->make_minor_nodes_small = false;
    options->is_major_node = false;
    options->max_nodes_to_output = 200;
    options->max_add_length = 10;
    options->max_node_depth = 10;
    options->starting_node = 0;
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
void graph_tools_free_state(GraphToolsState* state)
{
    if (state->marked_nodes) {
        free(state->marked_nodes);
    }
    state->number_of_marked_nodes = 0;

    if (state->node_ids) {
        free(state->node_ids);
    }
    state->max_node_ids = 0;
    state->number_of_node_ids = 0;

}

/*----------------------------------------------------------------------*
 * Function: add_and_flag_node                                          *
 * Purpose:  Add a node to the list of marked nodes and flag it has     *
 *           visited.                                                   *
 * Params:   node -> node to add                                        *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void graph_tools_add_and_flag_node(dBNode * node, GraphToolsState* state)
{	
    if (!db_node_check_flag_visited(node)) {
        if (state->number_of_marked_nodes < state->max_nodes_to_output) {
            state->marked_nodes[state->number_of_marked_nodes++] = node;
            db_node_action_set_flag_visited(node);
        }
    }
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
int graph_tools_get_marked_node_index(Element* e, GraphToolsState* state)
{
    int index = -1;
    int i;
    
    for (i=0; i<state->number_of_marked_nodes; i++) {
        if (state->marked_nodes[i] == e) {
            index = i;
            break;
        }
    }
    
    return index;
}

/*----------------------------------------------------------------------*
 * Function: walk_around                                                *
 * Purpose:  Given a starting node, walk around in all directions,      *
 *           marking nodes.                                             *
 * Params:   start_node -> starting node in graph                       *
 *           graph -> the hash table structure                          *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
int graph_tools_walk_around(GraphToolsOptions* options, GraphToolsState* state, dBGraph* graph)
{
    Queue* nodes_to_walk;
    dBNode* node;
    int nodes_examined = 0;
    int orientation;
    int depth;
    
    // Nucleotide iterator, used to walk all possible paths from a node
    void walk_if_exists(Nucleotide n) {
        // If there is an edge in any colour for this nucleotide...
        if (db_node_edge_exist_any_colour(node, n, orientation)) {
            // Get first node along this edge and check we've not already visited it...
            Orientation next_orientation;
            Nucleotide reverse_nucleotide;
            dBNode * next_node;
            next_node = db_graph_get_next_node(node, orientation, &next_orientation, n, &reverse_nucleotide, graph);
            if (!next_node) {
                printf("Error: Something went wrong with db_graph_get_next_node\n");
                exit(1);
            }
            
            // If not already visited the first node, walk it...
            if (!db_node_check_flag_visited(next_node)) {
                int i, loop_end;
                pathStep first_step;
                Path * new_path;
                dBNode* end_node; 
                
                // Get path				
                first_step.node = node;
                first_step.orientation = orientation;
                first_step.label = n;
                new_path = path_new(options->max_nodes_to_output, graph->kmer_size);
                if (!new_path) {
                    printf("ERROR: Not enough memory to allocate new path.\n");
                    exit(1);
                }
                
                db_graph_get_perfect_path_with_first_edge_all_colours(&first_step, &db_node_action_do_nothing, new_path, graph);
                
                // Mark path as visited
                if ((new_path->length-1) > options->max_add_length) {
                    loop_end = options->max_add_length;
                } else {
                    loop_end = new_path->length-1;
                }
                
                for (i=1; i<loop_end; i++) {
                    graph_tools_add_and_flag_node(new_path->nodes[i], state);
                }
                
                // Add end node to list of nodes to visit
                end_node = new_path->nodes[loop_end];
                if (!db_node_check_flag_visited(end_node)) {
                    graph_tools_add_and_flag_node(end_node, state);
                    if ((depth < options->max_node_depth) && (state->number_of_marked_nodes < options->max_nodes_to_output)) {
                        if (queue_push_node(nodes_to_walk, end_node, depth+1) == NULL) {
                            log_and_screen_printf("Queue too large. Ending.\n");
                            exit(1);
                        }
                    }
                }
                
                if (!options->only_major_nodes) {
                    // Free path
                    path_destroy(new_path);
                }
            }
        }
    }
    
    // Start a queue of nodes to walk
    nodes_to_walk = queue_new(options->max_nodes_to_output);
    if (!nodes_to_walk) {
        printf("Couldn't get memory for node queue.\n");
        exit(1);
    }
    
    // Add start node to list of nodes to visit
    graph_tools_add_and_flag_node(options->starting_node, state);
    if (queue_push_node(nodes_to_walk, options->starting_node, 0) == NULL) {
        log_and_screen_printf("Queue too large. Ending.\n");
        exit(1);
    }
    
    // Now keep visiting nodes and walking paths
    while ((nodes_to_walk->number_of_items > 0) && (state->number_of_marked_nodes < options->max_nodes_to_output)) {
        // Take top node from list
        node = queue_pop_node(nodes_to_walk, &depth);
        nodes_examined++;
        
        // Look at all paths out from here
        orientation = forward;
        nucleotide_iterator(&walk_if_exists);
        orientation = reverse;
        nucleotide_iterator(&walk_if_exists);				
    }
    
    queue_free(nodes_to_walk);
    
    return nodes_examined;
}


/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
int graph_tools_walk_subgraph_for_kmer(BinaryKmer* start_kmer, GraphToolsOptions* options, GraphToolsState* state, dBGraph* db_graph)
{
    BinaryKmer tmp_kmer;
    int r = 0;
    
    graph_tools_free_state(state);
    
    state->max_nodes_to_output = options->max_nodes_to_output;
    state->number_of_marked_nodes = 0;
    state->marked_nodes = calloc(state->max_nodes_to_output, sizeof(Element*));
    state->max_node_ids = 100000;
    state->number_of_node_ids = 0;
    state->node_ids = calloc(state->max_node_ids, sizeof(Element*));
    
    if (!state->marked_nodes) {
        log_and_screen_printf("Error: can't get memory for marked nodes array.\n");
        exit(1);
    }
    
    options->starting_node = hash_table_find(element_get_key(start_kmer, options->kmer_size, &tmp_kmer), db_graph);
    if (options->starting_node != NULL) {
        hash_table_traverse(&db_node_action_clear_flags, db_graph);
        r = graph_tools_walk_around(options, state, db_graph);
    } else {
        char seq[1024];
        binary_kmer_to_seq(start_kmer, db_graph->kmer_size, seq);     
        log_and_screen_printf("Error: can't find kmer %s\n", seq);        
    }
    
    return r;
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
void graph_tools_clear_node_ids(GraphToolsState* state)
{
    state->number_of_node_ids = 0;
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
int graph_tools_get_node_id(dBNode * node, GraphToolsState* state)
{
    int i;
    int id = -1;
    
    assert(state->number_of_node_ids < state->max_node_ids);
    
    for (i=0; i<state->number_of_node_ids; i++) {
        if (state->node_ids[i] == node) {
            id = i;
            break;
        }
    }
    
    if (id == -1) {
        id = state->number_of_node_ids;
        state->node_ids[state->number_of_node_ids] = node;
        state->number_of_node_ids = state->number_of_node_ids + 1;
    }
    
    return id;
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
void graph_tools_prepare_node(dBNode* node, GraphToolsOptions* options, GraphToolsState* state, GraphoutNode *gon)
{
    char seq1[options->kmer_size+1];
    BinaryKmer tmp_kmer;
    binary_kmer_assignment_operator(tmp_kmer, node->kmer);
    binary_kmer_to_seq(&tmp_kmer, options->kmer_size, seq1);
    strcpy(gon->seq, seq1);
    gon->node = node;
    gon->id = graph_tools_get_node_id(node, state);
}

/*----------------------------------------------------------------------*
 * Function: write_graph_file                                           *
 * Purpose:  Write a file of all marked nodes.                          *
 * Params:   filename -> output filename                                *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void graph_tools_write_graph_file(char * filename, GraphFileFormat* gff, GraphToolsOptions* options, GraphToolsState* state, dBGraph* db_graph)
{
    FILE* fp;
    int i;
    int edge_counter = 0;
    
    log_printf("In graph_tools_write_graph_file\n");    
    
    // Open file and write header
    fp = fopen(filename, "w");
    if (!fp)
    {
        printf("Can't open file %s\n", filename);
    } else {
        
        if (gff->write_header) {
            gff->write_header(fp, options);
        }
        
        // 'Pre' node writing
        if (gff->write_pre_node) {
            for (i=0; i<state->number_of_marked_nodes; i++) {
                dBNode* node = state->marked_nodes[i];
                
                // If only outputting major nodes, we don't output a node with only one edge in either direction
                if (options->only_major_nodes) {
                    if (db_node_edges_count_all_colours(node, forward) == 1 &&
                        db_node_edges_count_all_colours(node, reverse) == 1) {
                        node = 0;
                    }
                }
                
                // If we've found a node...
                if (node != NULL) {
                    GraphoutNode n;
                    graph_tools_prepare_node(node, options, state, &n);
                    gff->write_pre_node(fp, &n, options);
                }
            }
        }       
        
        log_printf("Number of marked nodes: %i\n", state->number_of_marked_nodes);
        
        // Loop through marked nodes
        for (i=0; i<state->number_of_marked_nodes; i++) {
            dBNode* node = state->marked_nodes[i];
            
            // If only outputting major nodes, we don't output a node with only one edge in either direction
            if (options->only_major_nodes) {
                if (db_node_edges_count_all_colours(node, forward) == 1 &&
                    db_node_edges_count_all_colours(node, reverse) == 1) {
                    node = 0;
                }
            }
            
            // If we've found a node...
            if (node != NULL) {
                BinaryKmer tmp_kmer, tmp_rev;
                short kmer_size = db_graph->kmer_size;
                char seq[kmer_size+1], seqNext[kmer_size+1], seq1[kmer_size+1];
                Orientation o;
                Edges all_edges;
                Edges ea[NUMBER_OF_COLOURS];
                int i;
                
                // Get kmer sequence
                binary_kmer_assignment_operator(tmp_kmer, node->kmer);
                binary_kmer_reverse_complement(&tmp_kmer, kmer_size, &tmp_rev);
                binary_kmer_to_seq(&tmp_kmer, kmer_size, seq1);
                
                // Write node
                if (gff->write_mid_node) {
                    GraphoutNode n;
                    graph_tools_prepare_node(node, options, state, &n);
                    gff->write_mid_node(fp, &n, options);
                }
                
                // Iterator used to write connections between nodes (edges)
                void hasEdge(Nucleotide n) {
                    BinaryKmer bk;
                    Key k = &bk;
                    boolean colour0edge = (((ea[0] >> n) & 1) == 1);
                    boolean colour1edge = (((ea[1] >> n) & 1) == 1);
                    char* label_colour = (colour0edge && colour1edge) ? "black" : (colour0edge ? "green" : "orange");
                    char edge_string[4096];
                    dBNode* node_to;
                    
                    if (colour0edge || colour1edge) {
                        // If outputting only major nodes, then we have to build a perfect path to the next major node
                        if (options->only_major_nodes) {
                            Path* new_path = path_new(options->max_nodes_to_output, db_graph->kmer_size);
                            pathStep first_step;
                            first_step.node = node;
                            first_step.orientation = o;
                            first_step.label = n;
                            
                            db_graph_get_perfect_path_with_first_edge_all_colours(&first_step, &db_node_action_do_nothing, new_path, db_graph);						
                            
                            binary_kmer_to_seq(&(node->kmer), kmer_size, seq);
                            binary_kmer_to_seq(&(new_path->nodes[new_path->length-1]->kmer), kmer_size, seqNext);
                            node_to = new_path->nodes[new_path->length-1];
                            
                            if (strlen(new_path->seq) <= (kmer_size*2)) {
                                strcpy(edge_string, new_path->seq);
                            } else {
                                strncpy(edge_string, new_path->seq, kmer_size);
                                strcpy(edge_string+kmer_size, "...");
                                strcat(edge_string, new_path->seq+strlen(new_path->seq)-kmer_size);
                            }							
                            
                            path_destroy(new_path);
                        }
                        // Otherwise we just print one label
                        else 
                        {
                            binary_kmer_modify_base(&tmp_kmer, n, kmer_size, 0);
                            element_get_key(&tmp_kmer, kmer_size, k);
                            node_to = hash_table_find(k, db_graph);
                            assert(node_to != NULL);
                            binary_kmer_to_seq(k, kmer_size, seqNext);
                            sprintf(edge_string, "%c", binary_nucleotide_to_char(n));
                        }
                        
                        if (gff->write_edge) {
                            if ((gff->can_handle_missing_nodes) ||
                                ((graph_tools_get_marked_node_index(node, state) != -1) && (graph_tools_get_marked_node_index(node_to, state) != -1)))
                            {
                                GraphoutEdge e;
                                strcpy(e.from, seq);
                                strcpy(e.to, seqNext);
                                strcpy(e.label, edge_string);
                                strcpy(e.label_colour, label_colour);
                                e.orientation = o;
                                e.id = edge_counter++;
                                e.id_from = graph_tools_get_node_id(node, state);
                                e.id_to = graph_tools_get_node_id(node_to, state);
                                gff->write_edge(fp, &e, options);
                            }
                        }
                    }
                }
                
                // Write forward edges
                binary_kmer_to_seq(&tmp_kmer, kmer_size, seq);
                binary_kmer_left_shift(&tmp_kmer, 2, kmer_size);
                o = forward;
                all_edges = db_node_get_edges_for_orientation_all_colours(node, o);
                for (i=0; i<NUMBER_OF_COLOURS; i++) {
                    ea[i] = db_node_get_edges_for_orientation_by_colour(node, o, i);
                }
                nucleotide_iterator(&hasEdge);
                
                // Write reverse edges
                binary_kmer_assignment_operator(tmp_kmer, tmp_rev);
                binary_kmer_left_shift(&tmp_kmer, 2, kmer_size);			
                o = reverse;
                all_edges = db_node_get_edges_for_orientation_all_colours(node, o);
                for (i=0; i<NUMBER_OF_COLOURS; i++) {
                    ea[i] = db_node_get_edges_for_orientation_by_colour(node, o, i);			
                }
                nucleotide_iterator(&hasEdge);			
            }
        }
        
        // Footer
        if (gff->write_footer) {
            gff->write_footer(fp, options);
        }
        
        // BioLayout node classes have to come at the end of file
        if (gff->write_post_node) {
            for (i=0; i<state->number_of_marked_nodes; i++) {
                dBNode* node = state->marked_nodes[i];
                
                // If only outputting major nodes, we don't output a node with only one edge in either direction
                if (options->only_major_nodes) {
                    if (db_node_edges_count_all_colours(node, forward) == 1 &&
                        db_node_edges_count_all_colours(node, reverse) == 1) {
                        node = 0;
                    }
                }
                
                // If we've found a node...
                if (node != NULL) {
                    GraphoutNode n;
                    graph_tools_prepare_node(node, options, state, &n);
                    gff->write_post_node(fp, &n, options);
                }
            }
        }
        
        // Close file
        fclose(fp);
    }
}

/*----------------------------------------------------------------------*
 * Function:                                                            *
 * Purpose:                                                             *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
void graph_tools_output_subgraph_for_kmer(char* filename, int graph_format, GraphToolsOptions* options, GraphToolsState* state, dBGraph* db_graph)
{
    GraphFileFormat gff;
    
    log_printf("In graph_tools_output_subgraph_for_kmer\n");
    
    graph_tools_clear_node_ids(state);
    
    switch(graph_format) {
        case GRAPH_FORMAT_GRAPHVIZ: setup_file_format_graphviz(&gff); break;
        case GRAPH_FORMAT_GRAPHML: setup_file_format_graphml(&gff); break;
        case GRAPH_FORMAT_BIOLAYOUT: setup_file_format_biolayout(&gff); break;
        case GRAPH_FORMAT_UBIGRAPH: setup_file_format_ubigraph(&gff); break;
        case GRAPH_FORMAT_GEXF: setup_file_format_gexf(&gff); break;
        case GRAPH_FORMAT_GML: setup_file_format_gml(&gff); break;
    }
    
    graph_tools_write_graph_file(filename, &gff, options, state, db_graph);
}


