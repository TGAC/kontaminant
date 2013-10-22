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

typedef struct {
    int kmer_size;
    boolean only_major_nodes;
    boolean circles_for_major_nodes;
    boolean make_minor_nodes_small;
    boolean is_major_node;
    int max_nodes_to_output;
    int max_add_length;
    int max_node_depth;
    dBNode* starting_node;
} GraphToolsOptions;

typedef struct {
    int number_of_marked_nodes;
    int max_nodes_to_output;
    Element** marked_nodes;
    int number_of_node_ids;
    int max_node_ids;
    Element** node_ids;
} GraphToolsState;

void graph_tools_initialise_options(GraphToolsOptions* options);
void graph_tools_add_and_flag_node(dBNode * node, GraphToolsState* state);
int graph_tools_get_marked_node_index(Element* e, GraphToolsState* state);
int graph_tools_walk_around(GraphToolsOptions* options, GraphToolsState* state, dBGraph* graph);
int graph_tools_walk_subgraph_for_kmer(BinaryKmer* start_kmer, GraphToolsOptions* options, GraphToolsState* state, dBGraph* db_graph);
void graph_tools_output_subgraph_for_kmer(char* filename, int graph_format, GraphToolsOptions* options, GraphToolsState* state, dBGraph* db_graph);
