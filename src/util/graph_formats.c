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
 * File:    graph_formats.c                                             *
 * Purpose: Graph file format functions for graphout                    *
 * Author:  Richard Leggett                                             *
 * History: 15-Feb-11: RML: Created                                     *
 *----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "binary_kmer.h"
#include "flags.h"
#include "element.h"
#include "dB_graph.h"
#include "graph_tools.h"
#include "graph_formats.h"

/*----------------------------------------------------------------------*
 * GraphViz (.dot)                                                      *
 *----------------------------------------------------------------------*/
void write_graphviz_header(FILE* fp, GraphToolsOptions* options)
{
    fprintf(fp, "digraph finite_state_machine {\n");
    fprintf(fp, "\trankdir=LR;\n");
    fprintf(fp, "\tsize=\"12,8\";\n");
    fprintf(fp, "\tfontsize=18;\n");
    fprintf(fp, "\tnode [shape=rectangle, color=grey];\n");
}

void write_graphviz_footer(FILE* fp, GraphToolsOptions* options)
{
	fprintf(fp, "}\n");
}

void write_graphviz_node(FILE* fp, GraphoutNode *n, GraphToolsOptions* options)
{
    char label_string[512];
    char reverse_seq1[options->kmer_size+1];
    char node_shape[64];
    char node_colour[64];
    int is_major_node = 0;    
        
    // Choose node colour - green for colour 0, orange for colour 1, black for both
    if ((n->node->coverage[0] > 0) && (n->node->coverage[1] == 0)) {
        strcpy(node_colour, "green");
    } else if ((n->node->coverage[0] == 0) && (n->node->coverage[1] > 0)) {
        strcpy(node_colour, "orange");
    } else {
        strcpy(node_colour, "black");
    }
    
    // Chose node shape	
    if (db_node_edges_count_all_colours(n->node, forward) > 1 ||
        db_node_edges_count_all_colours(n->node, reverse) > 1) {
        if (options->circles_for_major_nodes) {
            strcpy(node_shape, "circle");
        } else {
            strcpy(node_shape, "ellipse");
        }
        is_major_node = 1;
    } else {
        if (options->make_minor_nodes_small) {
            strcpy(node_shape, "point");
        } else {
            strcpy(node_shape, "ellipse");
        }
    }    
    
    seq_reverse_complement(n->seq, options->kmer_size, reverse_seq1);
    if (options->make_minor_nodes_small && !options->is_major_node) {
        sprintf(label_string, " ");
    } else {
        if (options->kmer_size < 8) {
            // Bodge - miss off the top line (total coverage), to make ellipses not look like circles
            sprintf(label_string, "<<table border=\"0\" cellpadding=\"0\" cellspacing=\"0\">"
                    "<tr><td><font color=\"blue\">%s</font></td></tr>"
                    "<tr><td><font color=\"red\">%s</font></td></tr>"
                    "<tr><td><font color=\"green\">%d</font> <font color=\"orange\">%d</font></td></tr>"
                    "</table>>",
                    n->seq, reverse_seq1, n->node->coverage[0], n->node->coverage[1]);
        } else {
            // Includes top line that represents total coverage
            sprintf(label_string, "<<table border=\"0\" cellpadding=\"0\" cellspacing=\"0\">"
                    "<tr><td><font color=\"black\">%d</font></td></tr>"
                    "<tr><td><font color=\"blue\">%s</font></td></tr>"
                    "<tr><td><font color=\"red\">%s</font></td></tr>"
                    "<tr><td><font color=\"green\">%d</font> <font color=\"orange\">%d</font></td></tr>"
                    "</table>>",
                    n->node->coverage[0] + n->node->coverage[1], n->seq, reverse_seq1, n->node->coverage[0], n->node->coverage[1]);
        }
    }
    
    if (n->node == options->starting_node) {
        fprintf(fp, "%s [label=%s, style=filled, fillcolor=grey90, shape=%s, color=%s]\n", n->seq, label_string, node_shape, node_colour);
    } else {
        fprintf(fp, "%s [label=%s, shape=%s, color=%s]\n", n->seq, label_string, node_shape, node_colour);
    }
}

void write_graphviz_edge(FILE* fp, GraphoutEdge *e, GraphToolsOptions* options)
{
    fprintf(fp,
            "%s -> %s [ label=<<font color=\"%s\">%s</font>>, color=\"%s\"];\n",
            e->from,
            e->to,
            e->label_colour,
            e->label,
            e->orientation==forward ? "blue" : "red");
}

GraphFileFormat* setup_file_format_graphviz(GraphFileFormat* gff)
{
    gff->write_header = &write_graphviz_header;
    gff->write_pre_node = 0;
    gff->write_mid_node = &write_graphviz_node;
    gff->write_edge = &write_graphviz_edge;
    gff->write_post_node = 0;
    gff->write_footer = &write_graphviz_footer;
    gff->can_handle_missing_nodes = true;
    return gff;
}

/*----------------------------------------------------------------------*
 * BioLayout                                                            *
 *----------------------------------------------------------------------*/
void write_biolayout_edge(FILE* fp, GraphoutEdge *e, GraphToolsOptions* options)
{
    fprintf(fp, "%s\t%s\t%s\n", e->from, e->to, e->label);
}

void write_biolayout_node_class(FILE* fp, GraphoutNode *n, GraphToolsOptions* options)
{
    int is_major_node=0;
    char class[64];
    int size=4;
    
    if (db_node_edges_count_all_colours(n->node, forward) > 1 ||
        db_node_edges_count_all_colours(n->node, reverse) > 1) {
        is_major_node = 1;        
    }
    
    if ((n->node->coverage[0] > 0) && (n->node->coverage[1] > 0)) {
        strcpy(class, "c0_and_c1");
    } else if (n->node->coverage[0] > 0) {
        strcpy(class, "c0");
    } else if (n->node->coverage[1] > 0) {
        strcpy(class, "c1");
    }
    
    if (n->node == options->starting_node) {
        strcat(class, "_start");
    }
    
    if (is_major_node) {
        size = size * 2;
    }
    
    fprintf(fp, "//NODECLASS\t%s\t%s\n", n->seq, class);
    fprintf(fp, "//NODESIZE\t%s\t%d\n", n->seq, size);
}

GraphFileFormat* setup_file_format_biolayout(GraphFileFormat* gff)
{
    gff->write_header = 0;
    gff->write_pre_node = 0;
    gff->write_mid_node = 0;
    gff->write_edge = &write_biolayout_edge;
    gff->write_post_node = &write_biolayout_node_class;
    gff->write_footer = 0;
    gff->can_handle_missing_nodes = true;
    return gff;
}

/*----------------------------------------------------------------------*
 * Ubigraph                                                             *
 *----------------------------------------------------------------------*/
void write_ubigraph_node(FILE* fp, GraphoutNode* n, GraphToolsOptions* options)
{
    int is_major_node;
    
    if (db_node_edges_count_all_colours(n->node, forward) > 1 ||
        db_node_edges_count_all_colours(n->node, reverse) > 1) {
        is_major_node = 1;        
    }
    
    fprintf(fp, "N %s,%d,%d,%d,%d\n",
            n->seq,
            n->node->coverage[0],
            n->node->coverage[1],
            is_major_node,
            n->node == options->starting_node ? 1:0);
}

void write_ubigraph_edge(FILE* fp, GraphoutEdge* e, GraphToolsOptions* options)
{
    fprintf(fp, "E %s,%s,%s\n", e->from, e->to, e->label);
}

GraphFileFormat* setup_file_format_ubigraph(GraphFileFormat* gff)
{
    gff->write_header = 0;
    gff->write_pre_node = 0;
    gff->write_mid_node = &write_ubigraph_node;
    gff->write_edge = &write_ubigraph_edge;
    gff->write_post_node = 0;
    gff->write_footer = 0;
    gff->can_handle_missing_nodes = true;    
    return gff;
}

/*----------------------------------------------------------------------*
 * GraphML                                                              *
 *----------------------------------------------------------------------*/
void write_graphml_header(FILE *fp, GraphToolsOptions* options)
{
    fprintf(fp, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(fp, "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\" ");
    fprintf(fp, "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" ");
    fprintf(fp, "xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns ");
    fprintf(fp, "http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n");
    fprintf(fp, "<graph id=\"graphout\" edgedefault=\"directed\">\n");
}

void write_graphml_footer(FILE* fp, GraphToolsOptions* options)
{
    fprintf(fp, "</graph>\n");
    fprintf(fp, "</graphml>\n");
}

void write_graphml_node(FILE* fp, GraphoutNode* n, GraphToolsOptions* options)
{
    fprintf(fp, "<node id=\"%s\"/>\n", n->seq);
}

void write_graphml_edge(FILE* fp, GraphoutEdge* e, GraphToolsOptions* options)
{
    fprintf(fp, "<edge id=\"e%d\" source=\"%s\" target=\"%s\" label=\"%s\"/>\n", e->id, e->from, e->to, e->label);
}

GraphFileFormat* setup_file_format_graphml(GraphFileFormat* gff)
{
    gff->write_header = &write_graphml_header;
    gff->write_pre_node = &write_graphml_node;
    gff->write_mid_node = 0;
    gff->write_edge = &write_graphml_edge;
    gff->write_post_node = 0;
    gff->write_footer = &write_graphml_footer;
    gff->can_handle_missing_nodes = true;
    return gff;
}

/*----------------------------------------------------------------------*
 * GEXF                                                                 *
 *----------------------------------------------------------------------*/
void write_gexf_header(FILE *fp, GraphToolsOptions* options)
{    
    fprintf(fp, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(fp, "<gexf xmlns=\"http://www.gexf.net/1.2draft\" ");
    fprintf(fp, "xmlns:viz=\"http://www.gexf.net/1.2draft/viz\" ");
    fprintf(fp, "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" ");
    fprintf(fp, "xsi:schemaLocation=\"http://www.gexf.net/1.2draft ");
    fprintf(fp, "http://www.gexf.net/1.2draft/gexf.xsd\" version=\"1.2\">\n");
    fprintf(fp, "<graph mode=\"static\" defaultedgetype=\"directed\">\n");
}

void write_gexf_footer(FILE* fp, GraphToolsOptions* options)
{
    fprintf(fp, "</edges>\n");
    fprintf(fp, "</graph>\n");
    fprintf(fp, "</gexf>\n");
}

void write_gexf_node(FILE* fp, GraphoutNode* n, GraphToolsOptions* options)
{
    if (n->id==0) {
        fprintf(fp, "<nodes>\n");
    }
    
    fprintf(fp, "<node id=\"%s\" label=\"%s\">\n", n->seq, n->seq);
    
    if ((n->node->coverage[0] > 0) && (n->node->coverage[1] > 0)) {
        fprintf(fp, "<viz:color r=\"0\" g=\"255\" b=\"0\" a=\"0.5\"/>\n");
    } else if (n->node->coverage[0] > 0) {
        fprintf(fp, "<viz:color r=\"255\" g=\"0\" b=\"0\" a=\"0.5\"/>\n");
    } else if (n->node->coverage[1] > 0) {
        fprintf(fp, "<viz:color r=\"0\" g=\"0\" b=\"255\" a=\"0.5\"/>\n");
    } else {
        fprintf(fp, "<viz:color r=\"128\" g=\"128\" b=\"128\" a=\"0.5\"/>\n");
    }
    
    fprintf(fp, "<viz:size value=\"20\"/>\n");

    fprintf(fp, "</node>\n");
}

void write_gexf_edge(FILE* fp, GraphoutEdge* e, GraphToolsOptions* options)
{
    if (e->id==0) {
        fprintf(fp, "</nodes>\n");
        fprintf(fp, "<edges>\n");
    }
    fprintf(fp, "<edge id=\"e%d\" source=\"%s\" target=\"%s\" label=\"%s\">\n", e->id, e->from, e->to, e->label);
    fprintf(fp, "<viz:color r=\"0\" g=\"0\" b=\"0\"/>\n");
    //fprintf(fp, "<viz:thickness value=\"3.0\"/>\n");
    //fprintf(fp, "<viz:shape value=\"solid\"/>\n");
    fprintf(fp, "</edge>\n");
}

GraphFileFormat* setup_file_format_gexf(GraphFileFormat* gff)
{
    gff->write_header = &write_gexf_header;
    gff->write_pre_node = &write_gexf_node;
    gff->write_mid_node = 0;
    gff->write_edge = &write_gexf_edge;
    gff->write_post_node = 0;
    gff->write_footer = &write_gexf_footer;
    gff->can_handle_missing_nodes = true;
    return gff;
}

/*----------------------------------------------------------------------*
 * GML                                                                 *
 *----------------------------------------------------------------------*/
void write_gml_header(FILE *fp, GraphToolsOptions* options)
{
    fprintf(fp, "graph [\n");
    fprintf(fp, "  comment \"Generated by graphout\"\n");
    fprintf(fp, "  directed 1\n");
    fprintf(fp, "  IsPlanar 1\n");
}

void write_gml_footer(FILE* fp, GraphToolsOptions* options)
{
    fprintf(fp, "]\n");
}

void write_gml_node(FILE* fp, GraphoutNode* n, GraphToolsOptions* options)
{
    fprintf(fp, "  node [\n");
    fprintf(fp, "    id %d\n", n->id);
    fprintf(fp, "    label \"%s\"\n", n->seq);
    fprintf(fp, "  ]\n");
}

void write_gml_edge(FILE* fp, GraphoutEdge* e, GraphToolsOptions* options)
{
    fprintf(fp, "  edge [\n");
    fprintf(fp, "    source %d\n", e->id_from);
    fprintf(fp, "    target %d\n", e->id_to);
    fprintf(fp, "    label \"%s\"\n", e->label);
    fprintf(fp, "  ]\n");
}

GraphFileFormat* setup_file_format_gml(GraphFileFormat* gff)
{
    gff->write_header = &write_gml_header;
    gff->write_pre_node = &write_gml_node;
    gff->write_mid_node = 0;
    gff->write_edge = &write_gml_edge;
    gff->write_post_node = 0;
    gff->write_footer = &write_gml_footer;
    gff->can_handle_missing_nodes = false;
    return gff;
}

