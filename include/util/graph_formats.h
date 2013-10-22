/*----------------------------------------------------------------------*
 * File:    graph_formats.h                                             *
 * Purpose: Graph file format functions for graphout                    *
 * Author:  Richard Leggett                                             *
 * History: 15-Feb-11: RML: Created                                     *
 *----------------------------------------------------------------------*/

#define MAX_KMER 95
#define MAX_NODES 10000

#define GRAPH_FORMAT_GRAPHVIZ  1
#define GRAPH_FORMAT_GRAPHML   2
#define GRAPH_FORMAT_BIOLAYOUT 4
#define GRAPH_FORMAT_UBIGRAPH  8
#define GRAPH_FORMAT_GEXF      16
#define GRAPH_FORMAT_GML       32

/*----------------------------------------------------------------------*
 * Structures                                                           *
 *----------------------------------------------------------------------*/
typedef struct {
    int id;
    char seq[MAX_KMER+1];
    dBNode* node;
} GraphoutNode;

typedef struct {
    int id;
    int id_from;
    int id_to;
    char from[MAX_KMER+1];
    char to[MAX_KMER+1];
    char label[MAX_KMER+1];
    char label_colour[128];
    Orientation orientation;
} GraphoutEdge;

typedef struct {
    void (* write_header)(FILE* fp, GraphToolsOptions* options);
    void (* write_pre_node)(FILE* fp, GraphoutNode* node, GraphToolsOptions* options);
    void (* write_mid_node)(FILE* fp, GraphoutNode* node, GraphToolsOptions* options);
    void (* write_edge)(FILE* fp, GraphoutEdge* edge, GraphToolsOptions* options);
    void (* write_post_node)(FILE* fp, GraphoutNode* node, GraphToolsOptions* options);
    void (* write_footer)(FILE* fp, GraphToolsOptions* options);
    boolean can_handle_missing_nodes;
} GraphFileFormat;

/*----------------------------------------------------------------------*
 * GraphViz (.dot)                                                      *
 *----------------------------------------------------------------------*/
void write_graphviz_header(FILE* fp, GraphToolsOptions* options);
void write_graphviz_footer(FILE* fp, GraphToolsOptions* options);
void write_graphviz_node(FILE* fp, GraphoutNode *n, GraphToolsOptions* options);
void write_graphviz_edge(FILE* fp, GraphoutEdge *e, GraphToolsOptions* options);
GraphFileFormat* setup_file_format_graphviz(GraphFileFormat* gff);

/*----------------------------------------------------------------------*
 * GraphML                                                              *
 *----------------------------------------------------------------------*/
void write_graphml_header(FILE *fp, GraphToolsOptions* options);
void write_graphml_footer(FILE* fp, GraphToolsOptions* options);
void write_graphml_node(FILE* fp, GraphoutNode* n, GraphToolsOptions* options);
void write_graphml_edge(FILE* fp, GraphoutEdge* e, GraphToolsOptions* options);
GraphFileFormat* setup_file_format_graphml(GraphFileFormat* gff);

/*----------------------------------------------------------------------*
 * BioLayout                                                            *
 *----------------------------------------------------------------------*/
void write_biolayout_edge(FILE* fp, GraphoutEdge *e, GraphToolsOptions* options);
void write_biolayout_node_class(FILE* fp, GraphoutNode *n, GraphToolsOptions* options);
GraphFileFormat* setup_file_format_biolayout(GraphFileFormat* gff);

/*----------------------------------------------------------------------*
 * Ubigraph                                                             *
 *----------------------------------------------------------------------*/
void write_ubigraph_node(FILE* fp, GraphoutNode* n, GraphToolsOptions* options);
void write_ubigraph_edge(FILE* fp, GraphoutEdge* e, GraphToolsOptions* options);
GraphFileFormat* setup_file_format_ubigraph(GraphFileFormat* gff);

/*----------------------------------------------------------------------*
 * GEXF                                                                 *
 *----------------------------------------------------------------------*/
void write_gexf_header(FILE *fp, GraphToolsOptions* options);
void write_gexf_footer(FILE* fp, GraphToolsOptions* options);
void write_gexf_node(FILE* fp, GraphoutNode* n, GraphToolsOptions* options);
void write_gexf_edge(FILE* fp, GraphoutEdge* e, GraphToolsOptions* options);
GraphFileFormat* setup_file_format_gexf(GraphFileFormat* gff);

/*----------------------------------------------------------------------*
 * GML                                                                  *
 *----------------------------------------------------------------------*/
void write_gml_header(FILE *fp, GraphToolsOptions* options);
void write_gml_footer(FILE* fp, GraphToolsOptions* options);
void write_gml_node(FILE* fp, GraphoutNode* n, GraphToolsOptions* options);
void write_gml_edge(FILE* fp, GraphoutEdge* e, GraphToolsOptions* options);
GraphFileFormat* setup_file_format_gml(GraphFileFormat* gff);
