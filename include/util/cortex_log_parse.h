/*----------------------------------------------------------------------*
 * File:    cortex_log_parse.h                                          *
 *                                                                      *
 * Purpose: Parse log files output by cortex con                        *
 *                                                                      *
 * Author:  Richard Leggett                                             *
 *          The Sainsbury Laboratory                                    *
 *          Norwich Research Park, Colney, Norwich, NR4 7UH, UK         *
 *          richard.leggett@tsl.ac.uk                                   * 
 *                                                                      *
 * History: 06-Apr-11: RML: Created                                     *
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 * Error codes                                                          *
 *----------------------------------------------------------------------*/
#define ERROR_CANTOPENFILE     1
#define ERROR_HEADER           2
#define ERROR_DATE             3
#define ERROR_FOFF             4
#define ERROR_MEMORY           5
#define ERROR_HASHTABLE        6
#define ERROR_FORMAT           7
#define ERROR_FILETOOSHORT     8
#define ERROR_FILELOAD         9
#define ERROR_OPTIONS          10
#define ERROR_CLIPTIPS         11
#define ERROR_REMOVEBUBBLES    12
#define ERROR_DUMPSUPERNODES   13
#define ERROR_DUMPGRAPH        14
#define ERROR_LOWCOVERAGEPATHS 15
#define ERROR_LOWCOVERAGENODES 16
#define ERROR_FINDBUBBLES      17

/*----------------------------------------------------------------------*
 * Other constants                                                      *
 *----------------------------------------------------------------------*/
#define MAX_SOURCE_FILES 32
#define MAX_HASH_TRIES 32

/*----------------------------------------------------------------------*
 * Structures                                                           *
 *----------------------------------------------------------------------*/
typedef enum{
	UNKNOWN = 0,
	CORTEX_CON = 1,
	CORTEX_VAR = 2,
	CORTEX_BUB = 3
} CortexVariant;

typedef struct {
    int kmer_size;
    int n;
    int b;
    long long int capacity;
    long long int unique_kmers;
    double occupied_pc;
    int tries[MAX_HASH_TRIES];
} HashTableInfo;

typedef struct {
    char* filename;
    time_t time_read_complete;
    long long int bad_reads;
    long long int seq_length;
    long long int total_seq_length;
    HashTableInfo hash_table_post_load;
} CortexFileInfo;

typedef struct {
    boolean enabled;
    long long int tips_clipped;
    long long int nodes_removed;
} ClipTipsInfo;

typedef struct {
    boolean enabled;
    long long int nodes_removed;
} LowCoveragePathsInfo;

typedef struct {
    boolean enabled;
    long long int nodes_removed;
} LowCoverageNodesInfo;

typedef struct {
    boolean enabled;
    long long int bubbles_removed;
    long long int nodes_removed;
} RemoveBubblesInfo;

typedef struct {
    boolean enabled;
    long long int contigs_output;
} FindBubblesInfo;

typedef struct {
    boolean enabled;
    char* filename;
    long long int nodes_visited;
    long long int singletons;
} DumpContigsInfo;

typedef struct {
    boolean enabled;
    char* filename;
    int kmers_dumped;
} DumpGraphInfo;

typedef struct {
    CortexVariant variant;
    time_t start_time;
    char* command_line;
    int max_k;
    int quality_score_offset;
    int quality_score_threshold;
    time_t fof_time;
    char* file_of_filenames;
    int kmer_size;
    int n;
    int b;
    int number_of_files;
    CortexFileInfo* input_file[MAX_SOURCE_FILES];
    ClipTipsInfo clip_tips;
    LowCoveragePathsInfo remove_low_coverage_paths;
    LowCoverageNodesInfo remove_low_coverage_nodes;
    RemoveBubblesInfo remove_bubbles;
    FindBubblesInfo find_bubbles;
    DumpGraphInfo dump_graph;
    DumpContigsInfo dump_contigs;
    time_t end_time;
} CortexLog;

/*----------------------------------------------------------------------*
 * Function prototypes                                                  *
 *----------------------------------------------------------------------*/
int parse_cortex_log(char* filename, CortexLog* cl);
void print_cortex_log(CortexLog* cl);

