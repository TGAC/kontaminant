/*----------------------------------------------------------------------*
 * File:    cortex_log_parse.c                                          *
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

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#define __USE_XOPEN
#include <time.h>
#include <string.h>
#include <locale.h>
#include <flags.h>
#include <binary_kmer.h>
#include <element.h>
#include <dB_graph.h>
#include <file_reader.h>
#include <perfect_path.h>
#include <branches.h>
#include <y_walk.h>
#include <logger.h>
#include "cortex_log_parse.h"

/*----------------------------------------------------------------------*
 * Macros                                                               *
 *----------------------------------------------------------------------*/
#define return_error(A) error_line=__LINE__; return A
#define get_line_or_error(A,B,C,D) if (!get_next_line(A, B, C)) { error_line=__LINE__; return D; }

/*----------------------------------------------------------------------*
 * Globals                                                              *
 *----------------------------------------------------------------------*/
int error_line = 0;

/*----------------------------------------------------------------------*
 * Function: line_begins_with                                           *
 * Purpose:  Check if a line begins with a given string.                *
 * Params:   line -> line to test                                       *
 *           b -> beginning string to test                              *
 * Returns:  true if line begins with string, false otherwise           *
 *----------------------------------------------------------------------*/
boolean line_begins_with(char* line, char* b)
{
    return strncmp(line, b, strlen(b)) == 0;
}

/*----------------------------------------------------------------------*
 * Function: clean_string                                               *
 * Purpose:  Remove characters less than ASCII space (32) from end of   *
 *           string.                                                    *
 * Params:   s -> string                                                *
 * Returns:  pointer to string passed in                                *
 *----------------------------------------------------------------------*/
char* clean_string(char* s)
{
    int i = strlen(s) - 1;
    
    while ((i >= 0) && (s[i] <= ' ')) {
        s[i--] = 0;
    }
    
    return s;
}

/*----------------------------------------------------------------------*
 * Function: get_next_line                                              *
 * Purpose:  Get next non-blank line fromfile.                          *
 * Params:   line -> buffer to put line                                 *
 *           line_size = size of buffer                                 *
 *           fp -> file handle                                          *
 * Returns:  true if got line, false if not                             *
 *----------------------------------------------------------------------*/
boolean get_next_line(char* line, int line_size, FILE* fp)
{
    boolean got_line = false;
    int i=0;
    char c;
    
    while (!feof(fp) && (!got_line)) {
        while (!feof(fp)) {
            c = fgetc(fp);
            if (c >= ' ') {
                line[i++] = c;
            } else {
                break;
            }
        }
        line[i] = 0;
                    
        clean_string(line);
        if (strlen(line) > 0) {
            got_line = true;
        }
    }
    
    return got_line;
}

/*----------------------------------------------------------------------*
 * Function: parse_date                                                 *
 * Purpose:  Convert Cortex log date string to a time_t                 *
 * Params:   s -> date string                                           *
 *           d -> time_t date                                           *
 * Returns:  0 for success, else error code                             *
 *----------------------------------------------------------------------*/
int parse_date(char* s, time_t* d)
{
    struct tm tm;

    if (strptime(s, "%a %b %d %H:%M:%S %Y", &tm) == NULL) {
        *d = 0;
        return_error(ERROR_DATE);
    }
    
    tm.tm_isdst = -1;
    *d = mktime(&tm);
    
    return 0;
}

/*----------------------------------------------------------------------*
 * Function: parse_date_separator                                       *
 * Purpose:  A date separator is either ----- followed by a newline and *
 *           the date, or ----- Date -----                              *
 * Params:   fp -> file handle                                          *
 *           d -> time_t date                                           *
 * Returns:  0 for success, else error code                             *
 *----------------------------------------------------------------------*/
int parse_date_separator(FILE* fp, time_t* d)
{
    char line[1024];
    int  rc = 0;
    
    get_line_or_error(line, 1024, fp, ERROR_DATE);
    if (!(line_begins_with(line, "-----"))) {
        return_error(ERROR_DATE);
    }
    
    if (strlen(line) > 10) {
        rc = parse_date(line+6, d);
    } else {
        get_line_or_error(line, 1024, fp, ERROR_DATE);
        rc = parse_date(line, d);
    }
    
    return rc;
}

/*----------------------------------------------------------------------*
 * Function: is_digit                                                   *
 * Purpose:  Check if character is a digit (0-9) or not.                *
 * Params:   c = character to check                                     *
 * Returns:  true if digit, false if other character                    *
 *----------------------------------------------------------------------*/
boolean is_digit(char c)
{
    if ((c >= '0') && (c <= '9')) {
        return true;
    } else {
        return false;
    }
}

/*----------------------------------------------------------------------*
 * Function: parse_int                                                  *
 * Purpose:  Read an integer value (perhaps separated by commas) from   *
 *           string.                                                    *
 * Params:   s -> pointer to string (may be preceeded with spaces)      *
 * Returns:  the integer                                                *
 *----------------------------------------------------------------------*/
long int parse_int(char *s)
{
    char is[64];
    int i=0;
    int j=0;
    
    /* Get rid of preceeding spaces */
    while ((i < 64) && (s[i] == ' ')) {
        i++;
    }
    
    while (is_digit(s[i]) || s[i]==',') {
        if (s[i] != ',') {
            is[j++] = s[i];
        }
        i++;
    }
    is[j] = 0;
    
    return atol(is);
}

/*----------------------------------------------------------------------*
 * Function: parse_header                                               *
 * Purpose:  Parse header section of log                                *
 * Params:   fp -> file handle                                          *
 *           cl -> CortexLog structure to populate                      *
 * Returns:  0 for no error, else error number                          *
 *----------------------------------------------------------------------*/
int parse_header(FILE* fp, CortexLog* cl)
{
    char line[1024];
    
    get_line_or_error(line, 1024, fp, ERROR_HEADER);
    
    if (!strstr(line, "Log started")) {
        return_error(ERROR_HEADER);
    }
    
    if (parse_date(line, &cl->start_time) != 0) {
        return_error(ERROR_HEADER);
    }

    get_line_or_error(line, 1024, fp, ERROR_HEADER);

    if (strcmp(line, "Cortex Con") == 0) {
        cl->variant = CORTEX_CON;
    } else if (strstr(line, "Cortex Bub")) {
        cl->variant = CORTEX_BUB;
    } else if (strstr(line, "Cortex Var")) {
        cl->variant = CORTEX_VAR;
    } else {
        cl->variant = UNKNOWN;
    }

    return 0;
}

/*----------------------------------------------------------------------*
 * Function: parse_parameter_section                                    *
 * Purpose:  Parse parameter section of log                             *
 * Params:   fp -> file handle                                          *
 *           cl -> CortexLog structure to populate                      *
 * Returns:  0 for no error, else error number                          *
 *----------------------------------------------------------------------*/
int parse_parameter_section(FILE* fp, CortexLog* cl)
{
    char line[1024];
    long file_position = ftell(fp);
    
    while (get_next_line(line, 1024, fp) && (line[0] != '-')) {
        if (line_begins_with(line, "Max k:")) {
            cl->max_k = atoi(line + 6);
        } else if (line_begins_with(line, "Command:")) {
            if (!(cl->command_line = malloc(strlen(line)))) {
                return_error(ERROR_MEMORY);
            }
            strcpy(cl->command_line, line+9);
        } else if (line_begins_with(line, "Quality score offset:")) {
            cl->quality_score_offset = atoi(line + strlen("Quality score offset:"));
        } else if (line_begins_with(line, "Quality score threshold:")) {
            cl->quality_score_threshold = atoi(line + strlen("Quality score threshold:"));
        }
        
        file_position = ftell(fp);
    }
    
    fseek(fp, file_position, SEEK_SET);
    
    return 0;
}

/*----------------------------------------------------------------------*
 * Function: parse_foff_sectiokn                                        *
 * Purpose:  Parse file of files section of log                         *
 * Params:   fp -> file handle                                          *
 *           cl -> CortexLog structure to populate                      *
 * Returns:  0 for no error, else error number                          *
 *----------------------------------------------------------------------*/
int parse_foff_section(FILE* fp, CortexLog* cl)
{
    char line[1024];
    long file_position = ftell(fp);
    
    if (parse_date_separator(fp, &cl->fof_time) != 0) {
        return_error(ERROR_FOFF);
    }
    
    while (get_next_line(line, 1024, fp) && (line[0] != '-')) {
        if (line_begins_with(line, "Input file of filenames: ")) {
            if (!(cl->file_of_filenames = malloc(strlen(line)))) {
                return_error(ERROR_MEMORY);
            }
            strcpy(cl->file_of_filenames, line+strlen("Input file of filenames: "));
        } else if (line_begins_with(line, "Kmer size:")) {
            char* p;
            cl->kmer_size = parse_int(line+strlen("Kmer size:"));
            p = strstr(line, "Hash table size (");
            if (!p) {
                return_error(ERROR_FOFF);
            }
            
            cl->n = parse_int(p+strlen("Hash table size ("));
            p = strstr(line, "Hash table bucket size:");
            if (!p) {
                return_error(ERROR_FOFF);
            }
            cl->b = parse_int(p+strlen("Hash table bucket size:"));

            file_position = ftell(fp);
        }
    }
    
    fseek(fp, file_position, SEEK_SET);
    
    return 0;
}

/*----------------------------------------------------------------------*
 * Function: parse_hash_table_info                                      *
 * Purpose:  Parse hash table information section of log                *
 * Params:   fp -> file handle                                          *
 *           cl -> CortexLog structure to populate                      *
 * Returns:  0 for no error, else error number                          *
 *----------------------------------------------------------------------*/
int parse_hash_table_info(FILE* fp, HashTableInfo* hti)
{
    char line[1024];
    long last_line;
    boolean keep_going = true;

    memset(hti, 0, sizeof(HashTableInfo));

    while (keep_going && (!feof(fp))) {
        char* p;

        last_line = ftell(fp);
        
        get_line_or_error(line, 1024, fp, ERROR_HASHTABLE);
        
        if (line_begins_with(line, "dBGraph:")) {
            /* Do nothing */
        } else if (line_begins_with(line, "Hash:")) {
        } else if ((p=strstr(line, "unique kmers:")) > 0) {
            hti->unique_kmers = parse_int(p+strlen("unique kmers:"));
        } else if ((p=strstr(line, "Capacity:")) > 0) {
            hti->capacity = parse_int(p+strlen("Capacity:"));
        } else if ((p=strstr(line, "Occupied:")) > 0) {
            hti->occupied_pc = atof(p+strlen("Occupied:"));
        } else if (strstr(line, "Collisions:")) {
            /* Do nothing */
        } else if ((p=strstr(line, "tries ")) > 0) {
            int n = parse_int(p+strlen("tries "));
            int t;
            if ((p=strstr(line, ":")) > 0) {
                t = parse_int(p+1);
                hti->tries[n] = t;
            }
        } else {
            fseek(fp, last_line, SEEK_SET);
            keep_going = false;
        }
    }
    
    return 0;
}

/*----------------------------------------------------------------------*
 * Function: parse_hash_table_created                                   *
 * Purpose:  Parse hash table creation section of log                   *
 * Params:   fp -> file handle                                          *
 *           cl -> CortexLog structure to populate                      *
 * Returns:  0 for no error, else error number                          *
 *----------------------------------------------------------------------*/
int parse_hash_table_created(FILE* fp, CortexLog* cl)
{
    char line[1024];
    HashTableInfo hti;
    
    get_line_or_error(line, 1024, fp, ERROR_HASHTABLE);

    if (!line_begins_with(line, "Table created")) {
        return_error(ERROR_HASHTABLE);
    }
    parse_hash_table_info(fp, &hti);
    
    return 0;
}

/*----------------------------------------------------------------------*
 * Function: parse_file                                                 *
 * Purpose:  Parse file section of log                                  *
 * Params:   fp -> file handle                                          *
 *           cl -> CortexLog structure to populate                      *
 * Returns:  0 for no error, else error number                          *
 *----------------------------------------------------------------------*/
int parse_file(FILE* fp, CortexLog* cl)
{
    char line[1024];
    char* p;
    
    get_line_or_error(line, 1024, fp, ERROR_FILELOAD);

    if (!line_begins_with(line, "Reading")) {
        return_error(ERROR_FILELOAD);
    }

    if (!(cl->input_file[cl->number_of_files] = malloc(strlen(line)))) {
        return_error(ERROR_MEMORY);
    }

    if (!(p=strchr(line, ':'))) {
        return_error(ERROR_FILELOAD);
    }

    if (!(cl->input_file[cl->number_of_files] = malloc(sizeof(CortexFileInfo)))) {     
        return_error(ERROR_MEMORY);
    }
    
    if (!(cl->input_file[cl->number_of_files]->filename = malloc(strlen(line)))) {
        return_error(ERROR_MEMORY);
    }

    strcpy(cl->input_file[cl->number_of_files]->filename, p+2);
    
    if (parse_date_separator(fp, &cl->input_file[cl->number_of_files]->time_read_complete) != 0) {
        return ERROR_FILELOAD;
    }
    
    get_line_or_error(line, 1024, fp, ERROR_FILELOAD);

    if (!line_begins_with(line, "Read of file")) {
        return_error(ERROR_FILELOAD);
    }
    
    /* Read of file 2 complete. Total kmers: 8,519,229 Bad reads: 1,998,510 Seq length: 255,681,576 Total seq length: 383,522,364 */
    if (!(p=strstr(line, "Bad reads: "))) {
        return_error(ERROR_FILELOAD);
    }
    cl->input_file[cl->number_of_files]->bad_reads = parse_int(p+strlen("Bad reads: "));

    if (!(p=strstr(line, "Seq length: "))) {
        return_error(ERROR_FILELOAD);
    }
    cl->input_file[cl->number_of_files]->seq_length = parse_int(p+strlen("Seq length: "));

    if (!(p=strstr(line, "Total seq length: "))) {
        return_error(ERROR_FILELOAD);
    }
    cl->input_file[cl->number_of_files]->total_seq_length = parse_int(p+strlen("Total seq length: "));
        
    parse_hash_table_info(fp, &cl->input_file[cl->number_of_files]->hash_table_post_load);
    
    cl->number_of_files++;
    
    return 0;
}

/*----------------------------------------------------------------------*
 * Function: parse_file_loads                                           *
 * Purpose:  Parse file load section of log                             *
 * Params:   fp -> file handle                                          *
 *           cl -> CortexLog structure to populate                      *
 * Returns:  0 for no error, else error number                          *
 *----------------------------------------------------------------------*/
int parse_file_loads(FILE* fp, CortexLog* cl)
{
    long file_position;
    char line[1024];
    boolean got_one = false;

    do {
        got_one = false;
        file_position = ftell(fp);
        get_line_or_error(line, 1024, fp, ERROR_FILELOAD);        
        fseek(fp, file_position, SEEK_SET);
        if (line_begins_with(line, "Reading")) {
            parse_file(fp, cl);
            got_one = true;
        }
    } while (got_one);

    return 0;
}

/*----------------------------------------------------------------------*
 * Function: parse_clip_tips                                            *
 * Purpose:  Parse clip tip section of log                                  *
 * Params:   fp -> file handle                                          *
 *           cl -> CortexLog structure to populate                      *
 * Returns:  0 for no error, else error number                          *
 *----------------------------------------------------------------------*/
int parse_clip_tips(FILE* fp, CortexLog* cl)
{
    char line[1024];
    
    get_line_or_error(line, 1024, fp, ERROR_CLIPTIPS);
    if (!line_begins_with(line, "Clip tips")) {
        return_error(ERROR_CLIPTIPS);
    }
    
    cl->clip_tips.enabled = true;
    get_line_or_error(line, 1024, fp, ERROR_CLIPTIPS);
    if (!line_begins_with(line, "clip_tip removed nodes ")) {
        return_error(ERROR_CLIPTIPS);
    }
    cl->clip_tips.nodes_removed = parse_int(line + strlen("clip_tip removed nodes "));

    get_line_or_error(line, 1024, fp, ERROR_CLIPTIPS);
    if (!strstr(line, "tips clipped")) {
        return_error(ERROR_CLIPTIPS);
    }
    cl->clip_tips.tips_clipped = parse_int(line);
    
    return 0;
}

/*----------------------------------------------------------------------*
 * Function: parse_remove_bubbles                                       *
 * Purpose:  Parse remove bubbles section of log                        *
 * Params:   fp -> file handle                                          *
 *           cl -> CortexLog structure to populate                      *
 * Returns:  0 for no error, else error number                          *
 *----------------------------------------------------------------------*/
int parse_remove_bubbles(FILE* fp, CortexLog* cl)
{
    char line[1024];
    
    get_line_or_error(line, 1024, fp, ERROR_REMOVEBUBBLES);
    if (!line_begins_with(line, "Removing bubbles")) {
        return_error(ERROR_REMOVEBUBBLES);
    }

    cl->remove_bubbles.enabled = true;

    get_line_or_error(line, 1024, fp, ERROR_REMOVEBUBBLES);
    if (!strstr(line, "bubbles removed")) {
        return_error(ERROR_REMOVEBUBBLES);
    }
    cl->remove_bubbles.bubbles_removed = parse_int(line);

    get_line_or_error(line, 1024, fp, ERROR_REMOVEBUBBLES);
    if (!strstr(line, "kmers removed")) {
        return_error(ERROR_REMOVEBUBBLES);
    }
    cl->remove_bubbles.nodes_removed = parse_int(line);
        
    return 0;
}    

/*----------------------------------------------------------------------*
 * Function: parse_low_coverage_paths                                   *
 * Purpose:  Parse low coverage path cleaning section of log            *
 * Params:   fp -> file handle                                          *
 *           cl -> CortexLog structure to populate                      *
 * Returns:  0 for no error, else error number                          *
 *----------------------------------------------------------------------*/
int parse_low_coverage_paths(FILE* fp, CortexLog* cl)
{
    char line[1024];
    
    get_line_or_error(line, 1024, fp, ERROR_LOWCOVERAGEPATHS);
    if (!line_begins_with(line, "Removing paths with low coverage")) {
        return_error(ERROR_LOWCOVERAGEPATHS);
    }
    
    cl->remove_low_coverage_paths.enabled = true;

    get_line_or_error(line, 1024, fp, ERROR_LOWCOVERAGEPATHS);
    if (!strstr(line, "nodes removed")) {
        return_error(ERROR_LOWCOVERAGEPATHS);
    }
    
    cl->remove_low_coverage_paths.nodes_removed = parse_int(line);
    
    return 0;
}

/*----------------------------------------------------------------------*
 * Function: parse_low_coverage_nodes                                   *
 * Purpose:  Parse low coverage node cleaning section of log            *
 * Params:   fp -> file handle                                          *
 *           cl -> CortexLog structure to populate                      *
 * Returns:  0 for no error, else error number                          *
 *----------------------------------------------------------------------*/
int parse_low_coverage_nodes(FILE* fp, CortexLog* cl)
{
    char line[1024];
    char* p;
    
    get_line_or_error(line, 1024, fp, ERROR_LOWCOVERAGENODES);
    if (!line_begins_with(line, "Remove low coverage nodes")) {
        return_error(ERROR_LOWCOVERAGENODES);
    }    
    
    cl->remove_low_coverage_nodes.enabled = true;
    
    get_line_or_error(line, 1024, fp, ERROR_LOWCOVERAGENODES);
    if (!(p=strstr(line, "Removed "))) {
        return_error(ERROR_LOWCOVERAGENODES);
    }
        
    cl->remove_low_coverage_nodes.nodes_removed = parse_int(p + strlen("Removed "));

    return 0;
}

/*----------------------------------------------------------------------*
 * Function: parse_dump_supernodes                                      *
 * Purpose:  Parse contig dumping section of log                        *
 * Params:   fp -> file handle                                          *
 *           cl -> CortexLog structure to populate                      *
 * Returns:  0 for no error, else error number                          *
 *----------------------------------------------------------------------*/
int parse_dump_supernodes(FILE* fp, CortexLog* cl)
{
    char line[1024];
    char *p;

    get_line_or_error(line, 1024, fp, ERROR_DUMPSUPERNODES);
    if (!line_begins_with(line, "Dumping supernodes: ")) {
        return_error(ERROR_DUMPSUPERNODES);
    }
    cl->dump_contigs.enabled = true;
    
    if (!(cl->dump_contigs.filename = malloc(strlen(line)))) {
        return_error(ERROR_MEMORY);
    }
    strcpy(cl->dump_contigs.filename, line+strlen("Dumping supernodes: "));
    
    get_line_or_error(line, 1024, fp, ERROR_DUMPSUPERNODES);
    if (!(p=strstr(line, "nodes visited ["))) {
        return_error(ERROR_DUMPSUPERNODES);
    }

    cl->dump_contigs.nodes_visited = parse_int(line);
    cl->dump_contigs.singletons = parse_int(p + strlen("nodes visited ["));

    return 0;
}

/*----------------------------------------------------------------------*
 * Function: parse_dump_graph                                           *
 * Purpose:  Parse dump graph section of log                            *
 * Params:   fp -> file handle                                          *
 *           cl -> CortexLog structure to populate                      *
 * Returns:  0 for no error, else error number                          *
 *----------------------------------------------------------------------*/
int parse_dump_graph(FILE* fp, CortexLog* cl)
{
    char line[1024];
    long file_position=0;

    get_line_or_error(line, 1024, fp, ERROR_DUMPGRAPH);
    if (!line_begins_with(line, "Dumping graph: ")) {
        return_error(ERROR_DUMPGRAPH);
    }
    cl->dump_graph.enabled = true;

    if (!(cl->dump_graph.filename = malloc(strlen(line)))) {
        return_error(ERROR_MEMORY);
    }
    strcpy(cl->dump_graph.filename, line+strlen("Dumping graph: "));    

    file_position = ftell(fp);
    
    get_line_or_error(line, 1024, fp, ERROR_DUMPGRAPH);
    if (strstr(line, "kmers dumped")) {
        cl->dump_graph.kmers_dumped = parse_int(line);
    } else {
        fseek(fp, file_position, SEEK_SET);
    }
        
    return 0;
}
    
/*----------------------------------------------------------------------*
 * Function: parse_find_bubbles                                         *
 * Purpose:  Parse the bubble finding part of the log                   *
 * Params:   fp -> file handle                                          *
 *           cl -> CortexLog structure to populate                      *
 * Returns:  0 for no error, else error number                          *
 *----------------------------------------------------------------------*/
int parse_find_bubbles(FILE* fp, CortexLog* cl)
{
    char line[1024];
    char* p;
   
    get_line_or_error(line, 1024, fp, ERROR_FINDBUBBLES);

    if (!((line_begins_with(line, "Finding bubbles to output...")) ||
        (line_begins_with(line, "Looking for branches to walk.")))) {
        return_error(ERROR_FINDBUBBLES);
    }
    
    cl->find_bubbles.enabled = true;

    get_line_or_error(line, 1024, fp, ERROR_FINDBUBBLES);
    if (!(p=strstr(line, "Finished walking - "))) {
        return_error(ERROR_FINDBUBBLES);
    }
    cl->find_bubbles.contigs_output = atoi(p + strlen("Finished walking - "));
    
    return 0;
}

/*----------------------------------------------------------------------*
 * Function: ignore_section                                             *
 * Purpose:  Ignore a section of log we don't understand                *
 * Params:   fp -> file handle                                          *
 *           cl -> CortexLog structure to populate                      *
 * Returns:  0 for no error, else error number                          *
 *----------------------------------------------------------------------*/
int ignore_section(FILE* fp)
{
    boolean got_new_section = false;
    long file_position;
    char line[1024];
    
    while ((!got_new_section) && (!feof(fp))) {
        file_position = ftell(fp);
        get_line_or_error(line, 1024, fp, ERROR_OPTIONS);
        if (strncmp(line, "-----", 5) == 0) {
            fseek(fp, file_position, SEEK_SET);
            got_new_section = true;
        }
    }
    
    return 0;
}

/*----------------------------------------------------------------------*
 * Function: parse_options                                              *
 * Purpose:  Parse the optional bits of the log - each one begins with  *
 *           '-----' then the date.                                     *
 * Params:   fp -> file handle                                          *
 *           cl -> CortexLog structure to populate                      *
 * Returns:  0 for no error, else error number                          *
 *----------------------------------------------------------------------*/
int parse_options(FILE* fp, CortexLog* cl)
{
    long file_position;
    char line[1024];
    time_t t;
    boolean done = false;
    int rc = 0;
    
    while ((!done) && (rc == 0) &&(!feof(fp))) {
        if (parse_date_separator(fp, &t) != 0) {
            return_error(ERROR_OPTIONS);
        }
        
        file_position = ftell(fp);
        get_line_or_error(line, 1024, fp, ERROR_OPTIONS);
        
        fseek(fp, file_position, SEEK_SET);
        if (line_begins_with(line, "Clip tips")) {
            rc = parse_clip_tips(fp, cl);
        } else if (line_begins_with(line, "Removing bubbles")) {
            rc = parse_remove_bubbles(fp, cl);
        } else if (line_begins_with(line, "Dumping supernodes")) {
            rc = parse_dump_supernodes(fp, cl);
        } else if (line_begins_with(line, "Removing paths with low coverage")) {
            rc = parse_low_coverage_paths(fp, cl);
        } else if (line_begins_with(line, "Remove low coverage nodes")) {
            rc = parse_low_coverage_nodes(fp, cl);
        } else if (line_begins_with(line, "Dumping graph")) {
            rc = parse_dump_graph(fp, cl);
        } else if ((line_begins_with(line, "Finding bubbles to output...")) ||
                    (line_begins_with(line, "Looking for branches to walk."))) {
            rc = parse_find_bubbles(fp, cl);
        } else if (line_begins_with(line, "DONE")) {
            cl->end_time = t;            
            done = true;
        } else {
            printf("Ignoring section [%s]\n", line);
            ignore_section(fp);
        }
    }
    
    return rc;
}

/*----------------------------------------------------------------------*
 * Function: print_cortex_log                                           *
 * Purpose:  Print log to structure to screen                           *
 * Params:   fp -> file handle                                          *
 *           cl -> CortexLog structure to populate                      *
 * Returns:  0 for no error, else error number                          *
 *----------------------------------------------------------------------*/
void print_cortex_log(CortexLog* cl)
{
    int i,j;
    
    printf("            Log started: %s", ctime(&cl->start_time));
    printf("                Variant: ");
    if (cl->variant == CORTEX_CON) {
        printf("Cortex Con\n");
    } else if (cl->variant == CORTEX_BUB) {
        printf("Cortex Bub\n");
    } else if (cl->variant == CORTEX_VAR) {
        printf("Cortex Var\n");
    } else {
        printf("Unknown\n");
    }
    printf("                  Max k: %d\n", cl->max_k);
    printf("           Command line: %s\n", cl->command_line);
    printf("   Quality score offset: %d\n", cl->quality_score_offset);
    printf("Quality score threshold: %d\n", cl->quality_score_threshold);
    printf("Input file of filenames: %s\n", cl->file_of_filenames);
    printf("              Kmer size: %d\n", cl->kmer_size);
    printf("                      n: %d\n", cl->n);
    printf("                      b: %d\n", cl->b);

    for (i=0; i<cl->number_of_files; i++) {
        printf("                   File: %s\n", cl->input_file[i]->filename);
        printf("              Bad reads: %lld\n", cl->input_file[i]->bad_reads);
        printf("             Seq length: %lld\n", cl->input_file[i]->seq_length);
        printf("       Total seq length: %lld\n", cl->input_file[i]->total_seq_length);
        printf(" Hash table kmers after: %lld\n", cl->input_file[i]->hash_table_post_load.unique_kmers);
        for (j=0; j<MAX_HASH_TRIES; j++) {
            if (cl->input_file[i]->hash_table_post_load.tries[j] > 0) {
                printf("                  Tries: %d = %d\n", j, cl->input_file[i]->hash_table_post_load.tries[j]);
            }
        }
    }
    
    if (cl->remove_low_coverage_nodes.enabled) {
        printf(" Low cov. nodes removed: %lld\n", cl->remove_low_coverage_nodes.nodes_removed);
    }
    
    if (cl->clip_tips.enabled) {
        printf("           Tips clipped: %lld\n", cl->clip_tips.tips_clipped);
        printf("Clip tips nodes removed: %lld\n", cl->clip_tips.nodes_removed);
    }

    if (cl->remove_low_coverage_paths.enabled) {
        printf("   Low cov. paths nodes: %lld\n", cl->remove_low_coverage_paths.nodes_removed);
    }
    
    if (cl->remove_bubbles.enabled) {
        printf("        Bubbles removed: %lld\n", cl->remove_bubbles.bubbles_removed);
        printf("  Bubbles nodes removed: %lld\n", cl->remove_bubbles.nodes_removed);
    }

    if (cl->dump_graph.enabled) {
        printf("    Dump graph filename: %s\n", cl->dump_graph.filename);
    }
    
    if (cl->dump_contigs.enabled) {
        printf(" Output contig filename: %s\n", cl->dump_contigs.filename);
        printf("          Nodes visited: %lld\n", cl->dump_contigs.nodes_visited);
        printf("             Singletons: %lld\n", cl->dump_contigs.singletons);
    }

    printf("              Log ended: %s", ctime(&cl->end_time));
    printf("         Execution time: %d mins\n", (int)(difftime(cl->end_time, cl->start_time) / 60));
}

/*----------------------------------------------------------------------*
 * Function: parse_cortex_log                                           *
 * Purpose:  Parse log file and build CortexLog structure               *
 * Params:   fp -> file handle                                          *
 *           cl -> CortexLog structure to populate                      *
 * Returns:  0 for no error, else error number                          *
 *----------------------------------------------------------------------*/
int parse_cortex_log(char* filename, CortexLog* cl)
{
    FILE* fp = fopen(filename, "r");
    int   rc = 0;

    memset(cl, 0, sizeof(CortexLog));
    
    if (!fp) {
        return_error(ERROR_CANTOPENFILE);
    }
    
    if (!feof(fp)) {
        rc = parse_header(fp, cl);
    }
    
    if (!feof(fp) && (rc==0)) {
        rc = parse_parameter_section(fp, cl);
    }
    
    if (!feof(fp) && (rc==0)) {
        rc = parse_foff_section(fp, cl);
    }
    
    if (!feof(fp) && (rc==0)) {        
        time_t t;
        rc = parse_date_separator(fp, &t);
    }
    
    if (!feof(fp) && (rc==0)) {
        rc = parse_hash_table_created(fp, cl);
    }
    
    if (!feof(fp) && (rc==0)) {
        rc = parse_file_loads(fp, cl);
    }

    if (!feof(fp) && (rc==0)) {
        rc = parse_options(fp, cl);
    }
    
    fclose(fp);

    if (rc != 0) {
        printf("\nEncountered error while reading log: Error code %d, source line %d\n\n", rc, error_line);
    }
    
    print_cortex_log(cl);
    
    return rc;
}
