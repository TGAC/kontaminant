/*----------------------------------------------------------------------*
 * File:    reformat_log.c                                              *
 *                                                                      *
 * Purpose: Test program for cortex_log_parse                           *
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

int main(int argc, char **argv)
{
    CortexLog cl;

    if (argc != 2) {
        printf("Syntax: reformat_log <filename>\n");
        exit(1);
    }
        
    parse_cortex_log(argv[1], &cl);
    
    return 0;
}

