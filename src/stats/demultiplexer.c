//
//  demultiplexer.c
//  Untitled
//
//  Created by Ricardo Ramirez-Gonzalez (TGAC) on 16/03/2012.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>
#include <locale.h>
#include <getopt.h>
#include <ctype.h>
#include <assert.h>
#include <err.h>
#include <global.h>

#include <flags.h>
#include <file_format.h>
#include <nucleotide.h>
#include <seq.h>
#include <seq_io.h>
#include <logger.h>

typedef struct
{
    char * input_reads;
    char * output_folder;
    FileFormat format;
    sequence_header_type header_type;
    int max_read_length;
    int max_name_length;
    
}DemultiplexerCmdLine;

static void demultiplexer_default_opts(DemultiplexerCmdLine * cmd_line){
    cmd_line->input_reads=NULL;
    cmd_line->output_folder=NULL;
    cmd_line->format = FASTQ;
    cmd_line->header_type = CASAVA_1_8;
    cmd_line->max_name_length = MAX_READ_LENGTH;
    cmd_line->max_read_length = MAX_READ_LENGTH;

}
static void print_help(){
    fprintf(stderr,  "demultiplexer [-h | --help ] [ -input | -i reads ] (-o | --output) output_prefix \n");
    fprintf(stderr, "\t(-o | --output) output_prefix The prefix of the files to demultiplex. The filenames will be suffixed with the TAG\n");
    fprintf(stderr, "\t[ -input | -i reads ]: a FASTQ file. If ommited, reads from stdin\n");
    fprintf(stderr, "\t[-h | --help ]: prints this help\n");
    
     exit(-2);
}


static DemultiplexerCmdLine parse_args(int argc, char ** argv){
    //TODO: Add Parser for the rest of the options. 
    DemultiplexerCmdLine cmd_line;
    demultiplexer_default_opts(&cmd_line);
    int opt;
    int longopt_index;
    
    static struct option long_options[] = {
        
        {"help", no_argument, NULL, 'h'},
        {"input", required_argument, NULL, 'i'},
        {"output", required_argument, NULL, 'o'},
        
    };
    while ((opt = getopt_long(argc, argv,"hi:o:", long_options, &longopt_index)) > 0){
        switch (opt) {
            case 'o':
                cmd_line.output_folder = optarg;
                
                break;
            case 'i':
                cmd_line.input_reads = optarg;
                if (access(optarg, R_OK) == -1) {
                    errx(1,"[-i | --input] filename [%s] cannot be accessed", optarg);
                }
                break;
            case 'h':
                print_help();
                break;
            default:
                break;
        }
    }
    if (cmd_line.output_folder == NULL) {
        errx(1, "(-o | --output) output_prefix is required ");
    }

    return cmd_line;
}

int main( int argc, char ** argv){
    log_and_screen_printf("Demultiplexer\n\n");
    log_and_screen_printf(SVN_VERSION);
	log_and_screen_printf(SVN_COMMIT_DATE);
	log_and_screen_printf("Compiled on %s at %s \n\n", __DATE__, __TIME__);
    if (argc < 2) {
        print_help();
    }
    log_write_timestamp(1);
    if (argc == 0) {
        print_help();
    }
    DemultiplexerCmdLine cmd = parse_args(argc, argv);
    FILE * in = stdin;
    if (cmd.input_reads != NULL) {
        in = fopen(cmd.input_reads, "r");
        if(in == NULL){
			log_and_screen_printf("Unable to open file %s\n", cmd.input_reads);
			exit(-1);
		}
    }
    Sequence * seq = sequence_new(cmd.max_read_length, cmd.max_name_length, 33);
    seq->header = new_sequence_header(CASAVA_1_8);
    header_function * f = (header_function *) seq->header;
    char * index = f->get_index(seq);
    size_t prefix_length = strlen(cmd.output_folder);
    char * output_file = calloc(prefix_length + MAX_FIELD_SIZE + 1, sizeof(char *));
    char * index_pointer = output_file + prefix_length;
    strcpy(output_file, cmd.output_folder);
    printf("prefix: %s\n", output_file);
    while (read_sequence_from_fastq(in, seq, cmd.max_read_length)) {
        strcpy(index_pointer, index);
//     printf("index: %s\n new output: %s\n", index, output_file);
        append_sequence(output_file, seq, FASTQ);
    }
    
    if (in != stdin) {
        fclose(in);
    }
    
    
    
    
    return 0;
}