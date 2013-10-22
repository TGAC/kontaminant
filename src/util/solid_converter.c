#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <execinfo.h>
#include <signal.h>
#include <unistd.h>
#include <inttypes.h>
#include "binary_kmer.h"
#include "seq.h"
#include "colour_space_seq.h"

#define MAX_SEQ_LENGTH  20000

typedef enum {BASE_TO_SOLID, SOLID_TO_BASE} converter;
struct arguments{
    char * input_file;
    char * output_file;
    converter conv;
    
};

void parse_opt (int key, char *arg,struct arguments *arguments )
{
    switch(key){//TODO: write all the pertinent validations if we plan to extend this. 
        case 'h':
            fprintf(stderr, "Usage: solid_converted [-h] [-i INPUT_FILE] [-o OUTPUT_FILE] [-c SOLID_TO_BASE | BASE_TO_SOLID] \n");
            exit(0) ;
        case 'i':
            arguments->input_file = arg;
            break;
        case 'o':
            arguments->output_file = arg;
            break;
        case 'c':
            if(strcmp("SOLID_TO_BASE", arg) == 0){
                arguments->conv = SOLID_TO_BASE;
            }else if(strcmp("BASE_TO_SOLID", arg) == 0){
                arguments->conv = BASE_TO_SOLID;
            }else{
                fprintf(stderr, "Invalid option for -c [SOLID_TO_BASE | BASE_TO_SOLID]\n");
                assert(false);
            }
            break;
            
    }
    
}

void converter_parse(int argc, char **argv, struct arguments * arg){
    opterr = 0;
	int c;
	while ((c = getopt (argc, argv, "hi:o:c:")) != -1)
	{
		
		parse_opt (c, optarg, arg);
	}
	
	
}

int main(int argc, char **argv)
{
    FILE * in;
	FILE * out;
    boolean full_entry;
    int seq_len = MAX_SEQ_LENGTH ;//To make this parametizable in the future. 
    //Sequence * sequence_new(int max_read_length, int max_name_length, char offset);
    Sequence * seq_in = sequence_new(seq_len, 300, 0); //TODO: make this parametizable. 
    Sequence * seq_out = sequence_new(seq_len, 300, 0);
    struct arguments arg;
    arg.input_file = "-";
    arg.output_file = "-";
    arg.conv = BASE_TO_SOLID;
    
    converter_parse(argc, argv, &arg);
    
    if(strcmp(arg.input_file, "-") == 0){
        in = stdin;
    }else{
        in = fopen(arg.input_file, "r");
        if(in == NULL){
            fprintf(stderr, "Unable to open file %s to read.", arg.input_file);
        }
    }
    
    if(strcmp(arg.output_file, "-") == 0){
        out = stdout;
    }else{
        out = fopen(arg.output_file, "w");
        if(out == NULL){
            if(in != stdin){
                fclose(in);
            }
            fprintf(stderr, "Unable to open file %s for writing", arg.output_file);
        }
    }
    ////this routine can read long sequences (eg full chromosomes) , this is implemented by reading the sequence in chunks
    //int read_sequence_from_fasta(FILE * fp, Sequence * seq, int max_chunk_length, boolean new_entry, boolean * full_entry, int offset);
    
    while (read_sequence_from_fasta(in, seq_in, seq_len,  true,&full_entry ,0) > 0){
        if(arg.conv == BASE_TO_SOLID){
            base_sequence_to_cs_sequence(seq_in, seq_out);
        }else if(arg.conv == SOLID_TO_BASE){
            fflush(stdout);
            fprintf(stderr, "input seq:\n");
            sequence_print_fasta(stderr, seq_in);
            fprintf(stderr, "out seq:\n");
            fflush(stderr);
            fflush(stdout);
            cs_sequence_to_base_sequence(seq_in, seq_out);
        }
        sequence_set_name(seq_in->name, seq_out);//TODO: use accessor
        sequence_print_fasta(out, seq_out);
    }
    return 0;
}