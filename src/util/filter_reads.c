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
#include <global.h>
#include <flags.h>
#include <binary_kmer.h>
#include <element.h>
#include <dB_graph.h>
#include <file_reader.h>
#include <perfect_path.h>
#include <branches.h>
#include <y_walk.h>
#include <logger.h>
#include <assert.h>

#define MAX_LINE_LENGTH 100000

char* file_of_filenames = 0;
char* reference_filename = 0;
int kmer_size = 0;
int min_matches = 3; 
int hash_n=17;
int hash_b=60;

/*----------------------------------------------------------------------*
 * Function: parse_string                                               *
 * Purpose:  Return a string parameter (command line parsing)           *
 * Params:   argc = number of arguments                                 *
 *           argv -> argument array                                     *
 *           i -> pointer to argument number counter                    *
 * Returns:  a string and updates i                                     *
 *----------------------------------------------------------------------*/
char* parse_string(int argc, char* argv[], int* i)
{
    char* token = 0;
    
    if (strlen(argv[*i]) > 2)
        token = argv[*i] + 2;
    else if (*i < (argc-1))	{
        *i = *i + 1;
        token = argv[*i];
    }
    
    if ((token) && (token[0] == '-'))
        token = 0;
    
    return token;
}

/*----------------------------------------------------------------------*
 * Function: parse_int                                                  *
 * Purpose:  Return an integer parameter (command line parsing)         *
 * Params:   argc = number of arguments                                 *
 *           argv -> argument array                                     *
 *           i -> pointer to argument number counter                    *
 * Returns:  a number and updates i                                     *
 *----------------------------------------------------------------------*/
int parse_int(int argc, char* argv[], int* i)
{
    int v = -1;
    char* token = parse_string(argc, argv, i);
    
    if (token)
        v = atoi(token);
    
    return v;
}

/*----------------------------------------------------------------------*
 * Function: parse_command_line_args                                    *
 * Purpose:  Deal with command line arguments                           *
 * Params:   argc = number of arguments                                 *
 *           argv -> argument array                                     *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void parse_command_line_args(int argc, char* argv[])
{
    int i = 1;
    
    if (argc < 5)
    {
        printf("Syntax: filterreads [-i filename] [-r filename] [-k int] [options]\n");
        printf("where [-i filename] specifies the name of a file containing a list of read files\n");
        printf("      [-r int] specifies the name of the reference FASTA file\n");
        printf("      [-k int] specifies kmer size\n");
        printf("      [-m int] specifies minimum number of matching kmers per read (default 3)\n");
        printf("      [-n int] specifies number of buckets in hash table (default 17)\n");
        printf("      [-b int] specifies size of hash table buckets (default 60)\n");
        printf("\n");
        exit(1);
    }
    
    while(i < argc)
    {
        char* parameter = argv[i];
        if (parameter[0] == '-')
        {
            char option = parameter[1];
            
            switch (option)
            {
                case 'i':
                    file_of_filenames = parse_string(argc, argv, &i);
                    if (!file_of_filenames) {
                        printf("ERROR: Invalid filename for -i parameter.\n");
                        exit(1);
                    }
                    break;
                case 'r':
                    reference_filename = parse_string(argc, argv, &i);
                    if (!reference_filename) {
                        printf("ERROR: Invalid filename for -r parameter.\n");
                        exit(1);
                    }
                    break;
                case 'b':
                    hash_b = parse_int(argc, argv, &i);
                    break;
                case 'k':
                    kmer_size = parse_int(argc, argv, &i);
                    break;
                case 'n':
                    hash_n = parse_int(argc, argv, &i);
                    break;
                case 'm':
                    min_matches = parse_int(argc, argv, &i);
                    break;
                default:
                    printf("ERROR: Invalid option -%c\n", option);
                    exit(1);
            }
        }
        
        i++;
    }
    
    if ((kmer_size < 1) || (kmer_size > 255)) {
        printf("ERROR: Invalid kmer size.\n");
        exit(1);
    }
    
    if (min_matches < 1) {
        printf("ERROR: Invalid minimum matches.\n");
        exit(1);        
    }
    
    if (!file_of_filenames) {
        printf("ERROR: You must specify an input file of files.\n");
        exit(1);
    }

    if (!reference_filename) {
        printf("ERROR: You must specify a reference file.\n");
        exit(1);
    }

}

/*----------------------------------------------------------------------*
 * Function: chomp                                                      *
 * Purpose:  Removing trailing unprintable characters                   *
 * Params:   s -> string to work on                                     *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void chomp(char* s) {
    char* ptr = s + strlen(s)-1;
    
    while (*ptr < ' ') {
        *ptr-- = 0;
    }
}

/*----------------------------------------------------------------------*
 * Function: count_kmers                                                *
 * Purpose:  Go through all kmers in a read and see how many are in the *
             hash table.                                                *
 * Params:   read -> read to check                                      *
 *           db_graph -> hash table                                     *
 * Returns:  Number of matching kmers (never higher than min_matches)   *
 *----------------------------------------------------------------------*/
int count_kmers(char* read, dBGraph* db_graph, int* reads_with_n)
{
    int i;
    int count = 0;
    char kstr[kmer_size + 1];
    BinaryKmer b;
    BinaryKmer tmp_kmer;
    dBNode* node;
    
    for (i=0; i<=(strlen(read) - kmer_size); i++) {
        strncpy(kstr, read + i, kmer_size);
        kstr[kmer_size] = 0;
        if (strchr(kstr, 'N') == 0) {
            seq_to_binary_kmer(kstr, kmer_size, &b);
            node = hash_table_find(element_get_key(&b, kmer_size, &tmp_kmer), db_graph);
            if (node != NULL) {
                count++;
                if (count == min_matches) {
                    break;
                }
            }
        } else {
            *reads_with_n = *reads_with_n + 1;
        }
    }
    
    return count;
}

/*----------------------------------------------------------------------*
 * Function: filter_fastq_file                                          *
 * Purpose:  Read FASTQ file and write filtered version of it           *
 * Params:   filename -> name of input file                             *
 *           db_graph -> hash table to check against                    *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void filter_fastq_file(char* filename, dBGraph* db_graph)
{
    FILE* fp_in;
    FILE* fp_out;
    char output_filename[MAX_LINE_LENGTH];
    char read_header[MAX_LINE_LENGTH];
    char read[MAX_LINE_LENGTH];
    char quality_header[MAX_LINE_LENGTH];
    char quality_scores[MAX_LINE_LENGTH];
    int reads_with_n = 0;
    int reads_input = 0;
    int reads_output = 0;
        
    sprintf(output_filename, "%s.filtered.k%dm%d", filename, kmer_size, min_matches);

    printf("Filtering %s to %s...\n", filename, output_filename);    
    printf("Min matching kmers required: %d\n", min_matches);
	fflush(stdout);
    
    fp_in = fopen(filename, "r");
    if (!fp_in) {
        printf("Error: couldn't open input file %s\n", filename);
        return;
    }
    
    fp_out = fopen(output_filename, "w");
    if (!fp_out) {
        fclose(fp_in);
        printf("Error: couldn't open output file %s\n", output_filename);
        return;
    }
    
    while (!feof(fp_in)) {
        int success = 0;
        success += fgets(read_header, MAX_LINE_LENGTH, fp_in) ? 1:0;
        success += fgets(read, MAX_LINE_LENGTH, fp_in) ? 1:0;
        success += fgets(quality_header, MAX_LINE_LENGTH, fp_in) ? 1:0;
        success += fgets(quality_scores, MAX_LINE_LENGTH, fp_in) ? 1:0;
   
        if (success == 4) {
            reads_input++;
            
            if (read_header[0] != '@') {
                printf("Error: something wrong with read header: %s\n", read_header);
                break;
            }
            
            if (quality_header[0] != '+') {
                printf("Error: something wrong with quality header: %s\n", quality_header);
                break;
            }

            chomp(read_header);
            chomp(read);
            chomp(quality_header);
            chomp(quality_scores);
            if (count_kmers(read, db_graph, &reads_with_n) <  min_matches) {
                fprintf(fp_out, "%s\n", read_header);
                fprintf(fp_out, "%s\n", read);
                fprintf(fp_out, "%s\n", quality_header);
                fprintf(fp_out, "%s\n", quality_scores);
                reads_output++;
            }
        }
    }
    
    fclose(fp_out);
    fclose(fp_in);
    
    printf("Finished file... %d of %d written, with %d reads containing N\n", reads_output, reads_input, reads_with_n);
    fflush(stdout);
}

/*----------------------------------------------------------------------*
 * main                                                                 *
 *----------------------------------------------------------------------*/
int main (int argc, char * argv[])
{
    HashTable* db_graph = NULL;
    long long bad_reads;
    long int seq_length = 0;
    FILE* fp_fnames;
    char filename[MAX_LINE_LENGTH];

    parse_command_line_args(argc, argv);

    printf("Creating hash table...\n");
	fflush(stdout);
    db_graph = hash_table_new(hash_n, hash_b, 20, kmer_size);
    if (!db_graph) {
        log_and_screen_printf("Error: Couldn't get memory for graph\n");
        exit(1);
    }    
    
    printf("Loading reference %s... ", reference_filename);
	fflush(stdout);
    seq_length += load_fasta_from_filename_into_graph(reference_filename, 0, &bad_reads, 5000, db_graph);  
    printf("Read %ld nucleotides\n", seq_length);
 	fflush(stdout);

    // Open file of file names
    fp_fnames= fopen(file_of_filenames, "r");
    if (!fp_fnames) {
        printf("Can't open file of files\n");
        exit(1);
    }
    
	// For each file 
	while (!feof(fp_fnames)) {
		short colour = 0;
		fscanf(fp_fnames, "%s %hd\n", filename, &colour);	
		filter_fastq_file(filename, db_graph);
	}
	
	fclose(fp_fnames);

    printf("DONE\n");
    
    return 0;
}
