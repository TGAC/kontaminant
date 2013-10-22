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
#include "flags.h"
#include "binary_kmer.h"
#include "element.h"
#include "open_hash/hash_table.h"
#include "seq.h"
#include "element.h"
#include "dB_graph.h"
#include "file_reader.h"

/*----------------------------------------------------------------------*
 * Constants                                                            *
 *----------------------------------------------------------------------*/
#define NUMBER_OF_COLOURS 2
#define MAX_SEARCH_KMERS 100 

/*----------------------------------------------------------------------*
 * Global variables                                                     *
 *----------------------------------------------------------------------*/
HashTable* db_graph = NULL;
int hash_key_bits = 17;
int bucket_size = 60;
int kmer_size = 0;
int number_of_search_kmers = 0;
char* search_kmers[MAX_SEARCH_KMERS];
char* file_of_filenames;
char* output_filename;

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
    
	if (argc < 2)
	{
		printf("\nkmerinfo [-i filename] [-f kmer] [-k int] [-n int] [-b int]\n");
		printf("where [-k int] specifies kmer size\n");
		printf("      [-i filename] specifies the name of an input file of files (CTX only).\n");
		printf("      [-o filename] specifies the name of the output file where the frequencies of the kmer will be printed.\n");
		printf("      [-f kmer] specifies an individual kmer to provide information on (multiple allowed).\n");
		printf("      [-n int] specifies number of buckets in hash table (default 17).\n");
		printf("      [-b int] specifies size of hash table buckets (default 60).\n");
		exit(1);
    }
	
	while(i < argc)
	{
		char* parameter = argv[i];
		if (parameter[0] == '-')
		{
			switch (parameter[1])
			{
				case 'i':
					file_of_filenames = parse_string(argc, argv, &i);
					if (!file_of_filenames) {
						printf("ERROR: Invalid filename for -i parameter.\n");
						exit(1);
					}
					break;
				case 'f':
					search_kmers[number_of_search_kmers] = parse_string(argc, argv, &i);
					if (!search_kmers[number_of_search_kmers]) {
						printf("ERROR: Invalid filename for -f parameter.\n");
						exit(1);
					}
					number_of_search_kmers++;
					break;
				case 'k':
					kmer_size = parse_int(argc, argv, &i);
					break;
				case 'o':
					output_filename = parse_string(argc, argv, &i);
					if (!output_filename) {
						printf("ERROR: Invalid filename for -o parameter.\n");
						exit(1);
					}
					break;
				case 'n':
					hash_key_bits = parse_int(argc, argv, &i);
					break;
				case 'b':
					bucket_size = parse_int(argc, argv, &i);
					break;
				default:
					printf("ERROR: Invalid parameter\n");
					exit(1);
			}
		}
		
		i++;
	}
	
	if ((kmer_size < 1) || (kmer_size > 255)) {
		printf("ERROR: Invalid kmer size.\n");
		exit(1);
	}
}

/*----------------------------------------------------------------------*
 * Function: count_kmers                                                *
 * Purpose:  Count non-cleaned kmers.                                   *
 * Params:   None                                                       *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
long unsigned int count_kmers(long unsigned int *tc)
{
    long unsigned int number_of_kmers = 0;
    long unsigned int total_coverage = 0;
	
    void count_if_not_removed(dBNode * node) {
        if (db_node_check_for_flag(node, PRUNED) == false) {
            int total_edges = db_node_edges_count_all_colours(node, forward) + db_node_edges_count_all_colours(node, reverse);
            if (total_edges > 0) {
                number_of_kmers++;
                total_coverage += element_get_coverage_all_colours(node);
            }
        }
    }
	
    hash_table_traverse(&count_if_not_removed, db_graph);
	
    *tc = total_coverage;
    
    return number_of_kmers;
}

/*----------------------------------------------------------------------*
 * Function: kmer_info_print_coverage                                   *
 * Purpose:  Print coverage count file                                  *
 * Params:   out -> output file handle                                  *
 *           db_graph -> hash table                                     *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void kmer_info_print_coverage(FILE* out, dBGraph* db_graph)
{
    int max_cov = 0;
    int *count;
    int tmp_cov;
    
    fprintf(stdout, "Getting all the coverages...");
    
    // Function to find maximum coverage
    void f_max_cov(dBNode * node) {
        tmp_cov = element_get_coverage_all_colours(node);
        if (max_cov < tmp_cov) {
            max_cov = tmp_cov;
        }
    }
	
    // Function to find coverage counts
    void f_sum_cov(dBNode * node) {
        if (db_node_check_for_flag(node, PRUNED) == false) {
            tmp_cov = element_get_coverage_all_colours(node);
            count[tmp_cov]++;
        }
    }
    
    hash_table_traverse(&f_max_cov, db_graph);
    
    fprintf(stdout, "Trying to calloc %d ints\n", max_cov);
    fflush(stdout);
    
    // Calloc ensures all counts initially zero
    count = calloc(max_cov+1, sizeof(int));
	
    if (count) {
        hash_table_traverse(&f_sum_cov, db_graph);
        
        for (tmp_cov = 0; tmp_cov <=max_cov; tmp_cov++) {
            fprintf(out, "%d\t%d\n", tmp_cov, count[tmp_cov]);
        }
    } else {
        fprintf(stdout, "Error: Couldn't get memory.\n");
    }
    
    fprintf(stdout, "Done\n");
}

/*----------------------------------------------------------------------*
 * Function: main                                                       *
 * Purpose:  Entry point to program.                                    *
 * Params:   argc = number of arguments                                 *
 *           argv -> array of arguments                                 *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
int main (int argc, char * argv[])
{
    BinaryKmer bk;
    BinaryKmer tmp_kmer;
    int c;
    int i;
    long unsigned int n_kmers;
    FILE * out;
    FILE * fp_fnames;
    long unsigned int total_coverage;
    double average;
    char filename[1024];
	
    output_filename = NULL;//Make sure it is null, so we can use it in the if
	
    // Parse command line arguments
    parse_command_line_args(argc, argv);
    
    printf("\nkmerinfo\n\n");
	
    // Create hash table
    db_graph = hash_table_new(hash_key_bits, bucket_size, 20, kmer_size); 
    
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
		printf("Loading file %s colour %d...\n", filename, colour);
		load_binary_from_filename_into_graph(filename, db_graph, colour, 0);
	}
	
	fclose(fp_fnames);
    
	// Print stats
	hash_table_print_stats(db_graph);
    
	// Count non-cleaned kmers
	printf("Counting non-cleaned kmers... ");
	n_kmers = count_kmers(&total_coverage);
    average = (double)total_coverage / (double)n_kmers;
	printf("%ld kmers found ", n_kmers);
    printf("with total coverage %ld and average coverage %lf.\n", total_coverage, average);
	
	// Find kmers
	for (i=0; i<number_of_search_kmers; i++) {
		printf("Searching for kmer %s... ",search_kmers[i]);
		seq_to_binary_kmer(search_kmers[i], kmer_size, &bk);
		Element *e = hash_table_find(element_get_key(&bk, kmer_size, &tmp_kmer), db_graph);
		if (e == NULL) {
			printf("Couldn't find kmer.\n");
		} else {
			printf("Kmer found.\n");
			for (c=0; c<NUMBER_OF_COLOURS; c++) {
				printf("  Coverage colour %d: %d\n", c, e->coverage[c]);
			}
		}
	}
	
    // Output coverage file
	if (output_filename != NULL) {
		out = fopen(output_filename, "w");
        if (out) {
            fprintf(out, "Coverage\tCount\n");		    
            kmer_info_print_coverage(out, db_graph);
		    fclose(out);
        } else {
            printf("Error: Can't open output file.\n");
        }
	}
    
    return 0;
}
