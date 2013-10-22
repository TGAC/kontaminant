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
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <locale.h>
#include <getopt.h>
#include <err.h>
#include <unistd.h>
#include <fcntl.h>

#include <open_hash/hash_table.h>
#include <dB_graph.h>
#include <logger.h>
#include <file_format.h>
#include <file_reader.h>
#include <kmer_hash.h>
#include <kmer_reader.h>



typedef struct
{
    char * input_filename;
    char * reference_kmers;
    char * output_filename;
    char * output_histogram;
    FileFormat format;
    int kmer_size;
    int number_of_buckets_bits;
    int bucket_size;
    int quality_score_offset;
    int quality_score_threshold;
}KmerStatsCmdLine;

static void kmer_stats_default_opts(KmerStatsCmdLine * c)
{
    c->input_filename = calloc(LENGTH_FILENAME+1, sizeof(char *));
    c->output_filename = NULL;
    c->reference_kmers = NULL;
    c->output_histogram = NULL;
    c->kmer_size = 21;
    c->number_of_buckets_bits = 16;
    c->bucket_size = 100;
    c->quality_score_offset = 33;
    c->input_filename[0]='-';
    c->format = FASTQ;
    c->quality_score_threshold = 0;
}
const char *usage =
"\nusage: kmer_contamination [-h] [--input file_of_sequences] [--reference kmers_file][--mem_height 16][--output stats_file] \n" \
"\t --mem_width = Size of hash table buckets (default 100).\n" \
"\t --help = This help message"\
"\t --input = File of reads in fastq format.\n"\
"\t --file_format = Format of the read. fasta or fastq. Default fastq\n"\
"\t --kmer_size = Size of the kmer in the reference file\n"\
"\t --mem_heigh = Number of buckets in hash table in bits (default 10, this is a power of 2, ie 2^mem_height).\n"
"\t --output = A tab separated file with the results of the analysis. The header is printed if the file doesn't exist. If it exists, the new stats are appeneded. We don't validate that previous entries are valid.\n"\
"\t --quality_score_offset = Offset for the fastq file. Default 33.\n"\
"\t --quality_score_threshold = Minimum quality across a kmer to be used\n"\
"\t --reference_kmers = File containing the kmers to compare. This are built with kmer_hash_build\n"\
"\t --histogram_output = Output file where the histogram of encounter kmers per read is to be written\n";

static void print_help(){
	fprintf(stderr, "%s \n", usage);
    exit(-1);
}

static KmerStatsCmdLine parse_cmdline(int argc, char *argv[], int unit_size)
{
    int i;
    printf("Command: ");
    for (i = 0; i < argc; i++) {
        printf("%s ", argv[i]);
    }
    printf("\n");
    printf("Unit size: %i\n", unit_size);

    KmerStatsCmdLine cmd_line;
    kmer_stats_default_opts(&cmd_line);
    int opt;
    int longopt_index;
     char * tmp_format;
     static struct option long_options[] = {
         {"mem_width", required_argument, NULL, 'b'},
         {"file_format", required_argument, NULL, 'f'},
         {"help", no_argument, NULL, 'h'},
         {"input", required_argument, NULL, 'i'},
         {"kmer_size", required_argument, NULL, 'k'},
         {"mem_height", required_argument, NULL, 'n'},
         {"output", required_argument, NULL, 'o'},
         {"quality_score_offset", required_argument, NULL, 'p'},
         {"quality_score_threshold", required_argument, NULL, 'q'},
         {"reference_kmers", required_argument, NULL, 'r'},
         {"histogram_output", required_argument, NULL, 't'}
     };
     while ((opt = getopt_long(argc, argv,"b:f:hi:k:n:o:p:q:r:t:", long_options, &longopt_index)) > 0){
         switch (opt) {
             case 'b':
                 if (optarg == NULL)
                     errx(1, "[-b | --mem_width] option requires int argument [hash table bucket size]");
                 cmd_line.bucket_size = atoi(optarg);
                 if (cmd_line.bucket_size == 0 || cmd_line.bucket_size > SHRT_MAX)	//check that -b is not bigger than max_short 
                     errx(1, "[-b | --mem_width] option requires 'short' argument bigger than 0");
                 break;
                 
             case 'f':
                 
                 tmp_format = calloc(20, sizeof(char));
                 strcpy(tmp_format, optarg);
                 cmd_line.format = string_to_file_format(tmp_format);
                 free(tmp_format);
                 break;
             case 'h':
				print_help();
				break;
             case 'i': //file of filenames
                 if (optarg == NULL)
                     errx(1,"[-i | --input] option requires a filename [kmers file]");
                 
                 if (strlen(optarg) < LENGTH_FILENAME) {
                     if(cmd_line.input_filename == NULL){
                         cmd_line.input_filename = calloc(LENGTH_FILENAME+1, sizeof(char *));
                     }
                     strcpy(cmd_line.input_filename, optarg);
                 } else {
                     errx(1, "[-i | --input] filename too long [%s]", optarg);
                 }
                 
                 if (access(optarg, R_OK) == -1) {errx(1,"[-i | --input] filename [%s] cannot be accessed", optarg);
                 }
                 break;
             case 'k':	//kmer size
                 if (optarg == NULL){
                     errx(1,"[-k | --kmer_size] option requires int argument [kmer size]");
                 }
                 cmd_line.kmer_size = atoi(optarg);
                 
                 if (cmd_line.kmer_size == 0){
                     errx(1,"[-k | --kmer_size] option requires int argument bigger than 0");
                 }
                 if(cmd_line.kmer_size % 2 == 0){
                     errx(1, "[-k | --kmer_size] option can't be even");
                 }
                 if (cmd_line.kmer_size > (NUMBER_OF_BITFIELDS_IN_BINARY_KMER * 32)) {
                     errx(1, "[-k | --kmer_size] The kmer size is greater than the maximum kmer_size. ");
                 }
                 break;
             case 'n':	//number of buckets
                 if (optarg == NULL)
                     errx(1,"[-n | --mem_height] option requires int argument [hash table number of buckets in bits]");
                 cmd_line.number_of_buckets_bits = atoi(optarg);
                 break;
             case 'p':	//quality_score_offset
                 if (optarg == NULL){
                     errx(1, "[-p | --quality_score_offset] option requires int argument");
                 }
                 cmd_line.quality_score_offset = atoi(optarg);
                 
                 if (cmd_line.quality_score_offset == 0){
                     errx(1,"[-p | --quality_score_offset] option requires int argument bigger than 0");
                 }
                 break;
                 
             case 'q':	//quality threshold
                 if (optarg == NULL){
                     errx(1, "[-q | --quality_score_threshold] option requires int argument [quality score threshold]");
                 }
                 cmd_line.quality_score_threshold = atoi(optarg);
                 
                 if (cmd_line.quality_score_threshold == 0){
                     errx(1, "[-q | --quality_score_threshold] option requires int argument bigger than 0");
                 }
                 break;  
              
              case 'o': //file of filenames
                 if (optarg == NULL)
                     errx(1,"[-o | --output] option requires a filename [output file]");

                 if (strlen(optarg) < LENGTH_FILENAME) {
                     if(cmd_line.output_filename == NULL){
                         cmd_line.output_filename = calloc(LENGTH_FILENAME+1, sizeof(char *));
                     }
                     strcpy(cmd_line.output_filename, optarg);
                 } else {
                     errx(1, "[-o | --output] filename too long [%s]", optarg);
                 }

                 
                  if (access(optarg, F_OK) == 0) {
                    if (access(optarg, W_OK) == -1) {
						errx(1,"[-o | --output] filename [%s] cannot be modified", optarg);
					}
                 }
                  
                  break; 
              case 'r': //file of filenames
                 if (optarg == NULL)
                     errx(1,"[-r | --reference_kmers] option requires a filename [kmers file]");

                 if (strlen(optarg) < LENGTH_FILENAME) {
                     if(cmd_line.reference_kmers == NULL){
                         cmd_line.reference_kmers = calloc(LENGTH_FILENAME+1, sizeof(char *));
                     }
                     strcpy(cmd_line.reference_kmers, optarg);
                 } else {
                     errx(1, "[-r | --reference_kmers] filename too long [%s]", optarg);
                 }

                  if (access(optarg, R_OK) == -1) {
                     errx(1,"[-r | --reference_kmers] filename [%s] cannot be accessed", optarg);
                  }
                  break;
             case 't': //file of filenames
                 if (optarg == NULL)
                     errx(1,"[-t | --histogram_output] option requires a filename [output file]");
                 
                 if (strlen(optarg) < LENGTH_FILENAME) {
                     if(cmd_line.output_histogram == NULL){
                         cmd_line.output_histogram = calloc(LENGTH_FILENAME+1, sizeof(char *));
                     }
                     strcpy(cmd_line.output_histogram, optarg);
                 } else {
                     errx(1, "[-t | --histogram_output] filename too long [%s]", optarg);
                 }
                 
                 
                 if (access(optarg, F_OK) == 0) {
                     if (access(optarg, W_OK) == -1) {
                         errx(1,"[-t | --histogram_output] filename [%s] cannot be modified", optarg);
                     }
                 }
                 
                 break; 
             default:
                 fprintf(stderr, "Unknown argument: %c \n", opt);
         }
         
     }
        
    if (cmd_line.input_filename == NULL) {
        errx(1,"[-i | --input] option requires a sample file [fastq file]");
    }
    
    if (cmd_line.reference_kmers == NULL) {
        errx(1,"[-r | --reference_kmers] option requires a filename [kmers file]");
    }
    return  cmd_line;
    
}

static KmerHash * load_kmer_table(KmerStatsCmdLine cmd_line){
    //TODO: Use a special format that contains the memory requirements. 
    
     
    log_and_screen_printf("\nHash table from file: %s\n", cmd_line.reference_kmers);
    KmerHash * kmer_hash = hash_table_new(cmd_line.number_of_buckets_bits,
                              cmd_line.bucket_size, 25, cmd_line.kmer_size);
    boolean all_entries_are_unique = true;
    
    log_and_screen_printf("\nReading kmers file: %s\n", cmd_line.reference_kmers);		
    fflush(stdout);
    load_binary_from_filename_into_kmers_hash(cmd_line.reference_kmers, kmer_hash, KMER_HASH_REFERENCE_INDEX, all_entries_are_unique);
    
    log_and_screen_printf("\nRead of file complete. Total kmers: %'lld\n", hash_table_get_unique_kmers(kmer_hash));
    hash_table_print_stats(kmer_hash);
    return kmer_hash;

}   
static long long load_reads_coverage_table(KmerStatsCmdLine cmd_line,  KmerHash * kmer_hash){
    
    
    log_and_screen_printf("\nLoading sample from: %s\n", cmd_line.input_filename);
    long long loaded_kmers = 0;
    if(cmd_line.format == KMERS){
        //boolean all_entries_are_unique = false;
        //loaded_kmers = load_binary_from_filename_into_kmers_hash(cmd_line.input_filename, kmer_hash, KMER_HASH_SAMPLE_INDEX, all_entries_are_unique);
        loaded_kmers = load_kmers_binary_from_filename_update_coverage(cmd_line.input_filename, kmer_hash, KMER_HASH_SAMPLE_INDEX);
    }else{
        KmerFileReaderArgs fra;
        fra.bad_reads = 0;
        fra.colour = KMER_HASH_SAMPLE_INDEX;
        fra.fastq_ascii_offset = cmd_line.quality_score_offset;
        KmerFileReaderInnerArgs fria;
        fria.format = cmd_line.format;
        fria.kmer_size = kmer_hash->kmer_size;
        fria.max_read_length = 1000;
        fria.new_entry = true;
        fra.filename = cmd_line.input_filename;
        fra.quality_cut_off = cmd_line.quality_score_threshold;
        fra.inner_args = &fria;
        fra.insert = false;
        fra.max_read_length = 1000;
        fra.maximum_ocupancy = 75;
        fra.KmerHash = kmer_hash;
           //fp, seq, max_read_length, full_entry, &full_entry
        
        loaded_kmers = load_seq_into_kmers_hash(&fra);
        log_and_screen_printf("Loaded %'lld kmers (bad reads %'lld)", loaded_kmers, fra.bad_reads);
        hash_table_print_stats(kmer_hash);
    }
   
    
    
    return loaded_kmers;
}

static void print_kmer_stats(KmerStatsCmdLine * cmd_line, KmerHash * kmer_hash){
    
        kmer_hash_print_kmer_stats(cmd_line->output_filename,cmd_line->reference_kmers, cmd_line->input_filename,kmer_hash);
    
    
}

static void print_contaminated_kmers_histogram(KmerStatsCmdLine * cmd_line, KmerHash * kmer_hash){
    kmer_hash_print_contaminated_kmers_histogram(cmd_line->output_histogram,  kmer_hash);
    
}


int main(int argc, char **argv)
{
	setlocale (LC_ALL, "");

    log_and_screen_printf("\nkmer_contamination.\n\n");    
	log_and_screen_printf(SVN_VERSION);
	log_and_screen_printf(SVN_COMMIT_DATE);
	log_and_screen_printf("Compiled on %s at %s \n\n", __DATE__, __TIME__);
    
    
    KmerStatsCmdLine cmd_line = parse_cmdline(argc, argv, sizeof(Element));
    
    
    KmerHash * kmer_hash = load_kmer_table(cmd_line);
    load_reads_coverage_table(cmd_line, kmer_hash);
    print_kmer_stats(&cmd_line, kmer_hash);
    print_contaminated_kmers_histogram(&cmd_line, kmer_hash);
    log_and_screen_printf("\nDONE");
    return 0;
    
}
