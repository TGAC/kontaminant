//
//  kmer_hash_build.c
//  Untitled
//
//  Created by Ricardo Ramirez-Gonzalez (TGAC) on 14/12/2011.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <stdio.h>

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <locale.h>
#include <getopt.h>
#include <err.h>
#include <unistd.h>

#include <open_hash/hash_table.h>
#include <element.h>
#include <dB_graph.h>
#include <logger.h>

#include <file_format.h>
#include <nucleotide.h>
#include <seq.h>
#include <kmer_hash.h>
#include <kmer_reader.h>


typedef struct
{
    char * input_filename;
    char * output_filename;
    char ** input_file_list;
    short kmer_size;
    int number_of_buckets_bits;
    int bucket_size;
    int quality_score_offset;
    int quality_score_threshold;
    int input_file_list_size;
    FileFormat input_format;
}KmerHashBuildCmdLine;


static  const char * usage = 
"\nusage: kmer_hash_build [-h] [--input file_of_sequences][--mem_height 16][--output kmers_file] \n" \
"\t --mem_width = Size of hash table buckets (default 100).\n" \
"\t --help = This help message"\
"\t --input = File of reads in fastq or fasta format. If not privided, STDIN is used\n"\
"\t --file_format = Format of the reads. fasta or fastq. Default fastq\n"\
"\t --kmer_size = Size of the kmer in the reference file\n"\
"\t --mem_heigh = Number of buckets in hash table in bits (default 10, this is a power of 2, ie 2^mem_height).\n"
"\t --output = A kmers file in binary format. To be used with the kmer_stats and kmer_contamination programs.\n"\
"\t --quality_score_offset = Offset for the fastq file. Default 33.\n"\
"\t --quality_score_threshold = Minimum quality across a kmer to be used\n";


static void print_help(){
	fprintf(stderr, "%s \n", usage);
    exit(-1);
}
static void timestamp() {
	time_t ltime = time(NULL);
	log_and_screen_printf("\n-----\n%s",asctime(localtime(&ltime)));	
	fflush(stdout);
}

static void kmer_hash_build_default_opts(KmerHashBuildCmdLine * c)
{
    c->input_filename = "-";
    c->output_filename = NULL;
    c->kmer_size = 21;
    c->number_of_buckets_bits = 16;
    c->bucket_size = 100;
    c->quality_score_offset = 33;
    c->quality_score_threshold = 0;
    c->input_format = UNSPECIFIED_FORMAT;
    c->input_file_list = NULL;
    c->input_file_list_size = 0;
}


static KmerHashBuildCmdLine parse_buid_hash_cmdline(int argc, char *argv[], int unit_size)
{
    int i;
    printf("Command: ");
    for (i = 0; i < argc; i++) {
        printf("%s ", argv[i]);
    }
    printf("\n");
    printf("Unit size: %i\n", unit_size);
    
    KmerHashBuildCmdLine cmd_line;
    kmer_hash_build_default_opts(&cmd_line);
    int opt;
    int longopt_index;
    char * tmp_format;
   
    static struct option long_options[] = {
        {"mem_width", required_argument, NULL, 'b'},
        {"format", required_argument, NULL, 'f'},
        {"help", no_argument, NULL, 'h'},
        {"input", required_argument, NULL, 'i'},
        {"kmer_size", required_argument, NULL, 'k'},
        {"mem_height", required_argument, NULL, 'n'},
        {"output", required_argument, NULL, 'o'}, 
        {"quality_score_offset", required_argument, NULL, 'p'},
        {"quality_score_threshold", required_argument, NULL, 'q'},
       
    };
    while ((opt = getopt_long(argc, argv,"b:f:i:k:n:o:p:q:", long_options, &longopt_index)) > 0){
        
        
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
                cmd_line.input_format = string_to_file_format(tmp_format);
                
                free(tmp_format);
                break;
            case 'h':
				print_help();
				break;
            case 'i': //file of filenames
                if (optarg == NULL)
                    errx(1,"[-i | --input] option requires a filename [kmers file]");
                cmd_line.input_filename = optarg;
                
                if (access(optarg, R_OK) == -1) {errx(1,"[-i | --input] filename [%s] cannot be accessed", optarg);
                }
                break;
            case 'k':	//kmer size
                if (optarg == NULL){
                    errx(1,"[-k | --kmer_size] option requires int argument [kmer size]");
                }
                cmd_line.kmer_size = (short) atoi(optarg);
                
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
            case 'o': //file of filenames
                if (optarg == NULL)
                    errx(1,"[-o | --output] option requires a filename [kmers file]");
                
                if (strlen(optarg) < LENGTH_FILENAME) {
                    if(cmd_line.output_filename == NULL){
                        cmd_line.output_filename = calloc(LENGTH_FILENAME+1, sizeof(char *));
                    }
                    strcpy(cmd_line.output_filename, optarg);
                } else {
                    errx(1, "[-o | --output] filename too long [%s]", optarg);
                }
                
                if (access(optarg, F_OK) == 0) {
                    errx(1, "[-o | --output] filename [%s] already existes", optarg);
                }
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
            
            default:
                fprintf(stderr, "Unknown argument: %c \n", opt);
        }
        
    }
    
    if (optind < argc) {
        cmd_line.input_file_list_size = argc - optind;
        cmd_line.input_file_list = &argv[optind];
        printf ("non-option ARGV-elements: ");
        while (optind < argc)
            printf ("%s ", argv[optind++]);
        printf ("\n");
    }
    
    
    if (cmd_line.input_filename == NULL  &&  cmd_line.input_file_list_size == 0) {
        errx(1,"[-i | --input] option requires a sample file [fastq file]");
    }
    
    
    if (cmd_line.output_filename == NULL) {
        errx(1,"[-o | --output_filename] option requires a file [kmer hash file]");
    }
    
    if (! (cmd_line.input_format == FASTQ || cmd_line.input_format == FASTA  || cmd_line.input_format == KMERS) ) {
        errx(1,"[-f | --format] option requires a suported format (FASTQ | FASTA | KMERS)");
    }
    
    return  cmd_line;
    
}

//At the moment, we are only using the traditional n and m. We will look to set a maximum amount of memory at some point
static KmerHash * alloc_kmer_hash(KmerHashBuildCmdLine * cmd_line){
    KmerHash * kmer_hash;
    
    int n = cmd_line->number_of_buckets_bits;
    int b =  cmd_line->bucket_size;
    //Here we can add calculations. 
    
    kmer_hash = hash_table_new(n,b, 25, cmd_line->kmer_size);
    
    return kmer_hash;
} 


static void load_reads_into_table(KmerHashBuildCmdLine * cmd_line,  KmerHash * kmer_hash){
    log_and_screen_printf("\nHash table from file: %s\n", cmd_line->input_filename);
    
    KmerFileReaderArgs fra;
    fra.bad_reads = 0;
    fra.colour = 0;
    fra.fastq_ascii_offset = cmd_line->quality_score_offset;
    KmerFileReaderInnerArgs fria;
    fria.format = cmd_line->input_format;
    fria.kmer_size = kmer_hash->kmer_size;
    fria.max_read_length = 1000;
    fria.new_entry = true;
    fria.full_entry = false;
    fra.filename = cmd_line -> input_filename;
    fra.quality_cut_off = cmd_line->quality_score_threshold;
    fra.inner_args = &fria;
    fra.insert = true;
    fra.max_read_length = 1000;
    fra.maximum_ocupancy = 75;
    fra.KmerHash = kmer_hash;
     //fp, seq, max_read_length, full_entry, &full_entry
    long long loaded_kmers = 0;
    loaded_kmers = load_seq_into_kmers_hash(&fra);
    
    printf("Loaded %'lld kmers (bad reads %'lld)", loaded_kmers, fra.bad_reads);
    
}

static void output_basic_info(KmerHashBuildCmdLine  cmd_line)
{
    log_and_screen_printf("Max k: %i\n", (NUMBER_OF_BITFIELDS_IN_BINARY_KMER*32)-1);
    if (cmd_line.input_format == FASTQ) { 
        log_and_screen_printf("Quality score offset: %i", cmd_line.quality_score_offset);
        if (cmd_line.quality_score_offset == 33) {
            log_and_screen_printf(" (Sanger format)");
        } else if (cmd_line.quality_score_offset == 64) {
            log_and_screen_printf(" (Solexa/Illumina format)");
        } else {
            log_and_screen_printf(" (Unknown format)");
        }
        log_and_screen_printf("\nQuality score threshold: %i\n", cmd_line.quality_score_threshold);
    }
}

static void dump_kmer_hash(KmerHashBuildCmdLine * cmd_line, KmerHash * kmer_hash){
    timestamp();
    log_and_screen_printf("\nDumping hash table to file: %s\n", cmd_line->output_filename);
    kmer_hash_dump_binary(cmd_line->output_filename, &db_node_check_flag_not_pruned, kmer_hash);
    fflush(stdout);
}

static void load_kmers_from_binary(KmerHashBuildCmdLine * cmd_line, KmerHash * kmer_hash){
    boolean all_entries_are_unique = true;
    int i;
    for(i = 0; i < cmd_line->input_file_list_size; i++){
        char * filename = cmd_line->input_file_list[i];
        log_and_screen_printf("\nReading kmers file: %s\n", filename);		
        fflush(stdout);
        load_binary_from_filename_into_kmers_hash(filename, kmer_hash, 0, all_entries_are_unique);
        all_entries_are_unique = false;
        
        log_and_screen_printf("\nRead of file complete. Total kmers: %'lld\n", hash_table_get_unique_kmers(kmer_hash));
        hash_table_print_stats(kmer_hash);
    }
}

int main(int argc, char **argv)
{
	setlocale (LC_ALL, "");
    
    log_and_screen_printf("\nkmer hash build.\n\n");    
	log_and_screen_printf(SVN_VERSION);
	log_and_screen_printf(SVN_COMMIT_DATE);
	log_and_screen_printf("Compiled on %s at %s \n\n", __DATE__, __TIME__);
    
    KmerHashBuildCmdLine cmd_line = parse_buid_hash_cmdline(argc, argv, sizeof(Element));
    output_basic_info(cmd_line);
    KmerHash * kmer_hash = alloc_kmer_hash(&cmd_line);
    
    if (cmd_line.input_format == KMERS) {
        load_kmers_from_binary(&cmd_line, kmer_hash);
    }else{
        load_reads_into_table(&cmd_line, kmer_hash);
    }
    
   /* KmerHash * kmer_hash = load_kmer_table(cmd_line);
    load_reads_coverage_table(cmd_line, kmer_hash);
    print_kmer_stats(kmer_hash);*/
    //print_kmer_stats(kmer_hash);
    hash_table_print_stats(kmer_hash);
    dump_kmer_hash(&cmd_line,kmer_hash);
    log_and_screen_printf("\nDONE\n");
    return 0;
    
}
