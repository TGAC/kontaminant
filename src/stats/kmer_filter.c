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
#include <assert.h>

#include <open_hash/hash_table.h>
#include <dB_graph.h>
#include <logger.h>
#include <file_format.h>
#include <file_reader.h>
#include <kmer_hash.h>
#include <kmer_reader.h>
#include <seq_io.h>

typedef struct
{
    char * input_read_1;
    char * input_read_2;
    char * reference_kmers;
    char * suffix_in;
    char * output_prefix;
    char * stats_file;
    char * output_folder;
    FileFormat format;
    int kmer_size;
    int number_of_buckets_bits;
    int bucket_size;
    int quality_score_offset;
    int quality_score_threshold;
    int max_read_length;
    int max_name_length;
    int threashold;
    int print_contaminated;
    
}KmerFilterCmdLine;


static void kmer_stats_default_opts(KmerFilterCmdLine * c)
{
    c->input_read_1 = NULL;
    c->input_read_2 = NULL;
    c->suffix_in = NULL;
    c->output_prefix = NULL;
    c->reference_kmers = NULL;
    c->stats_file = NULL;
    c->output_folder = NULL;
    c->kmer_size = 21;
    c->number_of_buckets_bits = 16;
    c->bucket_size = 100;
    c->quality_score_offset = 33;
    c->format = FASTQ;
    c->quality_score_threshold = 0;
    c->max_name_length = 2000;
    c->max_read_length = 2000;
    c->threashold = 0;
    c->print_contaminated = 0;
}

const char *usage =
"\nusage: kmer_filter [-h] [--input file_of_sequences] [--reference kmers_file][--mem_height 16][--output stats_file] \n" \
"\t --mem_width = Size of hash table buckets (default 100).\n" \
"\t --help = This help message\n"\
"\t --read_1 = File of reads in fastq format.\n"\
"\t --read_2 = File of reads in fastq format. Optional. Should be in the same reads in the same order than the reads in read_1\n"\
"\t --file_format = Format of the read. fasta or fastq. Default fastq\n"\
"\t --kmer_size = Size of the kmer in the reference file\n"\
"\t --mem_heigh = Number of buckets in hash table in bits (default 10, this is a power of 2, ie 2^mem_height).\n"
"\t --output_prefix = The prefix for the files with bins.\n"\
"\t --quality_score_offset = Offset for the fastq file. Default 33.\n"\
"\t --quality_score_threshold = Minimum quality across a kmer to be used\n"\
"\t --reference_kmers = File containing the kmers to compare. This are built with kmer_hash_build\n"\
"\t --print_contaminated = Print contaminated reads in a new fastq file\n"\
"\t --threashold = Number of kmers to be considered as a above or below the cutoff" \
"\t --output_folder = Folder to write the filtered reads. "
;
static void print_help(){
	fprintf(stderr, "%s \n", usage);
    exit(-1);
}

static KmerFilterCmdLine parse_cmdline(int argc, char *argv[], int unit_size)
{
    int i;
    printf("Command: ");
    for (i = 0; i < argc; i++) {
        printf("%s ", argv[i]);
    }
    printf("\n");
    printf("Unit size: %i\n", unit_size);
    
    KmerFilterCmdLine cmd_line;
    kmer_stats_default_opts(&cmd_line);
    int opt;
    int longopt_index;
    char * tmp_format;
    static struct option long_options[] = {
        {"read_1", required_argument, NULL, '1'},
        {"read_2",required_argument, NULL, '2'},
        {"mem_width", required_argument, NULL, 'b'},
        {"print_contaminated", no_argument, NULL, 'c'},
        {"file_format", required_argument, NULL, 'f'},
        {"help", no_argument, NULL, 'h'},
        {"kmer_size", required_argument, NULL, 'k'},
        {"mem_height", required_argument, NULL, 'n'},
        {"suffix_in", required_argument, NULL, 'i'},
        {"output_prefix", required_argument, NULL, 'o'},
        {"quality_score_offset", required_argument, NULL, 'p'},
        {"quality_score_threshold", required_argument, NULL, 'q'},
        {"reference_kmers", required_argument, NULL, 'r'},
        {"stats_file", required_argument, NULL, 's'},
        {"threashold", required_argument, NULL, 't'},
        {"output_folder", required_argument, NULL, 'O'},
    };
    while ((opt = getopt_long(argc, argv,"1:2:b:cf:hi:k:n:o:p:q:r:s:t:O:", long_options, &longopt_index)) > 0){
        switch (opt) {
            case '1': //file of filenames
                if (optarg == NULL)
                    errx(1,"[-1 | --read_1] option requires a filename [read 1]");
                
                if (strlen(optarg) < LENGTH_FILENAME) {
                    if(cmd_line.input_read_1 == NULL){
                        cmd_line.input_read_1 = calloc(LENGTH_FILENAME+1, sizeof(char *));
                    }
                    strcpy(cmd_line.input_read_1, optarg);
                } else {
                    errx(1, "[-1 | --read_1] filename too long [%s]", optarg);
                }
                
                if (access(optarg, R_OK) == -1) {
                    errx(1,"[-1 | --read_1] filename [%s] cannot be accessed", optarg);
                }
                break;
            case '2': //file of filenames
                if (optarg == NULL)
                    errx(1,"[-2 | --read_2] option requires a filename [read 2]");
                
                if (strlen(optarg) < LENGTH_FILENAME) {
                    if(cmd_line.input_read_2 == NULL){
                        cmd_line.input_read_2 = calloc(LENGTH_FILENAME+1, sizeof(char *));
                    }
                    strcpy(cmd_line.input_read_2, optarg);
                } else {
                    errx(1, "[-2 | --read_2] filename too long [%s]", optarg);
                }
                
                if (access(optarg, R_OK) == -1) {errx(1,"[-2 | --read_2] filename [%s] cannot be accessed", optarg);
                }
                break;
            case 'b':
                if (optarg == NULL)
                    errx(1, "[-b | --mem_width] option requires int argument [hash table bucket size]");
                cmd_line.bucket_size = atoi(optarg);
                if (cmd_line.bucket_size == 0 || cmd_line.bucket_size > SHRT_MAX)	//check that -b is not bigger than max_short 
                    errx(1, "[-b | --mem_width] option requires 'short' argument bigger than 0");
                break;
            case 'c':
                cmd_line.print_contaminated = 1;
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
            case 'i': //Sufix to prepend
                if (optarg == NULL)
                    errx(1,"[-i | --suffix_in] option requires a suffix [suffix]");
                
                if (strlen(optarg) < LENGTH_FILENAME) {
                    if(cmd_line.suffix_in == NULL){
                        cmd_line.suffix_in = calloc(LENGTH_FILENAME+1, sizeof(char *));
                    }
                    strcpy(cmd_line.suffix_in, optarg);
                } else {
                    errx(1, "[-i | --suffix_in] suffix too long [%s]", optarg);
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
                
            case 'o': //Sufix to prepend
                if (optarg == NULL)
                    errx(1,"[-o | --output_prefix] option requires a suffix [suffix]");
                
                if (strlen(optarg) < LENGTH_FILENAME) {
                    if(cmd_line.output_prefix == NULL){
                        cmd_line.output_prefix = calloc(LENGTH_FILENAME+1, sizeof(char *));
                    }
                    strcpy(cmd_line.output_prefix, optarg);
                } else {
                    errx(1, "[-o | --output_prefix] suffix too long [%s]", optarg);
                }
                
                break; 
            case 'r': 
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
            case 's': 
                if (optarg == NULL)
                    errx(1,"[-s | --stats_file] option requires a filename ");
                
                if (strlen(optarg) < LENGTH_FILENAME) {
                    if(cmd_line.stats_file == NULL){
                        cmd_line.stats_file = calloc(LENGTH_FILENAME+1, sizeof(char *));
                    }
                    strcpy(cmd_line.stats_file, optarg);
                } else {
                    errx(1, "[-s | --stats_file] filename too long [%s]", optarg);
                }
            case 't':
                if (optarg == NULL) {
                    errx(1, "[-t | --threashold] option requires a number float between 0 and 1");
                }
                cmd_line.threashold = atof(optarg);
                if ( cmd_line.threashold  < 0 ||  cmd_line.threashold  > MAX_READ_LENGTH) {
                    errx(1, "[-t | --threashold] option requires a number float between 0 and %d", MAX_READ_LENGTH);
                }
                
            case 'O': //Sufix to prepend
                if (optarg == NULL)
                    errx(1,"[-O | --output_folder] option requires a suffix [suffix]");
                
                if (strlen(optarg) < LENGTH_FILENAME) {
                    if(cmd_line.output_folder == NULL){
                        cmd_line.output_folder = calloc(LENGTH_FILENAME+1, sizeof(char *));
                    }
                    strcpy(cmd_line.output_folder, optarg);
                } else {
                    errx(1, "[-O | --output_folder] suffix too long [%s]", optarg);
                }   
                
                break; 
            default:
                fprintf(stderr, "Unknown argument: %c \n", opt);
        }
        
    }
    
    if (cmd_line.input_read_1 == NULL) {
        errx(1,"[-1 | --read_1] option requires a sample file [fastq file]");
    }
    
    if (cmd_line.reference_kmers == NULL) {
        errx(1,"[-r | --reference_kmers] option requires a filename [kmers file]");
    }
    return  cmd_line;
    
}



static KmerHash * load_kmer_table(KmerFilterCmdLine cmd_line){
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



typedef struct {
    char ** filenames;
    int capacity;
    FileFormat format;
    FILE ** file;
} kmer_bin_filenames;


static kmer_bin_filenames * init_filenames(char * prefix, short read ,int max_kmers, FileFormat format){
    kmer_bin_filenames * bin_files = calloc(1, sizeof(kmer_bin_filenames));
   
    int    bins = max_kmers + 3 ;
    bin_files->capacity = bins;
    bin_files->filenames = calloc(bins, sizeof(char *));
    bin_files->file = calloc(bins, sizeof(FILE *));
    bin_files->format = format;
    int i ;
    
    for (i = 0; i < max_kmers; i++) {
        bin_files->filenames[i] = calloc(strlen(prefix)+20, sizeof(char));
        sprintf(bin_files->filenames[i], "%s_R%d_KMERS%d.%s", prefix, read, i, file_format_to_string(format) );
    }   

    bin_files->filenames[i] = calloc(strlen(prefix)+20, sizeof(char));
    sprintf(bin_files->filenames[i], "%s_R%d_clean.%s", prefix, read,  file_format_to_string(format) );
    i++;
    bin_files->filenames[i] = calloc(strlen(prefix)+20, sizeof(char));
    sprintf(bin_files->filenames[i], "%s_R%d_contaminated.%s", prefix, read,  file_format_to_string(format) );
    
    return bin_files;
}

static void free_filenames(kmer_bin_filenames * bin_files){
    int i;
    for (i = 0; i < bin_files->capacity; i++) {
        
        free(bin_files->filenames[i]);
    }
    free( bin_files->filenames);
    free(bin_files);
}

/*
 * Wrapper function that avoids opening and closing the file. It only opens the file once. 
 */
static void print_sequnce_in_kmer_bin(Sequence * seq, int kmers, kmer_bin_filenames * bin_files){
    assert(kmers <= bin_files->capacity);
  //  int i = kmers;
    if ( bin_files->file[kmers] == NULL) {
        
        bin_files->file[kmers] = fopen(bin_files->filenames[kmers], "w");
        if(bin_files->file[kmers] == NULL){
            log_and_screen_printf("Unable to open file: %s\n", bin_files->filenames[kmers] );
            exit(-1);
        }
    }
    
    append_sequence_fh(bin_files->file[kmers], seq, bin_files->format);
    
}


static void filter_reads(KmerFilterCmdLine * cmd_line, KmerHash * kmer_hash){
    
    int length_1 = 0;
    int length_2 = 0;
    
    
    KmerFileReaderInnerArgs * read_1 = open_kmer_file_reader_wrapper(cmd_line->input_read_1, cmd_line->kmer_size, cmd_line->max_read_length, cmd_line->format, cmd_line->quality_score_offset);
    
    KmerFileReaderInnerArgs * read_2 = open_kmer_file_reader_wrapper(cmd_line->input_read_2, cmd_line->kmer_size, cmd_line->max_read_length, cmd_line->format, cmd_line->quality_score_offset);
    
    KmerSlidingWindowSet * kmers_1 = binary_kmer_sliding_window_set_new_from_read_length(cmd_line->kmer_size, cmd_line->max_read_length);
    

    
    
    
    kmer_bin_filenames * files_1 = NULL;
    kmer_bin_filenames * files_2 = NULL;
    
    if (cmd_line->output_prefix != NULL) {
        files_1 = init_filenames(cmd_line->output_prefix, 1, kmers_1->max_kmers, cmd_line->format);
        if(read_2) {
            files_2 = init_filenames(cmd_line->output_prefix, 2, kmers_1->max_kmers, cmd_line->format);
        }
    }
    
    
    int total_length = 0;
    int total_kmers_in_hash = 0;
    int current_kmers_in_path = 0;
    
    while((length_1 = file_reader_wrapper(read_1)) ){
        if (read_2) {
            length_2 = file_reader_wrapper(read_2);
        }
        
        total_length += get_sliding_windows(read_1->seq, cmd_line->quality_score_threshold, kmers_1);
        current_kmers_in_path = kmer_hash_present_kmers_in_sliding_window_set(kmers_1, kmer_hash);
        total_kmers_in_hash += current_kmers_in_path;
        
        
        if (read_1->full_entry) {
                       
            kmer_hash->contaminated_kmers_per_read[current_kmers_in_path < MAX_READ_LENGTH? current_kmers_in_path:MAX_READ_LENGTH]++; 
            
            if (cmd_line->threashold == 0 ) {
                //Beware, long fasta entries won't be printed properly. 
                if (files_1 != NULL) {
                    print_sequnce_in_kmer_bin(read_1->seq,current_kmers_in_path, files_1);
                }
                if (files_2) {
                    print_sequnce_in_kmer_bin(read_2->seq,current_kmers_in_path, files_2);
                }
            }else{
                if (current_kmers_in_path < cmd_line->threashold) {
                    if (files_1 != NULL) {
                        print_sequnce_in_kmer_bin(read_1->seq,kmers_1->max_kmers, files_1);
                    }
                    if (files_2) {
                        print_sequnce_in_kmer_bin(read_2->seq,kmers_1->max_kmers, files_2);
                    }
                }else{
                    if ( cmd_line->print_contaminated == 1) {
                        if (files_1 != NULL) {
                            print_sequnce_in_kmer_bin(read_1->seq,kmers_1->max_kmers+1, files_1);
                        }
                        if (files_2) {
                            print_sequnce_in_kmer_bin(read_2->seq,kmers_1->max_kmers+1, files_2);
                        }
                    }
                }
                
            }
            
           
            
            
            total_length = 0;
            total_kmers_in_hash = 0;
        }else{
            log_and_screen_printf("Read %s is longer than the maximum supported length (%d)", read_1->seq->name, MAX_READ_LENGTH);
            errx(-1, "Unable to process a sequence (too long)");
        }
        
    }
    if (cmd_line->stats_file != NULL && strlen(cmd_line->stats_file) > 0) {
        kmer_hash_print_kmer_stats(cmd_line->stats_file,cmd_line->reference_kmers, cmd_line->input_read_1, kmer_hash);
    }
    
    
    char * hist_file = calloc(strlen(cmd_line->output_prefix)+20, sizeof(char));
    
    sprintf(hist_file, "%s.hist", cmd_line->output_prefix );
    kmer_hash_print_contaminated_kmers_histogram(hist_file, kmer_hash);
    free(hist_file);
        
    
    
    
    close_kmer_file_reader_wrapper(&read_1);
    if (read_2) {
        
        close_kmer_file_reader_wrapper(&read_2);
    }
    
    if(files_1 != NULL){
        free_filenames(files_1);
    }
    
    if (files_2) {
        free_filenames(files_2);
    }
    
    
    
}



int main(int argc, char **argv)
{
	setlocale (LC_ALL, "");
    
    log_and_screen_printf("\nkmer_filter.\n\n");    
	log_and_screen_printf(SVN_VERSION);
	log_and_screen_printf(SVN_COMMIT_DATE);
	log_and_screen_printf("Compiled on %s at %s \n\n", __DATE__, __TIME__);
    
    
    KmerFilterCmdLine cmd_line = parse_cmdline(argc, argv, sizeof(Element));
    
    
    KmerHash * kmer_hash = load_kmer_table(cmd_line);
    filter_reads(&cmd_line, kmer_hash);
    
    log_and_screen_printf("\nDONE");
    return 0;
    
}
