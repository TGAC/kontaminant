/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
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
#include <string.h>
#include <limits.h>
#include <getopt.h>
#include <cmd_line.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <binary_kmer.h>
#include <element.h>
#include <path.h>
#include <err.h>
#include <logger.h>

#ifdef ENABLE_MARK_PAIR
#include <db_graph.h>
#include <mark_pair.h>
#endif

const char *usage =
"\nusage: cortex_con [-h] [--input file_of_files] [--mem_height 14] [--dump_binary bin_output] [--input_format fastq|fasta|binary] [--output_contigs contigs.fa] \n" \
"by M. Caccamo (mario.caccamo@bbsrc.ac.uk) (primary contact for cortex_con), and Z. Iqbal (zam@well.ox.ac.uk)\n" \
"\n" \
"   [-h | --help] = This help screen.\n" \
"   [--input FILENAME] = File of filenames to be processed (start and end read is optional, format <filename>  <start read index>  <end read index> ).\n" \
"   [--kmer_size INT] = Kmer size (default 21), it has to be an odd number.\n" \
"   [--mem_width INT] = Size of hash table buckets (default 100).\n" \
"   [--mem_height INT] = Number of buckets in hash table in bits (default 10, this is a power of 2, ie 2^mem_height).\n" \
"   [--tip_clip INT] = Clips the tips in the graph, the argument defines the max length for the tips.\n" \
"   [--quality_score_threshold INT] = Filter for quality scores in the input file, any k-mer wiht a base wiht quality in the threshold or smaller is not considered (default 0).\n" \
"   [--remove_low_coverage_kmers INT] = Filter for kmers with coverage in the threshold or smaller.\n" \
"   [--dump_binary FILENAME] = Dump binary for graph in file (after applying all specified actions on graph).\n" \
"   [--ouput_supernodes FILENAME] = Fasta file with all the supernodes (after applying all specified actions on graph).\n" \
"   [--ouput_contigs FILENAME] = Fasta file with all the contigs (after applying all specified actions on graph).\n" \
"   [--input_format FORMAT] = File format for input (binary | fasta | fastq | hash ).\n" \
"   [--output_coverages] = Print coverages for contigs/supernodes in a different file with _cov suffix.\n" \
"   [--remove_seq_errors] = remove sequence of kmers induced by errors. Equivalent to --remove_low_coverage_kmers 1\n"\
"   [--remove_bubbles] = Removes the bubbles in the graph.\n"\
"   [--max_read_len] = Maximum read length over all input files.\n"\
"   [--quality_score_offset] = Fastq quality offset. Default 33. Use 63 for illumina.\n"\
"   [--quality_score_threshold ] = The minimiun phred value for a base to be considered in an assembly.\n"\
"   [--hash_output_file ] = Dumps the whole graph into a file. Read wiht the input_format hash. The file stores the information required to restore the hash table, hence mem_height and mem_width don't have any effect."
"\n";

 int default_opts(CmdLine * c)
{
    //core parameters
    c->verbose = false;
    c->kmer_size = 21;
    c->bucket_size = 100;
    c->number_of_buckets_bits = 10;
    
    //----------------
    //actions on/off
    //----------------
    //input
    c->input_file = false; //it is not present until we see it
    c->input_file_format_known = false;
    c->health_check_binary = false; 
	c->input_reference_known = false;
    c->output_kmer_coverage_known = false;
    
    //cleaning
    c->tip_clip = false;
    c->remove_bubbles = false;
    c->low_coverage_path_clip = false;
    c->low_coverage_node_clip = false;
    c->remove_low_coverage_supernodes = false;
    c->remove_spurious_links = false;

    
    //output
    c->dump_binary = false;
    c->output_fasta = false;
    c->detect_bubbles = false;
    c->output_coverages = false;
    c->dump_hash = false;
    c->graphviz = false;
    c->print_uncertain_as_n = true;
    c->output_log = false;
    c->output_reference_coverage_file_known=false;
    
    //read pairs
#ifdef ENABLE_READ_PAIR
    c->read_pair_enabled=false;
#endif
    
    //-----------
    //parameters
    //-----------
    //threads 
    c->threads = 1;
    c->max_double_y_complexity = PATH_MAX_DOUBLE_Y_COMPLEXITY;
    
    //input
    //c->input_file_format required
    c->quality_score_threshold = 0;
    //c->input_filename required
    //c->qual_filename required
    c->quality_score_offset = 33;
    c->max_read_len = 1000;
    //c->binary_version //TODO
    
    //cleaning
    //c->node_coverage_threshold required
    c->tip_length = 100; //TODO
    c->tip_clip_iterations = 100;
    //c->remove_low_coverage_supernodes_threshold required;
    //c->bubble_max_depth required
    //c->bubble_max_length required
    c->remove_spurious_links_min_coverage = 40;
    c->remove_spurious_links_max_difference = 1;
    
    //output
    //c->output_ctx_filename required
    //c->output_fasta_filename required
    //c->output_graphviz_filename required
    c->singleton_length = 100;
#ifdef ENABLE_BUBBLEPARSE
	c->algorithm = BUBBLES;
#else
	c->algorithm = PERFECT_PATH;
#endif	
	c->max_length=200000;
    c->min_subgraph_size=0;
    
	//read pairs
#ifdef ENABLE_READ_PAIR
	c->read_pair_distance = 100;
	c->read_pair_coverage = 20;
	c->read_pair_tolerance = 0;
	c->read_pair_min_bits = 2;
	c->read_pair_start_length = 1000;
	//c->read_pair_filename[0] required 
	c->read_pair_max_paths = 16;
	c->read_pair_min_kmers = 81+1; //TODO: check this
#endif 

#ifdef SOLID
    c->output_base_space = false;
#endif
	return 1;
}

CmdLine parse_cmdline(int argc, char *argv[], int unit_size)
{
    int i;
    printf("Command: ");
    for (i = 0; i < argc; i++) {
        printf("%s ", argv[i]);
    }
    printf("\n");
    printf("Unit size: %i\n", unit_size);
    
    CmdLine cmd_line;
    default_opts(&cmd_line);
    
    int opt;
    int longopt_index;
    
    static struct option long_options[] = {
        {"remove_bubbles", no_argument, NULL, 'a'},
        {"mem_width", required_argument, NULL, 'b'},
        {"tip_clip", required_argument, NULL, 'c'},
        {"output_contigs", required_argument, NULL, 'd'},
        {"output_coverages", no_argument, NULL, 'e'},
        {"output_supernodes", required_argument, NULL, 'f'},
        {"help", no_argument, NULL, 'h'},
        {"input", required_argument, NULL, 'i'},
        {"health_check",no_argument,NULL,'j'},
        {"kmer_size", required_argument, NULL, 'k'},
        {"log_file", required_argument, NULL, 'l'},
        {"min_subgraph_size", required_argument, NULL, 'm'},
        {"mem_height", required_argument, NULL, 'n'},
        {"dump_binary", required_argument, NULL, 'o'},
        {"quality_score_offset", required_argument, NULL, 'p'},
        {"quality_score_threshold", required_argument, NULL, 'q'},
        {"remove_low_coverage_supernodes", required_argument, NULL,'s'},
        {"input_format", required_argument, NULL, 't'},
        {"remove_seq_errors", no_argument, NULL, 'u'}, 
        {"verbose", no_argument, NULL, 'v'},
        {"detect_bubbles", required_argument, NULL, 'w'},
        {"singleton_length", required_argument, NULL, 'x'},
        {"remove_low_coverage_kmers", required_argument, NULL, 'z'},
        {"algorithm",required_argument,NULL,'A'},
#ifdef ENABLE_READ_PAIR
        {"read_pair_min_bits", required_argument, NULL, 'B'},
        {"read_pair_coverage", required_argument, NULL, 'C'},
        {"read_pair_distance", required_argument, NULL, 'D'},
        {"read_pair_max_paths", required_argument, NULL, 'E'},        
#endif        
#ifdef SOLID
        {"output_base_space", no_argument, NULL, 'F'},
#endif
        {"graphviz", required_argument, NULL, 'G'},
   		{"input_reference", required_argument, NULL, 'H'},
#ifdef ENABLE_READ_PAIR
        {"pair_info", required_argument, NULL, 'I'},
#endif        
		{"output_kmer_coverage", required_argument, NULL, 'J'},
#ifdef ENABLE_READ_PAIR
        {"read_pair_min_kmers", required_argument, NULL, 'K'},
#endif        
        {"remove_spurious_links",required_argument,NULL,'L'},
		{"output_reference_coverage_file", required_argument, NULL, 'M'},
#ifdef ENABLE_READ_PAIR
        {"print_stack_as_n", required_argument, NULL, 'N'},
#endif        
        {"hash_output_file", required_argument, NULL, 'O'},
        {"tip_clip_iterations", required_argument, NULL, 'P'},
#ifdef ENABLE_READ_PAIR
        {"read_pair_start_length", required_argument, NULL, 'S'},
#endif        
        {"threads", required_argument, NULL, 'T'},
        {0, 0, 0, 0}
    };
    
    while ((opt = getopt_long(argc, argv,
                              "ab:c:d:ef:hi:jk:l:n:o:p:q:s:t:uvw:x:z:A:B:C:D:E:FG:H:I:J:K:L:M:N:O:P:S:TZ:",
                              long_options, &longopt_index)) > 0)
    {                
        //Parse the default options
        switch (opt) {
            case 'a':
                cmd_line.remove_bubbles = true;
                break;
                
            case 'b':	//bucket size
                if (optarg == NULL)
                    errx(1, "[-b | --mem_width] option requires int argument [hash table bucket size]");
                cmd_line.bucket_size = atoi(optarg);
                if (cmd_line.bucket_size == 0 || cmd_line.bucket_size > SHRT_MAX)	//check that -b is not bigger than max_short 
                    errx(1, "[-b | --mem_width] option requires 'short' argument bigger than 0");
                break;
                
            case 'c':	//clip tip
                if (optarg == NULL)
                    errx(1,
                         "[-c | --clip_tip] option requires int argument [max length of tips]");
                cmd_line.tip_length = atoi(optarg);
                
                if (cmd_line.tip_length <= 0)
                    errx(1,
                         "[-c | --clip_tip] option requires int argument bigger than 0");
                cmd_line.tip_clip = true;
                break;
                
            case 'd': //dump paths from the tree decomposion of the graph
                cmd_line.output_fasta = true;
                if (optarg==NULL)
                    errx(1,"[-d | --output_contigs] option requires a filename");
                
                if (strlen(optarg)<LENGTH_FILENAME)
                {
                    strcpy(cmd_line.output_fasta_filename,optarg);
                }
                else
                {
                    errx(1,"[-d | --output_contigs] filename too long [%s]",optarg);
                }
                
                if (access(optarg,F_OK)==0){
                    errx(1,"[-d | --output_contigs] filename [%s] exist!",optarg);
                }
                
                cmd_line.algorithm=Y_WALK;
                break; 
                
            case 'e':	//output coverages
                cmd_line.output_coverages = true;
                break;
                
            case 'f':	//output of supernodes (more conservative paths)
                cmd_line.output_fasta = true;
                if (optarg == NULL)
                    errx(1,
                         "[-f | --output_supernodes] option requires a filename");
                
                if (strlen(optarg) < LENGTH_FILENAME) {
                    strcpy(cmd_line.output_fasta_filename, optarg);
                } else {
                    errx(1,"[-f | --output_supernodes] filename too long [%s]",optarg);
                }
                
                if (access(optarg, F_OK) == 0) {
                    errx(1,"[-f | --output_supernodes] filename [%s] exist!",optarg);
                }
				
                break;
                
            case 'h':
                printf("%s", usage);
                exit(0);
                break;
                
                
            case 'i': //file of filenames
                if (optarg == NULL)
                    errx(1,"[-i | --input] option requires a filename [file of filenames]");
                
                if (strlen(optarg) < LENGTH_FILENAME) {
                    strcpy(cmd_line.input_filename, optarg);
                } else {
                    errx(1, "[-i | --input] filename too long [%s]", optarg);
                }
                
                if (access(optarg, R_OK) == -1) {errx(1,"[-i | --input] filename [%s] cannot be accessed", optarg);
                }
                cmd_line.input_file = true;
                break;
                
            case 'j'://health check of binary 
                cmd_line.health_check_binary = true;
                break ;
                
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
                
            case 'l':
                if (optarg == NULL) {
                    errx(1, "[-l | --log_file] option requires a filename");
                    exit(-1);
                }
                if (strlen(optarg) < LENGTH_FILENAME) {
                    cmd_line.output_log = true;            
                    strcpy(cmd_line.log_filename, optarg);
                    log_start(optarg);
                    printf("Log file: %s\n", optarg);
                } else {
                    errx(1, "[-l | --log_file] filename too long");
                }                
                break;
                
            case 'n':	//number of buckets
                if (optarg == NULL)
                    errx(1,"[-n | --mem_height] option requires int argument [hash table number of buckets in bits]");
                cmd_line.number_of_buckets_bits = atoi(optarg);
                break;
                
            case 'o':	//output of binary ctx
                cmd_line.dump_binary = true;
                if (optarg == NULL)
                    errx(1, "[-o | --output_ctx] option requires a filename");
                if (strlen(optarg) < LENGTH_FILENAME) {
                    strcpy(cmd_line.output_ctx_filename, optarg);
                } else {
                    errx(1, "[-o | --output_ctx] filename too long [%s]",optarg);
                }
                if (access(optarg, F_OK) == 0) {
                    errx(1, "[-o | --output_ctx] filename [%s] exists!",optarg);
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
                
            case 's': //remove_low_coverage_supernodes
                cmd_line.remove_low_coverage_supernodes=true;
                if (optarg==NULL)
                    errx(1, "[-s | --remove_low_coverage_supernodes] option requires int argument [coverage threshold]");
                cmd_line.remove_low_coverage_supernodes_threshold= atoi(optarg);
                
                if (cmd_line.remove_low_coverage_supernodes_threshold <= 0)
                    errx(1, "[-s | --remove_low_coverage_supernodes] option requires intint argument bigger than 0");
                break;  
                
            case 't':
                if (optarg == NULL){
                    errx(1, "[-t  | --input_format binary | fasta | fastq | hash ] option requires a file type");
                    exit(-1);
                }
                
                cmd_line.input_file_format_known=true;
                
                if(strcmp(optarg, "binary") == 0){
                    cmd_line.input_file_format = CTX;
                }else if(strcmp(optarg, "fasta") == 0){
                    cmd_line.input_file_format = FASTA;
                }else if(strcmp(optarg, "fastq") == 0){
                    cmd_line.input_file_format = FASTQ;
                }else if(strcmp(optarg, "hash") == 0){
                    cmd_line.input_file_format = HASH;
                }else  {
                    fprintf(stderr, "[-t  | --input_format  binary | fasta | fastq | hash ] invalid option %s ]", optarg);
                    exit(-1);
                }
                break;
                
            case 'u': //remove_seq_errors --> this is consistent with an option on cortex_var -- it removes low coverate supernodes fixing the length tp k+1 and the coverage at 1
                cmd_line.remove_low_coverage_supernodes=true;
                cmd_line.remove_low_coverage_supernodes_threshold=1;
                break;        
                
            case 'v':
                cmd_line.verbose = true;
                break;
            case 'x':
                if (optarg == NULL) {
                    errx(1, "[-x | --singleton_length] option requires argument");
                }
                cmd_line.singleton_length = atoi(optarg);
                if(optarg == 0){
                    errx(1,  "[-x | --singleton_length] requires an integer greater than 0");
                }
                break;
                    
#ifdef ENABLE_BUBBLEPARSE
            case 'w':
                if (optarg == NULL) {
                    errx(1, "[-w | --detect_bubbles] option requires argument [max depth,max length]");
                }
                sscanf(optarg, "%d,%d", &cmd_line.bubble_max_depth, &cmd_line.bubble_max_length);
                if ((cmd_line.bubble_max_depth < 0) || (cmd_line.bubble_max_depth > 10)) {
                    errx(1, "[-w | --detect_bubbles] max depth out of range (0-10).");
                }
                if ((cmd_line.bubble_max_length < 1) || (cmd_line.bubble_max_length > 5000)) {
                    errx(1, "[-w | --detect_bubbles] max length out of range (1-5000).");
                }
                cmd_line.detect_bubbles = true;
                cmd_line.algorithm = BUBBLES;
                break;	
#endif
                
            case 'z':	//node coverage threshold
                if (optarg == NULL){
                    errx(1, "[-z | --node_coverage_threshold] option requires int argument [node coverate cut off]");
                }
                cmd_line.node_coverage_threshold = atoi(optarg);
                cmd_line.low_coverage_node_clip = true;
                
                if (cmd_line.node_coverage_threshold == 0){
                    errx(1, "[-z | --remove_low_coverage_kmers] option requires int argument bigger than 0");
                }
                break;
                
            case 'A':
#ifdef ENABLE_BUBBLEPARSE
                cmd_line.algorithm = BUBBLES;	
                errx(1,"[-A  | --algorithm ] option is not used for coloured bubbles");
                exit(-1);
                break;
#else
                if (optarg == NULL){
                    errx(1,"[-A  | --algorithm  BRANCHES | BUBBLES | PERFECT_PATH | Y_WALK | METACORTEX ] argument required]");
                    exit(-1);
                }
                
                if (strcmp(optarg, "BRANCHES")==0) {
                    cmd_line.algorithm = BRANCHES;
                } else if (strcmp(optarg, "BUBBLES")==0) {
                    cmd_line.algorithm = BUBBLES;
                } else if (strcmp(optarg, "PERFECT_PATH")==0) {
                    cmd_line.algorithm = PERFECT_PATH;
                } else if (strcmp(optarg, "Y_WALK")==0) {
                    cmd_line.algorithm = Y_WALK;
                } else if (strcmp(optarg, "METACORTEX")==0) {
                    cmd_line.algorithm = METACORTEX;                
                } else {
                    fprintf(stderr,
                            "[-A  | --algorithm  BRANCHES | BUBBLES | PERFECT_PATH | Y_WALK | METACORTEX ] recieved an invalid option %s ]", optarg);
                    exit(-1);
                }
#endif
                break;
                
            case 'B':
#ifdef ENABLE_READ_PAIR
                cmd_line.read_pair_enabled = true;
                cmd_line.read_pair_min_bits =  atoi(optarg);
#else
                fprintf(stderr, "[-B | --read_pair_min_bits] Cortex is not compiled to use read pair information, please compile with make ENABLE_READ_PAIR=1]");
                exit(-1);
#endif
                break;
                
            case 'C':
#ifdef ENABLE_READ_PAIR
                cmd_line.read_pair_enabled = true;
                cmd_line.read_pair_coverage =  atoi(optarg);
#else
                fprintf(stderr, "[-C | --read_pair_coverage] Cortex is not compiled to use read pair information, please compile with make ENABLE_READ_PAIR=1]");
                exit(-1);
#endif
                break;
                
            case 'D':
#ifdef ENABLE_READ_PAIR
                cmd_line.read_pair_enabled = true;
                cmd_line.read_pair_distance =  atoi(optarg);
#else
                fprintf(stderr, "[-D | --read_pair_distance] Cortex is not compiled to use read pair information, please compile with make ENABLE_READ_PAIR=1]");
                exit(-1);
#endif
                break;
                
            case 'E':
#ifdef ENABLE_READ_PAIR
                cmd_line.read_pair_enabled = true;
                cmd_line.read_pair_max_paths = atoi(optarg);
#else
                fprintf(stderr, "[-E | --read_pair_max_paths] Cortex is not compiled to use read pair information, please compile with make ENABLE_READ_PAIR=1]");
                exit(-1);
#endif
                break;
#ifdef SOLID
            case 'F':
                cmd_line.output_base_space = true;
                break;
#endif  
            
            case 'G':
                if (optarg == NULL) {
                    errx(1, "[-G | --graphviz] option requires a filename [file of filenames]");
                    exit(-1);
                }
                if (strlen(optarg) < LENGTH_FILENAME) {
                    strcpy(cmd_line.output_graphviz_filename, optarg);
                } else {
                    errx(1, "[-G | --graphviz] filename too long [%s]", optarg);
                }
                cmd_line.graphviz = true;
                printf("Graphviz file: %s\n", optarg);
                break;
			case 'H':
               	if (optarg == NULL) {
                   errx(1, "[-H | --input_reference ] option requires a filename [fasta file reference]");
                   exit(-1);
               }
               if (strlen(optarg) < LENGTH_FILENAME) {
                   strcpy(cmd_line.input_reference, optarg);
               } else {
                   errx(1, "[-H | --input_reference] filename too long [%s]", optarg);
               }
               if (access(optarg, R_OK) == -1) {
					errx(1,"[-H | --input_reference] filename [%s] cannot be accessed", optarg);
               }
               printf("Reference file: %s\n", optarg);
			   cmd_line.input_reference_known = true;
               break;
            case 'I':
#if defined (ENABLE_MARK_PAIR) || defined ( ENABLE_READ_PAIR)
                cmd_line.read_pair_enabled = true;
                if (strlen(optarg) < LENGTH_FILENAME) {
                    strcpy(cmd_line.read_pair_filename, optarg);
                } else {
                    errx(1, "[-I | --pair_info] filename too long [%s]", optarg);
                }
                
#ifdef ENABLE_MARK_PAIR
                mark_pair_validate_file(optarg);
#endif
#else
                fprintf(stderr, "[-I | --pair_info] Cortex is not compiled to use read pair information, please compile with make ENABLE_READ_PAIR=1]");
                exit(-1);
#endif
                break; 
            case 'J':
				cmd_line.output_kmer_coverage_know = true;
	            if (strlen(optarg) < LENGTH_FILENAME) {
	            	strcpy(cmd_line.output_kmer_coverage, optarg);
	             } else {
	                 errx(1, "[-J | --output_kmer_coverage] filename too long [%s]", optarg);
	             }
			break;
                
            case 'K':
#ifdef ENABLE_READ_PAIR
				
	            if (strlen(optarg) < LENGTH_FILENAME) {
	            	cmd_line.read_pair_min_kmers =  atoi(optarg);
	            } else {
	                errx(1, "[-K | --read_pair_min_kmers] filename too long [%s]", optarg);
	            }
                cmd_line.read_pair_enabled = true;
               
#else
                fprintf(stderr, "[-K | --read_pair_min_kmers] Cortex is not compiled to use read pair information, please compile with make ENABLE_READ_PAIR=1]");
                exit(-1);
#endif
                break;
                
            case 'L':
                if (optarg == NULL) {
                    errx(1, "[-L | --remove_spurious_links] option requires argument [max difference,min coverage]");
                } else {
                    sscanf(optarg, "%d,%d", &cmd_line.remove_spurious_links_max_difference, &cmd_line.remove_spurious_links_min_coverage);
                    if (cmd_line.remove_spurious_links_max_difference < 1) {
                        errx(1, "[-L | --remove_spurious_links] max difference must be 1 or more.");
                    } else if (cmd_line.remove_spurious_links_min_coverage < 1) {
                        errx(1, "[-L | --remove_spurious_links] min coverage must be 1 or more.");
                    } else {
                        log_and_screen_printf("Warning: Spurious links algorithm incomplete.\n");
                        cmd_line.remove_spurious_links = true;
                    }
                }
                break;
            case 'M':
				cmd_line.output_reference_coverage_file_known = true;
		        if (strlen(optarg) < LENGTH_FILENAME) {
		           	strcpy(cmd_line.output_reference_coverage_file, optarg);
		        } else {
               		errx(1, "[-M | --output_reference_coverage_file] filename too long [%s]", optarg);
		        }
			break;
				
            case 'N':
                cmd_line.print_uncertain_as_n = false;
                break;
                
            case 'O':
                cmd_line.dump_hash = true;
                if (strlen(optarg) < LENGTH_FILENAME) {
                    strcpy(cmd_line.output_hash_filename, optarg);
                } else {
                    errx(1, "[-O | --hash_output_file] filename too long [%s]", optarg);
                }
                break;
            case  'P':
                if (optarg == NULL)
                    errx(1,
                         "[-P | --tip_clip_iterations] option requires int argument ");
                cmd_line.tip_clip_iterations = atoi(optarg);
                
                if (cmd_line.tip_length <= 0)
                    errx(1,
                         "[-P | --tip_clip_iterations] option requires int argument bigger than 0");
                
                break;
            case 'Q':
                errx(1, "454 quality files not implemented yet. Use the regular fasta format. \n");
                if (optarg == NULL){
                    errx(1, "[-Q  | --qual_file] option requires a filename [file of filenames]");
                    exit(-1);
                }
                
                printf("Qual file %s\n", optarg);
                if (strlen(optarg) < LENGTH_FILENAME) {
                    strcpy(cmd_line.qual_filename, optarg);
                } else {
                    errx(1, "[-Q | --qual_file] filename too long [%s]", optarg);
                }
                
                if (access(optarg, R_OK) == -1) {
                    errx(1, "[-Q | --qual_file ] filename [%s] cannot be accessed", optarg);
                }
                break;
                
            case 'R':
                errx(1, "454 quality files not implemented yet. Use the regular fasta format. \n");
                cmd_line.input_file_format_known = true;
                cmd_line.input_file_format = ROCHE;
                printf ("Doing 454\n");
                break;
                
            case 'S':
#ifdef ENABLE_READ_PAIR
                cmd_line.read_pair_enabled = true;
                cmd_line.read_pair_start_length =  atoi(optarg);
#else
                fprintf(stderr, "[-S | --read_pair_start_length] Cortex is not compiled to use read pair information, please compile with make ENABLE_READ_PAIR=1]");
                exit(-1);
#endif
                break;
                
            case 'T':	//Number of threads
                if (optarg == NULL) {
                    errx(1, "[-T | --threads INT] option requires int argumen");
                }
                cmd_line.threads = atoi(optarg);
                if (cmd_line.threads <= 0) {
                    errx(1, "[-T | --threads INT] option requires int argument bigger than 0");
                }                
                int tmp = cmd_line.threads  - 1;
                int check = cmd_line.threads & tmp;	     
                if (check != 0 ) {
                    errx(1, "[-T | --threads INT] has to be a power of two. ");
                }                
                break;
                
            case 'Z'://max read length --> use by the parser to read fasta/fastq files.                 
                if (optarg==NULL) {
                    errx(1,"[-Z | --max_read_len] option requires int argument [maximum read length in input]");
                }
                
                cmd_line.max_read_len = atoi(optarg);
                
                if (cmd_line.max_read_len <= 0)
                    errx(1,"[-Z | --max_read_len] option requires int argument bigger than 0");
                break;

            default:
                // Error already displayed by getopt_long
                exit(1);
                break;
        }
    }

    // We've already sent this to the screen, but the log wasn't open then...
    if (cmd_line.output_log) {
#ifdef ENABLE_BUBBLEPARSE
        log_printf("Cortex Con Bubble Finding\n\n");
#else
        log_printf("\nCortex Con\n\n");
#endif
        log_printf("Command: ");
        for (i = 0; i < argc; i++) {
            log_printf("%s ", argv[i]);
        }
        log_printf("\nUnit size: %i\n", unit_size);
    }    
    
    //check if input file format is known
    if (!cmd_line.input_file_format_known) {
        errx(1, "file format not defined [ --input_format fasta | fastq | binary | hash]");
    }
    //check input file
    if (!cmd_line.input_file) {
        errx(1, "input file required [option -i | --input_file]");
    }
    
    if (cmd_line.number_of_buckets_bits < 0) {
        printf("args error -b %i -h %i \n", cmd_line.bucket_size,cmd_line.number_of_buckets_bits);
        errx(1, "memory configuration erorr - revise options -b -h and/or -m");
    }
    
    if(cmd_line.input_file_format == ROCHE){
        if(strlen(cmd_line.qual_filename) == 0){
            printf("454 format requires a quality file [option -Q | --qual_file FILENAME] \n");
        }
    }
    
    //check kmer_size is odd
    if (cmd_line.kmer_size % 2 == 0) {
        errx(1, "[-k | --kmer_size] is even [%i]!", cmd_line.kmer_size);
    }
    
#ifdef ENABLE_BUBBLEPARSE
    if ((cmd_line.detect_bubbles != true) || (cmd_line.algorithm != BUBBLES)) {
        errx(1, "Cortex_bub requires the -w option.");
    }
#endif

    
    //check that if we are printing coverages then we are dumping contigs
    //RHRG: We dont need to do t his since we may want to just print the coverages of the full graph. 
    //if (cmd_line.output_coverages == true){
    //   if (cmd_line.output_fasta == false){
    //       errx(1, "[-e] cannot output coverages without contigs\n");
    //    } 
    //}
    //check that if we are going to print the coverage  of a reference, the reference file exists;
	if(cmd_line.output_reference_coverage_file_known){
		if(!cmd_line.input_reference_known){
			errx(1, "[-M | --output_reference_coverage ] <filename> requires to give a reference file with [-H | --input_reference ] <filename>\n");
		}
	}

    return cmd_line;
}
