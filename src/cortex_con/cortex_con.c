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

/*
  dB_graph.c - implementation
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <locale.h>
#ifdef THREADS
#include <pthread.h>
#endif
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

#ifdef ENABLE_READ_PAIR
#include <binary_tree.h>
#include <read_pair.h>
#endif

#ifdef ENABLE_MARK_PAIR
#include <binary_tree.h>
#include <mark_pair.h>
#endif


#ifdef ENABLE_BUBBLEPARSE
#include <bubble_find.h>
#endif

#include <cleaning.h>
#include <cmd_line.h>
#include "metacortex.h"


typedef struct {
    char* filename;
    int colour;
    int pair;
} ReadFileDescriptor;

#define FILE_LIST_SIZE 1024

void write_graphviz_file(char *filename, dBGraph * db_graph);
void timestamp();

char *remove_file_extension(char *filename)
{
	int end = strlen(filename) - 1;
	int start = end - 6;
	int i;

	if (DEBUG)
		printf("[remove_file_extension] In:  %s\n", filename);

	if (start < 0)
		start = 0;

	for (i = end; i >= start; i--) {
		if (filename[i] == '.') {
			filename[i] = 0;
			break;
		}
	}

	if (DEBUG)
		printf("[remove_file_extension] Out: %s\n", filename);

	return filename;
}

#ifdef ENABLE_READ_PAIR
void load_pair_data(char* filename, ReadPairDescriptorArray * rpda, int fastq_ascii_offset,dBGraph * db_graph) {
    FILE *fp_fof = fopen(filename, "r");
    FILE *fp_1;
    FILE *fp_2;
    char file_1[1000];
    char file_2[1000];
    int read_length = 0;
    int read_tolerance = 0;
    int insert_size = 0;
    int allowed_misses = 0;
    int min_pair_coverage = 0;
    int max_pair_coverage = 0;
    int depth;
    
    if (!fp_fof) {
        printf("Error: can't open file %s\n", filename);
        exit(-1);
    }

    while (!feof(fp_fof)) {
        fscanf(fp_fof, "%s %s %d %d %d %d %d %d %d\n", file_1, file_2, &insert_size, &read_length, &read_tolerance, &allowed_misses, &min_pair_coverage, &max_pair_coverage, &depth);

        if (strcmp(file_1, file_2) == 0) {
            log_and_screen_printf("Error: Both files in the pair are the same!\n");
            fclose(fp_fof);
            exit(-1);
        }
        if (insert_size<(read_length*2)) {
            log_and_screen_printf("Error: Bad insert size (%'d).\n", insert_size);
            fclose(fp_fof);
            exit(-1);
        }
        
        if (read_length < 1) {
            log_and_screen_printf("Error: Bad read length (%'d).\n", read_length);
            fclose(fp_fof);
            exit(-1);
        }
        if (min_pair_coverage < 1) {
            log_and_screen_printf("Error: Bad minimum pair coverage(%'d).\n", min_pair_coverage);
            fclose(fp_fof);
            exit(-1);
        }
        if (max_pair_coverage < 1) {
            log_and_screen_printf("Error: Bad maximum read pair coverage (%'d).\n", max_pair_coverage);
            fclose(fp_fof);
            exit(-1);
        }
        
        if (depth < 1) {
            log_and_screen_printf("Error: Bad depth (%'d).\n", depth);
            fclose(fp_fof);
            exit(-1);
        }

        if ((read_tolerance < 0) || (read_tolerance > insert_size)) {
            log_and_screen_printf("Error: Bad tolerance (%'d).\n", read_tolerance);
            fclose(fp_fof);
            exit(-1);
        }
        
        if ((allowed_misses < 0) || (allowed_misses > read_length)) {
            log_and_screen_printf("Error: Bad allowed misses (%'d).\n", allowed_misses);
            fclose(fp_fof);
            exit(-1);
        }

       

        fp_1 = fopen(file_1, "r");
        if (!fp_1) {
            log_and_screen_printf("Error: can't open file %s\n", file_1);
            fclose(fp_fof);
            exit(-1);
        }
        
        fp_2 = fopen(file_2, "r");
        if (!fp_2) {
            log_and_screen_printf("Error: can't open file %s\n", file_2);
            fclose(fp_fof);
            fclose(fp_1);
            exit(-1);
        }
        
        read_pair_descriptor_array_add_pair(0, 0, 0, insert_size, read_length, read_tolerance, min_pair_coverage, max_pair_coverage,  allowed_misses, depth,  &simple_iterated_score_paths, rpda);


        fflush(stdout);
	
        fclose(fp_1);
        fclose(fp_2);        
    }
#ifdef THREADS
    hash_table_threaded_traverse(&db_node_action_clear_flags, db_graph);
#else
    hash_table_traverse(&db_node_action_clear_flags, db_graph);
#endif
    read_pair_mark_contiguous_perfect_paths(read_pair_get_maximum_insert_size(rpda) , db_graph);//TODO: bring back several read pairsr
    
    rewind(fp_fof);
    while (!feof(fp_fof)) {
        fscanf(fp_fof, "%s %s %d %d %d %d %d %d %d\n", file_1, file_2, &insert_size, &read_length, &read_tolerance, &allowed_misses,&min_pair_coverage, &max_pair_coverage, &depth);
        
        printf("Enriching pair: %s and %s, insert size %'d read length %'d tolerance %'d allowed misses %'d\n", file_1, file_2, insert_size, read_length, read_tolerance, allowed_misses);
       
        
        fp_1 = fopen(file_1, "r");
        if (!fp_1) {
            log_and_screen_printf("Error: can't open file %s\n", file_1);
            fclose(fp_fof);
            exit(-1);
        }
        
        fp_2 = fopen(file_2, "r");
        if (!fp_2) {
            log_and_screen_printf("Error: can't open file %s\n", file_2);
            fclose(fp_fof);
            fclose(fp_1);
            exit(-1);
        }
        
        
        
        fflush(stdout);
		
        //read_pair_enrich_graph(fp_1, fp_2, rpda->pair[0], fastq_ascii_offset,db_graph);
        read_pair_search_enrich_graph(fp_1, fp_2, rpda->pair[0], fastq_ascii_offset,db_graph);
        
        fclose(fp_1);
        fclose(fp_2);        
    }
    
    
    
    db_graph->rpda = rpda;

    fclose(fp_fof);
}

#elif ENABLE_MARK_PAIR

void load_pair_data(char* filename,  int fastq_ascii_offset,dBGraph * db_graph) {
    FILE *fp_fof = fopen(filename, "r");
    FILE *fp_1;
    FILE *fp_2;
    char file_1[1000];
    char file_2[1000];
        
    if (!fp_fof) {
        printf("Error: can't open file %s\n", filename);
        exit(-1);
    }
    
    
#ifdef THREADS
    hash_table_threaded_traverse(&db_node_action_clear_flags, db_graph);
#else
    hash_table_traverse(&db_node_action_clear_flags, db_graph);
#endif
    mark_pair_mark_contiguous_perfect_paths(db_graph);
    
    rewind(fp_fof);
    while (!feof(fp_fof)) {
        fscanf(fp_fof, "%s %s\n", file_1, file_2);
        
        printf("Enriching pair: %s and %s", file_1, file_2);
        
        
        fp_1 = fopen(file_1, "r");
        if (!fp_1) {
            log_and_screen_printf("Error: can't open file %s\n", file_1);
            fclose(fp_fof);
            exit(-1);
        }
        
        fp_2 = fopen(file_2, "r");
        if (!fp_2) {
            log_and_screen_printf("Error: can't open file %s\n", file_2);
            fclose(fp_fof);
            fclose(fp_1);
            exit(-1);
        }
        
        
        
        fflush(stdout);
		
        //read_pair_enrich_graph(fp_1, fp_2, rpda->pair[0], fastq_ascii_offset,db_graph);
       // read_pair_search_enrich_graph(fp_1, fp_2, rpda->pair[0], fastq_ascii_offset,db_graph);
        
        fclose(fp_1);
        fclose(fp_2);        
    }
    
    
    fclose(fp_fof);
}

#endif

boolean add_read_file(char* filename, int colour, int pair, int* n_files, ReadFileDescriptor* list[])
{
    boolean already_seen = false;
    FILE *fp;
    int i;
    
    // Check if file already in list
    for (i=0; i<*n_files; i++) {
        if (strcmp(filename, list[i]->filename) == 0) {
            already_seen = true;
            printf("Error: file %s is included more than once.\n", filename);
            exit(-1);
        }
    }
    
    // Check file exists
    fp = fopen(filename, "r");
    if (fp) {
        fclose(fp);
    } else {
        printf("Error: can't open file %s\n", filename);
        exit(-1);
    }
    
    // Add to list if not already seen
    if (!already_seen) {
        if (*n_files == FILE_LIST_SIZE) {
            printf("Error: maximum file list size reached. Continuing might be dangerous.\n");
            exit(-1);
        } else {            
            list[*n_files] = malloc(sizeof(ReadFileDescriptor));
            if (!list[*n_files]) {
                printf("Error: Couldn't allocate memory for file descriptor\n");
                exit(-1);
            }
            
            list[*n_files]->filename = malloc(strlen(filename)+1);
            if (!list[*n_files]->filename) {
                printf("Error: Couldn't allocate memory for filename\n");
                exit(-1);
            }
            strcpy(list[*n_files]->filename, filename);
            list[*n_files]->colour = colour;
            list[*n_files]->pair = pair;
            
            *n_files = *n_files + 1;
        }
    }
    
    return already_seen;
}

void output_basic_info(CmdLine cmd_line)
{
    log_and_screen_printf("Max k: %i\n", (NUMBER_OF_BITFIELDS_IN_BINARY_KMER*32)-1);
    if (cmd_line.input_file_format == FASTQ) { 
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
               

int main(int argc, char **argv)
{
	setlocale (LC_ALL, "");
    FILE *fp_fnames;
    int hash_key_bits;
    dBGraph *db_graph = NULL;
    short kmer_size;
    int bucket_size = 0;
    CmdLine cmd_line;
    boolean using_colours = false;
    ReadFileDescriptor* file_list[FILE_LIST_SIZE];
    long long seq_length = 0;
    int n_file_list = 0;
    int i;
    
#ifdef ENABLE_BUBBLEPARSE
    log_and_screen_printf("Cortex Con Bubble Finding\n\n");
#else
    log_and_screen_printf("\nCortex Con\n\n");
#endif
        
	log_and_screen_printf(SVN_VERSION);
	log_and_screen_printf(SVN_COMMIT_DATE);
	log_and_screen_printf("Compiled on %s at %s \n\n", __DATE__, __TIME__);
    //command line arguments
    cmd_line = parse_cmdline(argc, argv, sizeof(Element));

    output_basic_info(cmd_line);
    
#ifdef ENABLE_READ_PAIR
    log_and_screen_printf("Read pair enabled, %lu signatures on %i segments\n", READ_PAIR_LENGTH, READ_PAIR_SEGMENTS);
#endif
    
#ifdef ENABLE_MARK_PAIR
    log_and_screen_printf("Read pair enabled\n");
#endif
        
    if(cmd_line.input_file_format == HASH) {
        db_graph = hash_table_read_dumped_memory(cmd_line.input_filename);        
    } else {
        fp_fnames = fopen(cmd_line.input_filename, "r");	//open file of file names
        kmer_size = cmd_line.kmer_size;
        hash_key_bits = cmd_line.number_of_buckets_bits;	//number of buckets: 2^hash_key_bits
        bucket_size = cmd_line.bucket_size;
       // DEBUG = cmd_line.verbose;
        
        timestamp();
        log_and_screen_printf("\nInput file of filenames: %s\n", cmd_line.input_filename);
        log_and_screen_printf("Kmer size: %'d Hash table size (%'d bits): %'d Hash table bucket size: %'d Total size: %'qd\n",
                              cmd_line.kmer_size,
                              cmd_line.number_of_buckets_bits,
                              1 << cmd_line.number_of_buckets_bits,
                              cmd_line.bucket_size,
                              ((long long)1 << cmd_line.number_of_buckets_bits) * cmd_line.bucket_size);
        fflush(stdout);
        
        //Create the de Bruijn graph/hash table
        db_graph = hash_table_new(cmd_line.number_of_buckets_bits,
                                  cmd_line.bucket_size, 25, cmd_line.kmer_size);
        
        if (db_graph == NULL) {
            printf("Please free up memory and re-run.\n");
            exit(-1);
        }
        
        //TODO
        if (cmd_line.threads > 0) {
        	db_graph->number_of_threads = cmd_line.threads;
        }
		
        //Zone to initalise the buffers
        path_array_initialise_buffers(kmer_size);
        timestamp();
        log_and_screen_printf("\nTable created: %'d\n", 1 << cmd_line.number_of_buckets_bits);
        db_graph_print_status(db_graph);
        fflush(stdout);
        int count_file = 0;
        long long total_length = 0;	//total sequence length
        
        //Go through all the files, loading data into the graph
        
        boolean all_entries_are_unique = false;
        
        if (cmd_line.input_file_format == CTX) {
            all_entries_are_unique = true;
        }
        
        long long bad_reads = 0;
        
        while (!feof(fp_fnames)) {
            short colour = 0;
            short pair = 0;
            char filename[1024];
            char line[1024];
            count_file++;

            line[0] = 0;
            if (fgets(line, 1024, fp_fnames)) {
                if (strlen(line) > 1) {
                    sscanf(line, "%s %hd %hd\n", filename, &colour, &pair);
                    add_read_file(filename, colour, pair, &n_file_list, file_list);
                }
            }            
        }
        
        fclose(fp_fnames);
        
        for (i=0; i<n_file_list; i++) {
            short colour = file_list[i]->colour;
            char* filename = file_list[i]->filename;
            

            if (DEBUG) {
                log_and_screen_printf("\nNew file: %s colour %'d\n", filename, colour);
                fflush(stdout);
            }
            
            fflush(stdout);
            
            if ((colour < 0) || (colour >= NUMBER_OF_COLOURS)) {
                printf("Colour is out of allowable range (maximum number of colours: %'d).", NUMBER_OF_COLOURS);
                fflush(stdout);
                exit(-1);
            }
            
            if (colour > 0) {
                using_colours = true;
            }
            
            switch (cmd_line.input_file_format) {
                    
                case CTX:                    
                    log_and_screen_printf("\nReading ctx file %'d: %s\n", i+1, filename);		
				    fflush(stdout);
                    seq_length = load_binary_from_filename_into_graph(filename,
                                                         db_graph,
                                                         colour,
                                                         all_entries_are_unique);
                    
                    all_entries_are_unique = false;
                    break;
                    
                case FASTQ:
                    if (colour != 0) {
                        printf("\n*** Warning: Colour should be 0 for FASTQ files.\n");
                    }
                    log_and_screen_printf("\nReading fastq file %'d: %s\n", i+1, filename);		
				    fflush(stdout);
                    seq_length = load_fastq_from_filename_into_graph(filename,
                                                        colour,
                                                        &bad_reads,
                                                        cmd_line.quality_score_threshold,
                                                        5000, cmd_line.quality_score_offset, db_graph);
                    break;
                    
                case FASTA:
                    if (colour != 0) {
                        printf("\n*** Warning: Colour should be 0 for FASTA files.\n");
                    }
                    log_and_screen_printf("\nReading fasta file %'d: %s\n", i+1, filename);		
				    fflush(stdout);
                    seq_length = load_fasta_from_filename_into_graph(filename,
                                                        colour,
                                                        &bad_reads,
                                                        5000, db_graph);
                    break;
                default:
                    printf("Unknown file format. ");
                    fflush(stdout);
                    break;
            }

            total_length += seq_length;
            timestamp();
            log_and_screen_printf("\nRead of file %'d complete. Total kmers: %'lld Bad reads: %'qd Seq length: %'qd Total seq length: %'qd\n\n",
                                  i+1, hash_table_get_unique_kmers(db_graph), bad_reads, seq_length, total_length);
            hash_table_print_stats(db_graph);
            
            fflush(stdout);
            
                       
        }
    }
	path_array_initialise_buffers(db_graph->kmer_size);
    db_graph->max_double_y_complexity = cmd_line.max_double_y_complexity;
	fflush(stdout);
		
	if (cmd_line.low_coverage_node_clip) {
		timestamp();
		log_and_screen_printf("\nRemove low coverage nodes (<=%'d) \n", cmd_line.node_coverage_threshold);
		fflush(stdout);
		cleaning_remove_low_coverage(cmd_line.node_coverage_threshold, db_graph);
        //		db_graph_remove_low_coverage_nodes
        //		    (cmd_line.node_coverage_threshold, db_graph);
        hash_table_print_stats(db_graph);
		
	}
    
    if (cmd_line.remove_low_coverage_supernodes) {
		timestamp();
		log_and_screen_printf("\nRemoving paths with low coverage...\n");
		fflush(stdout);
		int p = cleaning_prune_low_coverage_path(cmd_line.remove_low_coverage_supernodes_threshold, cmd_line.tip_length, db_graph);
		log_and_screen_printf("%'d nodes removed\n", p);
        hash_table_print_stats(db_graph);
	}
	
	if (cmd_line.tip_clip) {
		timestamp();
		log_and_screen_printf("\nClip tips\n");
		fflush(stdout);
        log_and_screen_printf("%'d tips clipped\n", cleaning_remove_tips(cmd_line.tip_length, cmd_line.tip_clip_iterations ,db_graph));
        hash_table_print_stats(db_graph);
	}
    
    if (cmd_line.remove_spurious_links) {
		timestamp();
		log_and_screen_printf("\nRemoving spurious links\n");
		fflush(stdout);
		int links_removed = cleaning_remove_spurious_links(cmd_line.remove_spurious_links_max_difference, cmd_line.remove_spurious_links_min_coverage, db_graph);
		log_and_screen_printf("%'d links removed\n", links_removed);
        hash_table_print_stats(db_graph);
        if (cmd_line.tip_clip) {
            timestamp();
            log_and_screen_printf("\nClip tips\n");
            fflush(stdout);
            log_and_screen_printf("%'d tips clipped\n", cleaning_remove_tips(cmd_line.tip_length, cmd_line.tip_clip_iterations ,db_graph));
            hash_table_print_stats(db_graph);
        }
    }
    	
	if (cmd_line.remove_bubbles) {
		timestamp();
		log_and_screen_printf("\nRemoving bubbles\n");
		fflush(stdout);
		cmd_line.bubble_max_length = db_graph->kmer_size * 10 + 1;
		//TODO: we should defined the bubble depth here too
		int kmers_removed_1 = cleaning_remove_bubbles(cmd_line.bubble_max_length, db_graph);
		log_and_screen_printf("%'d kmers removed\n", kmers_removed_1);	
        hash_table_print_stats(db_graph);
        if (cmd_line.tip_clip) {
            timestamp();
            log_and_screen_printf("\nClip tips\n");
            fflush(stdout);
            log_and_screen_printf("%'d tips clipped\n", cleaning_remove_tips(cmd_line.tip_length, cmd_line.tip_clip_iterations ,db_graph));
            hash_table_print_stats(db_graph);
        }
    }
    
        
    if (cmd_line.output_kmer_coverage_known) {
        timestamp();
        char fname[1024];
        strcpy(fname, cmd_line.output_kmer_coverage);
        strcat(fname, ".tab");
        
        FILE *out = fopen(fname, "w");
        db_graph_print_coverage(out, db_graph);
        
        strcpy(fname, cmd_line.output_kmer_coverage);
        strcat(fname, "_kmer.tab");
        
        FILE *out2 = fopen(fname, "w");
        db_graph_print_kmer_coverage(out2, db_graph);
        fclose(out2);
        fclose(out);
        fflush(stdout);
    }

	if(cmd_line.output_reference_coverage_file_known){
		timestamp();
		log_and_screen_printf("\nPrinting coverages from fasta\n");
		fflush(stdout);
	}
    
    if(cmd_line.dump_hash){
		timestamp();
		log_and_screen_printf("\nDumping hash table to file: %s\n", cmd_line.output_hash_filename);
		hash_table_dump_memory(cmd_line.output_hash_filename, db_graph);
		fflush(stdout);
	}
    hash_table_print_stats(db_graph);
    if (cmd_line.dump_binary) {
        timestamp();
        log_and_screen_printf("\n");
        // If using colours, we dump a separate file for each colour. Otherwise, it's just one.
        if (!using_colours) {
            log_and_screen_printf("Dumping graph: %s\n", cmd_line.output_ctx_filename);
            db_graph_dump_binary(cmd_line.output_ctx_filename, &db_node_check_flag_not_pruned, db_graph);
            fflush(stdout);
        } else {
            int c;
            for (c = 0; c < NUMBER_OF_COLOURS; c++) {
                char fname[1024];
                char suffix[64];
                char *extension;
                
                sprintf(suffix, "_c%d.ctx", c);
                strcpy(fname, cmd_line.output_ctx_filename);
                extension = strstr(fname, ".ctx");
                if (!extension) {
                    extension = strstr(fname, ".CTX");
                }
                
                if (extension) {
                    strcpy(extension, suffix);
                } else {
                    strcat(fname, suffix);
                }
                
                log_and_screen_printf("Dumping graph %s - ", fname);
                db_graph_dump_binary_by_colour(fname, &db_node_check_flag_not_pruned, c, db_graph);
            }
            fflush(stdout);
        }
    }
    
    if (cmd_line.output_fasta) {
        
#ifdef SOLID
        db_graph->output_base_space = cmd_line.output_base_space;
#endif
      
#ifdef ENABLE_BUBBLEPARSE
        // Do bubble detection and output
         timestamp();
        log_and_screen_printf("\nFinding bubbles to output...\n");
        fflush(stdout);
        int max_length = cmd_line.bubble_max_length + 1000;
        
#ifdef DEBUG_CLEANUP
        db_graph_cleanup_graph(db_graph);
#endif
        
        // Identify branches
        db_graph_identify_branches(max_length, db_graph);
        
        // Check for .fa, .fasta or .fastq at the end of filename and remove
        remove_file_extension(cmd_line.output_fasta_filename);
        
        assert(cmd_line.bubble_max_length > 0);
        
        // Walk branches
        db_graph_walk_branches(cmd_line.output_fasta_filename,
                               max_length,
                               cmd_line.bubble_max_length,
                               cmd_line.bubble_max_depth,
                               db_graph);
#elif ENABLE_MARK_PAIR
        if (cmd_line.read_pair_filename[0] != 0) {
            // Load pair file and enrich graph
            
            timestamp();
            fflush(stdout);
            //load_pair_data(cmd_line.read_pair_filename, cmd_line.quality_score_offset,db_graph);
            mark_pair_load_pair_data(cmd_line.read_pair_filename, cmd_line.quality_score_offset,db_graph);
            fflush(stdout);
            log_and_screen_printf("Cleaning flags...\n");
#ifdef THREADS
            hash_table_threaded_traverse(&db_node_action_clear_flags, db_graph);
#else
            hash_table_traverse(&db_node_action_clear_flags, db_graph);
#endif
            mark_pair_print_paths(cmd_line.output_fasta_filename, cmd_line.max_length, cmd_line.output_coverages, cmd_line.graphviz, db_graph);
        }   
        
//        read_pair_print_paths_single_walk(cmd_line.output_fasta_filename, cmd_line.max_length, cmd_line.output_coverages, cmd_line.graphviz,rpda, db_graph);
        log_and_screen_printf("Finished printing read_pair paths.\n");
#elif ENABLE_READ_PAIR
        if (cmd_line.read_pair_enabled) {
            ReadPairDescriptorArray * rpda = new_read_pair_descriptor_array(10,
                                                                            cmd_line.read_pair_distance,
                                                                            cmd_line.read_pair_coverage,
                                                                            cmd_line.read_pair_min_bits,
                                                                            cmd_line.read_pair_start_length,
                                                                            cmd_line.read_pair_max_paths,
                                                                            cmd_line.read_pair_min_kmers,
                                                                            cmd_line.print_uncertain_as_n);
            
            fflush(stdout);
            db_graph->rpda = rpda;
            if (cmd_line.read_pair_filename[0] != 0) {
                // Load pair file and enrich graph
                
                timestamp();
                fflush(stdout);
                load_pair_data(cmd_line.read_pair_filename, rpda, cmd_line.quality_score_offset,db_graph);
                fflush(stdout);
            }
            
            if (rpda->number_of_pairs > 0) {
                timestamp();
                fflush(stdout);
                log_and_screen_printf("\nDumping contigs from read pairs to %s\n", cmd_line.output_fasta_filename);
                clear_visited_count();
               // read_pair_print_paths(cmd_line.output_fasta_filename, cmd_line.max_length, cmd_line.output_coverages, cmd_line.graphviz,rpda, db_graph); //TODO: have different calls for the different algorithms
#ifdef THREADS
                hash_table_threaded_traverse(&db_node_action_clear_flags, db_graph);
#else
                hash_table_traverse(&db_node_action_clear_flags, db_graph);
#endif
                read_pair_print_paths_single_walk(cmd_line.output_fasta_filename, cmd_line.max_length, cmd_line.output_coverages, cmd_line.graphviz,rpda, db_graph);
                log_and_screen_printf("Finished printing read_pair paths.\n");
            } else {
                // Should probably just display a warning and call y_walk, but useful for now
                printf("Error: no read pair descriptors. Can't continue.\n");
                exit(-1);
            }
            
            destroy_read_pair_descriptor_array(&rpda);
        }
#else
        switch (cmd_line.algorithm) {
            case PERFECT_PATH:
                log_and_screen_printf("\nDumping supernodes: %s\n", cmd_line.output_fasta_filename);
                fflush(stdout);
                
                perfect_path_print_paths(cmd_line.output_fasta_filename,
                                         cmd_line.max_length, cmd_line.singleton_length, 
                                         cmd_line.output_coverages,
                                         db_graph);
                
                break;
            
            case Y_WALK:
                log_and_screen_printf("\nDumping supernodes (Y_WALK): %s\n", cmd_line.output_fasta_filename);
                fflush(stdout);
                y_walk_print_paths(cmd_line.output_fasta_filename,
                                   cmd_line.max_length,cmd_line.singleton_length, 
                                   cmd_line.output_coverages, cmd_line.graphviz, db_graph);
                log_and_screen_printf("Dump complete\n");
                break;	
            case BRANCHES:
                log_and_screen_printf("\nDumping supernodes (branches): %s\n", cmd_line.output_fasta_filename);
                fflush(stdout);
                //branches_get_path(cmd_line.output_fasta_filename, cmd_line.max_length, cmd_line.output_coverages, db_graph);
                break;
            case METACORTEX:
                log_and_screen_printf("\nDumping subgraph consensus contigs: %s\n", cmd_line.output_fasta_filename);
                metacortex_find_subgraphs(db_graph, cmd_line.output_fasta_filename, cmd_line.min_subgraph_size);
                break;
            default:
                log_and_screen_printf("Algorithm not implmeneted \n");
                break;
        }
#endif
    }
    
    if (cmd_line.graphviz) {
        timestamp();
        // Write graphviz file
        
        write_graphviz_file(cmd_line.output_graphviz_filename, db_graph);
        fflush(stdout);
	}
	
    // Stop program terminating, so XCode leaks tool can report!
    //printf("press char to continue");
    //char c = getchar();       
    //printf("Character pressed %c\n", c);
    
    timestamp(); 
    log_and_screen_printf("\n\nDONE\n\n");
    
    return 0;
}

void print_graph(dBGraph * db_graph)
{
	void print_graphviz(dBNode * node) {
		if (node != NULL) {
			BinaryKmer tmp, co;

			short kmer_size = db_graph->kmer_size;
			binary_kmer_assignment_operator(tmp, node->kmer);
			binary_kmer_reverse_complement(&tmp, kmer_size, &co);
			char seq[kmer_size], seqNext[kmer_size],
			    seq1[kmer_size];
			binary_kmer_to_seq(&tmp, kmer_size, seq1);
			char *print = db_node_check_for_any_flag(node,
								 STARTING_FORWARD
								 |
								 BRANCH_NODE_FORWARD
								 |
								 BRANCH_NODE_REVERSE
								 |
								 END_NODE_FORWARD
								 |
								 END_NODE_REVERSE
								 | X_NODE) ?
			    "ellipse" : "point";

			char *color = db_node_check_for_any_flag(node,
								 STARTING_FORWARD
								 )
			    ? "turquoise4" : db_node_check_for_any_flag(node,
									X_NODE)
			    ? "red" : db_node_check_for_any_flag(node,
								 BRANCH_NODE_FORWARD
								 |
								 BRANCH_NODE_REVERSE)
			    ? "orange" : db_node_check_for_any_flag(node,
								    END_NODE_FORWARD
								    |
								    END_NODE_REVERSE)
			    ? "royalblue" : "black";
			printf("%s [label=%s, shape=%s, color=%s]\n", seq1,
			       seq1, print, color);
			binary_kmer_to_seq(&tmp, kmer_size, seq);
			binary_kmer_left_shift(&tmp, 2, kmer_size);
			Orientation o = forward;
			Edges e = db_node_get_edges_for_orientation(node, o);
			Edges e2 =
			    db_node_get_edges_for_orientation(node, reverse);
			/*
			 */
			void hasEdge1(Nucleotide n) {
				Edges et = e;
				BinaryKmer bk;
				Key k = &bk;
				if (((et >> n) & 1) == 1) {
					binary_kmer_modify_base(&tmp, n,kmer_size, 0);
					binary_kmer_to_seq(element_get_key(&tmp, kmer_size, k), kmer_size, seqNext);
					printf("%s -> %s [ label = \"%c\", color=\"%s\"];\n",
					       seq, seqNext,
					       binary_nucleotide_to_char(n),
					       "blue");
				}
			}
			nucleotide_iterator(&hasEdge1);
			binary_kmer_assignment_operator(tmp, co);
			binary_kmer_left_shift(&tmp, 2, kmer_size);
			o = reverse;
			e = e2;
			void hasEdge2(Nucleotide n) {
				Edges et = e;
				BinaryKmer bk;
				Key k = &bk;
				if (((et >> n) & 1) == 1) {
					binary_kmer_modify_base(&tmp, n,
								kmer_size, 0);
					binary_kmer_to_seq(element_get_key
							   (&tmp, kmer_size, k),
							   kmer_size, seqNext);
					printf
					    ("%s -> %s [ label = \"%c\", color=\"%s\"];\n",
					     seq, seqNext,
					     binary_nucleotide_to_char(n),
					     "red");
				}
			}
			nucleotide_iterator(&hasEdge2);
			
		}
	}
	hash_table_traverse(&print_graphviz, db_graph);
}

void write_graphviz_file(char *filename, dBGraph * db_graph)
{
	db_graph_write_graphviz_file(filename, db_graph);
}

void timestamp() {
	time_t ltime = time(NULL);
#ifdef ENABLE_BUBBLEPARSE
	char *timestring = asctime(localtime(&ltime));
	char  timestringcopy[64];
	strcpy(timestringcopy, timestring);
	timestringcopy[strlen(timestring)-1] = 0;
	log_and_screen_printf("\n----- %s -----\n", timestringcopy);
#else
	log_and_screen_printf("\n-----\n%s",asctime(localtime(&ltime)));	
#endif
	fflush(stdout);
}
