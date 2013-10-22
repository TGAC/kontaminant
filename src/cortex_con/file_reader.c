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
#include <err.h>
#include <assert.h>
#include <unistd.h>

#ifdef THREADS
#include <pthread.h>
#endif
#include <nucleotide.h>
#include <seq.h>
#include <binary_kmer.h>
#include <flags.h>
#include <element.h>
#include <dB_graph.h>
#include <file_format.h>
#include <file_reader.h>
#include <element.h>
#include <string.h>
#include <err.h>
#include <logger.h>


static void db_graph_load_sliding_windows(Element *previous_node,dBGraph * db_graph,short colour,boolean prev_full_entry,DBGraphFileReaderArgs *fra,short kmer_size, KmerSlidingWindowSet *windows) {
    Element *current_node = NULL;
    BinaryKmer tmp_kmer;
    Orientation current_orientation = forward;
    Orientation previous_orientation = reverse;
    int i, j;
    for (i = 0; i < windows->nwindows; i++) {	//for each window
        KmerSlidingWindow *current_window = &(windows->window[i]);
        for (j = 0; j < current_window->nkmers; j++) {	//for each kmer in window
            boolean found;
            Key key = element_get_key(&(current_window->kmer[j]),kmer_size, &tmp_kmer);
            if (fra->insert) {
                current_node = hash_table_find_or_insert(key, &found, db_graph);
            }else{
                current_node =	hash_table_find(key, db_graph);
            }
            if(current_node != NULL){
                current_orientation = db_node_get_orientation(&(current_window->kmer[j]),current_node, kmer_size);
                if (!(i == 0 && j == 0 && prev_full_entry == false && current_node == previous_node)) {	//otherwise is the same old last entry
                    element_update_coverage(current_node, colour, 1);
                }
            }
            
            if (j > 0) {
                if(previous_node != NULL){//RHRG. This was asserting before, now we validate so we can ignore the first node in case of solid. We could make a more robust code but I don't think is worth a complex validation here.
                    //fprintf(stderr, "Adding edges between nodes: \n");
                    //db_node_print_to_stream(previous_node, db_graph, stderr);
                    //db_node_print_to_stream(current_node, db_graph, stderr);
                    
                    db_node_add_edge(previous_node, current_node, previous_orientation, current_orientation, kmer_size, colour);
                    
                    //fprintf(stderr, "Now: \n");
                    //db_node_print_to_stream(previous_node, db_graph, stderr);
                    //db_node_print_to_stream(current_node, db_graph, stderr);
                    
                    
                }else{
                    fprintf(stderr, "Unable to connect nodes!\n");
                }
            }
			previous_node = current_node;
			previous_orientation = current_orientation;
            
        }
    }
    
}

static int file_reader(FILE * fp, Sequence * seq, int max_read_length,
                       boolean new_entry, boolean * full_entry) {
    *full_entry = true;
    
    if (new_entry != true) {
        puts("new_entry has to be true for fastq\n");
        exit(1);
    }
    
    return read_sequence_from_fastq(fp, seq, max_read_length);
}

static int file_reader_wrapper(DBGraphFileReaderInnerArgs * fria){
    int length = 0;
    int offset;
    assert(fria->fp != NULL);
    switch(fria->format){
        case FASTA:
            
    		offset = 0;
    		if (fria->new_entry == false) {
                shift_last_kmer_to_start_of_sequence(fria->seq, fria->seq->length, fria->kmer_size);
    			offset = fria->kmer_size;
    		}
    		length =  read_sequence_from_fasta(fria->fp, fria->seq, fria->max_read_length, fria->new_entry, &fria->full_entry, offset);
            if (!fria->full_entry) {
                fria->new_entry = false;
            }
            break;
        case FASTQ:
    		length = read_sequence_from_fastq(fria->fp, fria->seq, fria->max_read_length);
            fria->full_entry = true;
            break;
        default:
            fprintf(stderr, "Format not supported yet:%s\n", file_format_to_string(fria->format));
		    exit(1);
            break;
    }
    if (fria->full_entry) {
        fria->new_entry = true;
    }else{
        fria->new_entry = false;
    }
    return length;
}


static long long load_seq_into_db_graph(DBGraphFileReaderArgs * fra){
    
    assert(fra != NULL);
    assert(fra->inner_args != NULL);
    assert(fra->db_graph != NULL);
    if (strcmp (fra->filename, "-") == 0 ) {
        fra->inner_args->fp = stdin;
    }else{
        fra->inner_args->fp =  fopen(fra->filename, "r");
    }
    if (fra->inner_args->fp == NULL) {
		fprintf(stderr, "cannot open file:%s\n", fra->filename);
        assert(access(fra->filename, R_OK) == 0);
		exit(1);	//TODO - prefer to print warning and skip file and return an error code?
	}
    
    long long *bad_reads = &fra->bad_reads;
    int fastq_ascii_offset = fra->fastq_ascii_offset;
    char quality_cut_off = fra->quality_cut_off;
    int max_read_length = fra->max_read_length;
    short colour = fra->colour;
    dBGraph * db_graph = fra->db_graph;
    
	//----------------------------------
	// preallocate the memory used to read the sequences
	//----------------------------------
	Sequence *seq = malloc(sizeof(Sequence));
	if (seq == NULL) {
		fputs("Out of memory trying to allocate Sequence\n", stderr);
		exit(1);
	}
    fra->inner_args->seq = seq;
	alloc_sequence(seq, max_read_length, LINE_MAX, fastq_ascii_offset);
    
	long long seq_length = 0;
	long long seq_count = 0;
    
	short kmer_size = db_graph->kmer_size;
    
    //int max_windows;
    //int max_kmers;
    KmerSlidingWindowSet * windows = binary_kmer_sliding_window_set_new_from_read_length(kmer_size,max_read_length);
    
    
	int entry_length;
    
	boolean prev_full_entry = true;
    
	Element *previous_node = NULL;
    boolean keep_reading = true;
    
    // print_file_reader_args(fra);
	while ((entry_length = file_reader_wrapper(fra->inner_args)) && keep_reading)
	{
		seq_count++;
		
		seq_length += (long long)(entry_length -  (prev_full_entry == false ? kmer_size : 0));
		int nkmers = get_sliding_windows_from_sequence(seq->seq, seq->qual,
                                                       entry_length,
                                                       quality_cut_off,
                                                       kmer_size,
                                                       windows, windows->max_nwindows,
                                                       windows->max_kmers,false,0);
		if (nkmers == 0) {
			(*bad_reads)++;
		}else{
            db_graph_load_sliding_windows(previous_node,db_graph,colour,prev_full_entry,fra,kmer_size, windows);
        }
        if (fra->inner_args->full_entry == false) {
            shift_last_kmer_to_start_of_sequence(seq, entry_length, kmer_size);
        }else{
            hash_table_add_number_of_reads(1, colour, db_graph);
        }
        prev_full_entry = fra->inner_args->full_entry ;
        if (prev_full_entry) {
            fra->inner_args->new_entry = true;
        }else{
            fra->inner_args->new_entry = false;
        }
        if (seq_count % 10000) {
            if (hash_table_percentage_occupied(db_graph) > fra->maximum_ocupancy) {
                keep_reading = false;
            }
        }
    }
    
    free_sequence(&seq);
    fra->inner_args->seq = NULL;
    binary_kmer_free_kmers_set(&windows);
    return seq_length;
}


long long load_fastq_from_filename_into_graph(char *filename, short colour,
                                              long long *bad_reads,
                                              char quality_cut_off,
                                              int max_read_length, int fastq_ascii_offset, dBGraph * db_graph) {
    
    
    log_and_screen_printf("\nLoading reads from: %s\n", filename);
    long long loaded_kmers = 0;
    DBGraphFileReaderArgs fra;
    fra.bad_reads = 0;
    fra.colour = colour;
    fra.fastq_ascii_offset = fastq_ascii_offset;
    DBGraphFileReaderInnerArgs fria;
    fria.format = FASTQ;
    fria.kmer_size = db_graph->kmer_size;
    fria.max_read_length = 1000;
    fria.new_entry = true;
    fra.filename = filename;
    fra.quality_cut_off = quality_cut_off;
    fra.inner_args = &fria;
    fra.insert = true;
    fra.max_read_length = 1000;
    fra.maximum_ocupancy = 100;
    fra.db_graph = db_graph;
    //fp, seq, max_read_length, full_entry, &full_entry
    
    loaded_kmers = load_seq_into_db_graph(&fra);
    printf("Loaded %'lld kmers (bad reads %'lld)\n", loaded_kmers, fra.bad_reads);
    hash_table_print_stats(db_graph);
    return loaded_kmers;
}





#if 0

long long load_fastq_from_filename_into_graph(char *filename, short colour,
                                              long long *bad_reads,
                                              char quality_cut_off,
                                              int max_read_length, int fastq_ascii_offset, dBGraph * db_graph)

{
    
    FILE * fp = NULL;
	
    
    
	
    if (strcmp("-", filename) == 0) {
        fp = stdin;
    }else{
        fp = fopen(filename, "r");
        
    }
	if (fp == NULL) {
		fprintf(stderr, "cannot open file:%s\n", filename);
		exit(1);	//TODO - prefer to print warning and skip file and return an error code?
	}
	
    
	long long ret = load_seq_into_graph(fp, &file_reader, bad_reads, fastq_ascii_offset, quality_cut_off, max_read_length, colour, db_graph);
    
    if (fp != stdin) {
        fclose(fp);
    }
	return ret;
}

#endif
long long load_fastq_cov_from_filename_into_graph(char *filename, short colour, long long *bad_reads, char quality_cut_off, int max_read_length, int fastq_ascii_offset, dBGraph * db_graph)
{
    
	//int file_reader(FILE * fp, Sequence * seq, int max_read_length,
    //                boolean new_entry, boolean * full_entry) {
	//	*full_entry = true;
    
	//	if (new_entry != true) {
	//		puts("new_entry has to be true for fastq\n");
	//		exit(1);
	//	}
    
	//	return read_sequence_from_fastq(fp, seq, max_read_length);
	//}
    
    FILE * fp = NULL;
    if (strcmp("-", filename) == 0) {
        fp = stdin;
    }else{
        fp = fopen(filename, "r");
        
    }
    if (fp == NULL) {
        fprintf(stderr, "cannot open file:%s\n", filename);
        exit(1);	//TODO - prefer to print warning and skip file and return an error code?
    }
    
    
	long long ret = load_seq_cov_into_graph(fp, &file_reader, bad_reads, fastq_ascii_offset, quality_cut_off, max_read_length, colour, db_graph);
    
    if (fp != stdin) {
        fclose(fp);
    }
	
	return ret;
}


long long load_fasta_from_filename_into_graph(char *filename, short colour,
                                              long long *bad_reads,
                                              int max_chunk_length, dBGraph * db_graph){
    
    log_and_screen_printf("\nLoading reads from: %s\n", filename);
    long long loaded_kmers = 0;
    DBGraphFileReaderArgs fra;
    fra.bad_reads = 0;
    fra.colour = colour;
    fra.fastq_ascii_offset = 0;
    DBGraphFileReaderInnerArgs fria;
    fria.format = FASTA;
    fria.kmer_size = db_graph->kmer_size;
    fria.max_read_length = 1000;
    fria.new_entry = true;
    fria.full_entry = false;
    fra.filename = filename;
    fra.quality_cut_off = 0;
    fra.inner_args = &fria;
    fra.insert = true;
    fra.max_read_length = 1000;
    fra.maximum_ocupancy = 95;
    fra.db_graph = db_graph;
    //fp, seq, max_read_length, full_entry, &full_entry
    
    loaded_kmers = load_seq_into_db_graph(&fra);
    printf("Loaded %'lld kmers (bad reads %'lld)\n", loaded_kmers, fra.bad_reads);
    hash_table_print_stats(db_graph);
    return loaded_kmers;
    
}

//this routine supports big fasta entries (chromosome length for example)
#if 0
long long load_fasta_from_filename_into_graph(char *filename, short colour,
                                              long long *bad_reads,
                                              int max_chunk_length, dBGraph * db_graph){
    printf("MAX CHUNK LENGTH: %d\n", max_chunk_length);
	int file_reader(FILE * fp, Sequence * seq, int max_read_length,
                    boolean new_entry, boolean * full_entry) {
        printf("MAX CHUNK LENGTH INSIDE: %d", max_chunk_length);
		long long r;
		int offset = 0;
		if ( new_entry == false ) {
			offset = db_graph->kmer_size;
		}
		r = read_sequence_from_fasta(fp, seq, max_read_length,  new_entry, full_entry, offset);
        
		return r;
	}
    
	FILE * fp = NULL;
    if (strcmp("-", filename) == 0) {
        fp = stdin;
    }else{
        fp = fopen(filename, "r");
    }
	if (fp == NULL) {
		fprintf(stderr, "cannot open file:%s\n", filename);
		exit(1);	//TODO - prefer to print warning and skip file and return an error code?
	}
    
    
    long long ret = 0;
    
    ret = load_seq_into_graph(fp, &file_reader, bad_reads, 0, 0, max_chunk_length, colour, db_graph);
    
    
    if (fp != stdin) {
       	fclose(fp);
    }
    
    
	return ret;
}
#endif
//returns length of sequence loaded

long long load_seq_into_graph(FILE * fp,
                              int (*file_reader) (FILE * fp, Sequence * seq,
                                                  int max_read_length,
                                                  boolean new_entry,
                                                  boolean * full_entry),
                              long long *bad_reads, int fastq_ascii_offset, char quality_cut_off,
                              int max_read_length, short colour, dBGraph * db_graph)
{
	//----------------------------------
	// preallocate the memory used to read the sequences
	//----------------------------------
	Sequence *seq = malloc(sizeof(Sequence));
	if (seq == NULL) {
		fputs("Out of memory trying to allocate Sequence\n", stderr);
		exit(1);
	}
    
	alloc_sequence(seq, max_read_length, LINE_MAX, fastq_ascii_offset);
    
	long long seq_length = 0;
	long long seq_count = 0;
    
	short kmer_size = db_graph->kmer_size;
    
	//max_read_length/(kmer_size+1) is the worst case for the number of sliding windows, ie a kmer follow by a low-quality/bad base
	int max_windows = max_read_length / (kmer_size + 1);
    
	//number of possible kmers in a 'perfect' read
	int max_kmers = max_read_length - kmer_size + 1;
    
	//----------------------------------
	//preallocate the space of memory used to keep the sliding_windows. NB: this space of memory is reused for every call -- with the view
	//to avoid memory fragmentation
	//NB: this space needs to preallocate memory for orthogonal situations:
	//    * a good read -> few windows, many kmers per window
	//    * a bad read  -> many windows, few kmers per window
	//----------------------------------
	KmerSlidingWindowSet *windows = malloc(sizeof(KmerSlidingWindowSet));
	if (windows == NULL) {
		fputs("Out of memory trying to allocate a KmerArraySet",
		      stderr);
		exit(1);
	}
	//allocate memory for the sliding windows
	binary_kmer_alloc_kmers_set(windows, max_windows, max_kmers);
    
	int entry_length;
    
	boolean full_entry = true;
	boolean prev_full_entry = true;
    
	Element *previous_node = NULL;
    
    //RHRG: Variables to deal with solid.
#ifdef SOLID
    NucleotideBaseSpace first_base = Undef;
    NucleotideBaseSpace first_real_base = Undef;
    Nucleotide first_transition = Undefined;
#endif
    
	while ((entry_length =
            file_reader(fp, seq, max_read_length, full_entry, &full_entry)))
	{
		seq_count++;
        
		//printf("entry %s %i %i %s\n",seq->name,seq->start,seq->end,prev_full_entry?"not connect":"connect");
        
		if (DEBUG) {
			printf
            ("\nsequence %s - kmer size: %i - entry length: %i - max kmers:%i\n",
             seq->seq, db_graph->kmer_size, entry_length,
             max_kmers);
		}
        
		int i, j;
		seq_length +=
        (long long)(entry_length -
                    (prev_full_entry ==
                     false ? db_graph->kmer_size : 0));
        
		//printf("Length %qd %i %s %i %i\n",seq_length,entry_length, seq->name, seq->start, seq->end);
        
		int nkmers = get_sliding_windows_from_sequence(seq->seq, seq->qual,
                                                       entry_length,
                                                       quality_cut_off,
                                                       db_graph->kmer_size,
                                                       windows, max_windows,
                                                       max_kmers,false,0);
        
		//printf("number kmers:%i %i\n",nkmers,windows->nwindows);
        
		if (nkmers == 0) {
			(*bad_reads)++;
		} else {
			Element *current_node = NULL;
			BinaryKmer tmp_kmer;
            
			Orientation current_orientation = undefined;
			Orientation previous_orientation = forward;
            
			for (i = 0; i < windows->nwindows; i++) {	//for each window
				KmerSlidingWindow *current_window = &(windows->window[i]);
                
				for (j = 0; j < current_window->nkmers; j++) {	//for each kmer in window
					boolean found = false;
                    boolean search = true;
                    //Dont try to a
#ifdef SOLID
                    if ((i == 0 && j == 0 && prev_full_entry == true )){//Only happens on the first chunk of the fasta file. We don't expect them to be too big, anyway. Just to be consistent and in the case we have  genome in color space
                        search = false;
                        first_base = char_to_binary_nucleotide_base_space(sequence_get_base(0, seq));
                    }
                    
#endif
                    if(search){//This only makes sense when working with solid, this way we don't add tho the graph the first base.
                        current_node =	hash_table_find_or_insert(element_get_key(&(current_window->kmer[j]),db_graph->kmer_size, &tmp_kmer), &found, db_graph);
                        long long int k_read = hash_table_get_unique_kmers(db_graph);
                        if (k_read % 1000000 == 0){
                            printf("  Read %'lld kmers\n", k_read);
                            fflush(stdout);
                        }
                        
                        assert(current_node != NULL);
                        
                        if (!(i == 0 && j == 0 && prev_full_entry == false && current_node == previous_node)) {	//otherwise is the same old last entry
                            element_update_coverage(current_node, colour, 1);
                        }
                        current_orientation = db_node_get_orientation(&(current_window->kmer[j]),current_node,db_graph->kmer_size);
#ifdef SOLID
                        //RHRG: Here we add the first base to the first
                        if(first_base != Undef && i == 0 && j == 1 ){
                            first_transition = char_to_binary_nucleotide(sequence_get_base(1, seq));
                            first_real_base = binary_nucleotide_base_space_get_next_base(first_base, first_transition);
                            db_node_add_starting_base(first_real_base,current_orientation ,current_node);
                            first_real_base = Undef;
                            first_base = Undef;
                        }
#endif
                    }
                    
                    
                    
                    
                    
					
#ifdef __DEBUG
					if (DEBUG) {
						char kmer_seq[db_graph->kmer_size];
						char kmer_seq2[db_graph-> kmer_size];
                        
						printf
                        ("kmer i:%i j:%i:  %s %s %i\n",
                         i, j,
                         binary_kmer_to_seq(&
                                            (current_window->kmer
                                             [j]),
                                            db_graph->kmer_size,
                                            kmer_seq),
                         binary_kmer_to_seq
                         (binary_kmer_reverse_complement
                          (&
                           (current_window->kmer
							[j]),
                           db_graph->kmer_size,
                           &tmp_kmer),
                          db_graph->kmer_size,
                          kmer_seq2),
                         element_get_coverage_all_colours
                         (current_node));
					}
#endif
                    
					if (j > 0) {
                        if(previous_node != NULL){//RHRG. This was asserting before, now we validate so we can ignore the first node in case of solid. We could make a more robust codem but I don't think is worth a complex validation here.
                            db_node_add_edge(previous_node, current_node, previous_orientation, current_orientation, db_graph->kmer_size, colour);
                        }
                    }
					previous_node = current_node;
					previous_orientation = current_orientation;
                    
				}
			}
            
		}
        
		if (full_entry == false) {
			shift_last_kmer_to_start_of_sequence(seq, entry_length,
                                                 db_graph->kmer_size);
		}
        
		prev_full_entry = full_entry;
        
	}
    
	free_sequence(&seq);
	binary_kmer_free_kmers_set(&windows);
	return seq_length;
}

//returns the number of sequences
long long load_seq_cov_into_graph(FILE * fp,
                                  int (*file_reader) (FILE * fp, Sequence * seq,
                                                      int max_read_length,
                                                      boolean new_entry,
                                                      boolean * full_entry),
                                  long long *bad_reads, int fastq_ascii_offset, char quality_cut_off,
                                  int max_read_length, short colour, dBGraph * db_graph)

{
	//----------------------------------
	// preallocate the memory used to read the sequences
	//----------------------------------
	Sequence *seq = malloc(sizeof(Sequence));
	if (seq == NULL) {
		fputs("Out of memory trying to allocate Sequence\n", stderr);
		exit(1);
	}
    
	alloc_sequence(seq, max_read_length, LINE_MAX, fastq_ascii_offset);
    
	long long seq_length = 0;
	long long seq_count = 0;
    
	short kmer_size = db_graph->kmer_size;
    
	//max_read_length/(kmer_size+1) is the worst case for the number of sliding windows, ie a kmer follow by a low-quality/bad base
	int max_windows = max_read_length / (kmer_size + 1);
    
	//number of possible kmers in a 'perfect' read
	int max_kmers = max_read_length - kmer_size + 1;
    
	//----------------------------------
	//preallocate the space of memory used to keep the sliding_windows. NB: this space of memory is reused for every call -- with the view
	//to avoid memory fragmentation
	//NB: this space needs to preallocate memory for orthogonal situations:
	//    * a good read -> few windows, many kmers per window
	//    * a bad read  -> many windows, few kmers per window
	//----------------------------------
	KmerSlidingWindowSet *windows = malloc(sizeof(KmerSlidingWindowSet));
	if (windows == NULL) {
		fputs("Out of memory trying to allocate a KmerArraySet",
		      stderr);
		exit(1);
	}
	//allocate memory for the sliding windows
	binary_kmer_alloc_kmers_set(windows, max_windows, max_kmers);
    
	int entry_length;
    
	boolean full_entry = true;
	boolean prev_full_entry = true;
    
	Element *previous_node = NULL;
    
    //RHRG: Variables to deal with solid.
#ifdef SOLID
    NucleotideBaseSpace first_base = Undef;
    NucleotideBaseSpace first_real_base = Undef;
    Nucleotide first_transition = Undefined;
#endif
    
	while ((entry_length =file_reader(fp, seq, max_read_length, full_entry, &full_entry)))
	{
		seq_count++;
        
		//printf("entry %s %i %i %s\n",seq->name,seq->start,seq->end,prev_full_entry?"not connect":"connect");
        
		if (DEBUG) {
			printf
            ("\nsequence %s - kmer size: %i - entry length: %i - max kmers:%i\n",
             seq->seq, db_graph->kmer_size, entry_length,
             max_kmers);
		}
        
		int i, j;
		seq_length +=
        (long long)(entry_length -
                    (prev_full_entry ==
                     false ? db_graph->kmer_size : 0));
        
		//printf("Length %qd %i %s %i %i\n",seq_length,entry_length, seq->name, seq->start, seq->end);
        
		int nkmers = get_sliding_windows_from_sequence(seq->seq, seq->qual,
                                                       entry_length,
                                                       quality_cut_off,
                                                       db_graph->kmer_size,
                                                       windows, max_windows,
                                                       max_kmers,false,0);
        
		//printf("number kmers:%i %i\n",nkmers,windows->nwindows);
        
		if (nkmers == 0) {
			(*bad_reads)++;
		} else {
			Element *current_node = NULL;
			BinaryKmer tmp_kmer;
#ifdef KMER_TOOLS
            hash_table_add_number_of_reads(1, colour, db_graph);
#endif
            //			Orientation current_orientation = undefined;
            //			Orientation previous_orientation = forward;
            
			for (i = 0; i < windows->nwindows; i++) {	//for each window
				KmerSlidingWindow *current_window = &(windows->window[i]);
                
				for (j = 0; j < current_window->nkmers; j++) {	//for each kmer in window
					
                    boolean search = true;
                    //Dont try to a
#ifdef SOLID
                    if ((i == 0 && j == 0 && prev_full_entry == true )){//Only happens on the first chunk of the fasta file. We don't expect them to be too big, anyway. Just to be consistent and in the case we have  genome in color space
                        search = false;
                        first_base = char_to_binary_nucleotide_base_space(sequence_get_base(0, seq));
                    }
                    
#endif
                    if(search){//This only makes sense when working with solid, this way we don't add tho the graph the first base.
                        current_node =	hash_table_find(element_get_key(&(current_window->kmer[j]),db_graph->kmer_size, &tmp_kmer), db_graph);
                        long long int k_read = hash_table_get_unique_kmers(db_graph);
                        if (k_read % 1000000 == 0){
                            printf("  Read %'lld kmers\n", k_read);
                            fflush(stdout);
                        }
                        
                        if(current_node != NULL){
                            if (!(i == 0 && j == 0 && prev_full_entry == false && current_node == previous_node)) {	//otherwise is the same old last entry
                                element_update_coverage(current_node, colour, 1);
                            }
                            
#ifdef SOLID
                            //RHRG: Here we add the first base to the first
                            if(first_base != Undef && i == 0 && j == 1 ){
                                first_transition = char_to_binary_nucleotide(sequence_get_base(1, seq));
                                first_real_base = binary_nucleotide_base_space_get_next_base(first_base, first_transition);
                                db_node_add_starting_base(first_real_base,current_orientation ,current_node);
                                first_real_base = Undef;
                                first_base = Undef;
                            }
#endif
                        }
                        
                    }
                    
					previous_node = current_node;
                    //	previous_orientation = current_orientation;
                    
				}
			}
            
		}
        
		if (full_entry == false) {
			shift_last_kmer_to_start_of_sequence(seq, entry_length,
                                                 db_graph->kmer_size);
		}
        
		prev_full_entry = full_entry;
        
	}
    
	free_sequence(&seq);
	binary_kmer_free_kmers_set(&windows);
	return seq_count;
}

//returns number of sequence loaded (ie all the kmers concatenated)
//count_kmers returns the number of new kmers


long long load_binary_from_filename_into_graph(char *filename, dBGraph * db_graph, short colour, boolean all_entries_are_unique)

{
    
	FILE *fp_bin = fopen(filename, "r");
    
	if (fp_bin == NULL) {
		printf("cannot open file:%s\n", filename);
		exit(1);	//TODO - prefer to print warning and skip file and return an error code?
	}
    
	uint32_t num_colours_in_binary=1;
	uint32_t mean_read_len;
	uint64_t total_seq;
	if ( !check_binary_signature(fp_bin, db_graph->kmer_size, BINVERSION, &num_colours_in_binary, &mean_read_len, &total_seq)){
        errx(1,"binary version or kmer_size are inconsistent");
	}
    
	dBNode node_from_file;
	boolean found;
	long long count = 0;
	BinaryKmer tmp_kmer;
    
	//Go through all the entries in the binary file
	while (db_node_read_binary
	       (fp_bin, db_graph->kmer_size, &node_from_file)) {
		count++;
        
		dBNode *current_node = NULL;
		if (!all_entries_are_unique) {
			current_node =  hash_table_find_or_insert(element_get_key (element_get_kmer(&node_from_file), db_graph->kmer_size, &tmp_kmer), &found, db_graph);
		} else {
			current_node = hash_table_insert(element_get_key(element_get_kmer (&node_from_file), db_graph->kmer_size, &tmp_kmer), db_graph);
		}
        
		db_node_set_edges(current_node, colour,
                          db_node_get_edges(&node_from_file));
        
		element_update_coverage(current_node, colour,
                                element_get_coverage_all_colours
                                (&node_from_file));
	}
    
	fclose(fp_bin);
    
	return (long long)db_graph->kmer_size * count;
    
}

// This is just like get_sliding_windows_from_sequence is seq.c, but this one breaks the window also when a kmer is not in the graph.
//The first argument - seq - is a C string in A,C,G,T format
//The second argument - quality - is a string of qualities for the sequence, one byte per base.
//quality cutoff argument defines the threshold for quality
//return total number of kmers read
//The third argument - length - is the length in bases of the sequence.
//return total number of kmers read
int
get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph
(char *seq, char *qualities, int length, char quality_cut_off,
 KmerSlidingWindowSet * windows, int max_windows, int max_kmers,
 dBGraph * db_graph) {
	short kmer_size = db_graph->kmer_size;
    
	char first_kmer[kmer_size + 1];
	first_kmer[kmer_size] = '\0';
    
	BinaryKmer tmp_bin_kmer;
	BinaryKmer tmp_bin_kmer2;
    
	int i = 0;		//current index
	int count_kmers = 0;
    
	if (seq == NULL) {
		fputs
        ("in get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph, seq is NULL\n",
         stderr);
		exit(1);
	}
    
	if (length < kmer_size || max_windows == 0 || max_kmers == 0) {
		return 0;
	}
    
	int index_windows = 0;
    
	//loop over the bases in the sequence
	//index i is the current position in input sequence -- it nevers decreases.
    
	do {
        
		//build first kmer, ie a stretch of kmer_size good qualities bases
		int j = 0;	//count how many good bases
        
		while ((i < length) && (j < kmer_size)) {
            
			//collects the bases in the first kmer
			first_kmer[j] = seq[i];
            
			if ((char_to_binary_nucleotide(seq[i]) == Undefined) ||
			    (quality_cut_off != 0
			     && qualities[i] <= quality_cut_off)) {
                    j = 0;	//restart the first kmer
                } else {
                    j++;
                }
            
			i++;
		}
        
		if (j == kmer_size) {	//ie we did not parse the entire sequence looking for a single good kmer, the first kmer
            
			if (hash_table_find
			    (element_get_key
			     (seq_to_binary_kmer
			      (first_kmer, kmer_size, &tmp_bin_kmer), kmer_size,
			      &tmp_bin_kmer2), db_graph) == NULL) {
                     j = 0;
                     i = i - (kmer_size - 1);	// first kmer may be bad because of the first base. Start again from just after that
                     continue;	//want to restart building a first kmer
                 }
            
			count_kmers++;
            
			//new sliding window
			if (index_windows >= max_windows) {
				fputs
                ("number of windows is bigger than max_windows",
                 stderr);
				exit(1);
			}
            
			KmerSlidingWindow *current_window =
            &(windows->window[index_windows]);
            
			int index_kmers = 0;
			//do first kmer
			seq_to_binary_kmer(first_kmer, kmer_size,
                               &tmp_bin_kmer);
			binary_kmer_assignment_operator(current_window->kmer
                                            [index_kmers],
                                            tmp_bin_kmer);
            
			//do the rest --
			index_kmers++;
            
			while (i < length) {
                
				if (index_kmers >= max_kmers) {
					fputs
                    ("number of kmers is bigger than max_kmers\n",
                     stderr);
					exit(1);
				}
                
				Nucleotide current_base =
                char_to_binary_nucleotide(seq[i]);
				BinaryKmer tmp_next_kmer;
                
				//set to previous kmer, then shift an add new base
				binary_kmer_assignment_operator(tmp_next_kmer,
                                                current_window->kmer
                                                [index_kmers -
                                                 1]);
				binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end
                (&tmp_next_kmer, current_base, kmer_size);
                
				if ((current_base == Undefined) || (quality_cut_off != 0
                                                    && qualities[i] <= quality_cut_off) ) {
                    break;
                } else if (hash_table_find (element_get_key (&tmp_next_kmer, kmer_size, &tmp_bin_kmer2), db_graph) == NULL) {
                    i++;
                    break;
                }
				binary_kmer_assignment_operator
                (current_window->kmer[index_kmers],
                 tmp_next_kmer);
				index_kmers++;
				count_kmers++;
				i++;
			}
            
			current_window->nkmers = index_kmers;
			index_windows++;
            
		}
	}
	while (i < length);
    
	windows->nwindows = index_windows;
    
	return count_kmers;
    
}

void read_fastq_and_print_reads_that_lie_in_graph(FILE * fp, FILE * fout,
                                                  int (*file_reader) (FILE * fp,
                                                                      Sequence *
                                                                      seq,
                                                                      int
                                                                      max_read_length,
                                                                      boolean
                                                                      new_entry,
                                                                      boolean *
                                                                      full_entry),
                                                  long long *bad_reads,
                                                  int max_read_length,
                                                  dBGraph * db_graph,
                                                  boolean is_for_testing, char
                                                  **for_test_array_of_clean_reads,
                                                  int *for_test_index)
{
	//----------------------------------
	// preallocate the memory used to read the sequences
	//----------------------------------
	Sequence *seq = malloc(sizeof(Sequence));
	if (seq == NULL) {
		fputs("Out of memory trying to allocate Sequence\n", stderr);
		exit(1);
	}
	alloc_sequence(seq, max_read_length, LINE_MAX, 0);
    
	int seq_length = 0;
	short kmer_size = db_graph->kmer_size;
    
	//max_read_length/(kmer_size+1) is the worst case for the number of sliding windows, ie a kmer follow by a low-quality/bad base
	int max_windows = max_read_length / (kmer_size + 1);
    
	//number of possible kmers in a 'perfect' read
	int max_kmers = max_read_length - kmer_size + 1;
    
	//----------------------------------
	//preallocate the space of memory used to keep the sliding_windows. NB: this space of memory is reused for every call -- with the view
	//to avoid memory fragmentation
	//NB: this space needs to preallocate memory for orthogonal situations:
	//    * a good read -> few windows, many kmers per window
	//    * a bad read  -> many windows, few kmers per window
	//----------------------------------
	KmerSlidingWindowSet *windows = malloc(sizeof(KmerSlidingWindowSet));
	if (windows == NULL) {
		fputs("Out of memory trying to allocate a KmerArraySet",
		      stderr);
		exit(1);
	}
	//allocate memory for the sliding windows
	binary_kmer_alloc_kmers_set(windows, max_windows, max_kmers);
    
	char tmpseq[db_graph->kmer_size + 1];
	tmpseq[db_graph->kmer_size] = '\0';
	tmpseq[0] = '\0';
    
	BinaryKmer tmp_kmer;
	boolean full_entry = true;
	boolean prev_full_entry = true;
    
	int entry_length;
    
	while ((entry_length =
            file_reader(fp, seq, max_read_length, full_entry, &full_entry)))
	{
        
		if (DEBUG) {
			printf("\nsequence %s\n", seq->seq);
		}
        
		int i, j;
		seq_length +=
        (long long)(entry_length -
                    (prev_full_entry ==
                     false ? db_graph->kmer_size : 0));
        
		//use quality cutoff of 0, arg 4 below
		// int nkmers = get_sliding_windows_from_sequence(seq->seq,seq->qual,entry_length,0,db_graph->kmer_size,windows,max_windows, max_kmers);
		int nkmers =
        get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph
        (seq->seq, seq->qual, entry_length, 0,
         windows, max_windows, max_kmers, db_graph);
        
		if (nkmers == 0) {
			(*bad_reads)++;
		} else {
			Element *current_node = NULL;
			//we will throw out any window that does not have ALL kmers in the graph
			for (i = 0; i < windows->nwindows; i++) {	//for each window
				KmerSlidingWindow *current_window =
                &(windows->window[i]);
                
				boolean all_kmers_in_this_window_are_in_graph =
                true;
                
				for (j = 0; j < current_window->nkmers; j++) {	//for each kmer in window
                    
					current_node =
                    hash_table_find(element_get_key
                                    (&
                                     (current_window->kmer
                                      [j]),
                                     db_graph->kmer_size,
                                     &tmp_kmer),
                                    db_graph);
					if (current_node == NULL) {
						all_kmers_in_this_window_are_in_graph
                        = false;
					}
				}
                
				if (all_kmers_in_this_window_are_in_graph ==
				    true) {
					//print out this window as a "read". If this read has many windows, we will print each as a separate read (provided they lie in the graph)
					if (is_for_testing == false) {
						fprintf(fout, "> read\n");
						fprintf(fout, "%s",
                                binary_kmer_to_seq(&
                                                   (current_window->kmer
                                                    [0]),
                                                   db_graph->kmer_size,
                                                   tmpseq));
						for (j = 1;
						     j < current_window->nkmers;
						     j++) {
							fprintf(fout, "%c",
                                    binary_nucleotide_to_char
                                    (binary_kmer_get_last_nucleotide
                                     (&
                                      (current_window->kmer
                                       [j]))));
						}
						fprintf(fout, "\n");
					} else {
						for_test_array_of_clean_reads
                        [*for_test_index][0] = '\0';
						strcat
                        (for_test_array_of_clean_reads
                         [*for_test_index],
                         binary_kmer_to_seq(&
                                            (current_window->kmer
                                             [0]),
                                            db_graph->kmer_size,
                                            tmpseq));
						for (j = 1;
						     j < current_window->nkmers;
						     j++) {
							char tmp_ch[2];
							tmp_ch[0] =
                            binary_nucleotide_to_char
                            (binary_kmer_get_last_nucleotide
                             (&
                              (current_window->kmer
                               [j])));
							tmp_ch[1] = '\0';
							strcat
                            (for_test_array_of_clean_reads
                             [*for_test_index],
                             tmp_ch);
						}
						*for_test_index =
                        *for_test_index + 1;
					}
				}
				//else
				//  {
				//  }
                
			}
		}
        
		if (full_entry == false) {
			shift_last_kmer_to_start_of_sequence(seq, entry_length,
                                                 db_graph->kmer_size);
		}
        
		prev_full_entry = full_entry;
        
	}
    
	free_sequence(&seq);
	binary_kmer_free_kmers_set(&windows);
    
}

void print_binary_signature(FILE * fp, uint32_t kmer_size, uint32_t num_cols, uint32_t mean_read_len, uint64_t total_seq)
{
    char magic_number[6];
    uint32_t version = BINVERSION;
    uint32_t num_bitfields = NUMBER_OF_BITFIELDS_IN_BINARY_KMER;
    
    magic_number[0]='C';
    magic_number[1]='O';
    magic_number[2]='R';
    magic_number[3]='T';
    magic_number[4]='E';
    magic_number[5]='X';
    
    fwrite(magic_number,sizeof(char),6,fp);
    fwrite(&version,sizeof(uint32_t),1,fp);
    fwrite(&kmer_size,sizeof(uint32_t),1,fp);
    fwrite(&num_bitfields, sizeof(uint32_t),1,fp);
    fwrite(&num_cols, sizeof(uint32_t), 1, fp);
    fwrite(&mean_read_len, sizeof(uint32_t), 1, fp);
    fwrite(&total_seq, sizeof(uint64_t), 1, fp);
    fwrite(magic_number,sizeof(char),6,fp);
}



//return yes if signature is consistent
//third argument is ignored for now.
boolean check_binary_signature(FILE * fp, uint32_t kmer_size, uint32_t binary_version, uint32_t* number_of_colours_in_binary, uint32_t* mean_read_len, uint64_t* total_seq)
{
    size_t read;
    char magic_number[6];
    boolean ret = false;
    
    read = fread(magic_number, sizeof(char),6,fp);
    
    if (read>0 &&
        magic_number[0]=='C' &&
        magic_number[1]=='O' &&
        magic_number[2]=='R' &&
        magic_number[3]=='T' &&
        magic_number[4]=='E' &&
        magic_number[5]=='X' )
    {
        uint32_t version;
        read = fread(&version,sizeof(uint32_t),1,fp);
        if (read>0 && version==BINVERSION)
        {
            uint32_t kmer_size2;
            read = fread(&kmer_size2,sizeof(uint32_t),1,fp);
            if ((read>0) && (kmer_size2 == kmer_size) )
            {
                uint32_t num_bitfields;
                read = fread(&num_bitfields,sizeof(uint32_t),1,fp);
                
                if ( (read>0) && (num_bitfields==NUMBER_OF_BITFIELDS_IN_BINARY_KMER) )
                {
                    
                    uint32_t num_cols;
                    read = fread(&num_cols,sizeof(uint32_t),1,fp);
                    
                    if ( (read>0) && (num_cols==1)  )
                    {
                        *number_of_colours_in_binary = num_cols;
                        
                        read = fread(mean_read_len,sizeof(uint32_t),1,fp);
                        if (read>0)
                        {
                            read = fread(total_seq,sizeof(uint64_t),1,fp);
                            if (read>0)
                            {
                                magic_number[0]='\0';
                                magic_number[1]='\0';
                                magic_number[2]='\0';
                                magic_number[3]='\0';
                                magic_number[4]='\0';
                                magic_number[5]='\0';
                                read = fread(magic_number,sizeof(char),6,fp);
                                if ((read>0) &&
                                    magic_number[0]=='C' &&
                                    magic_number[1]=='O' &&
                                    magic_number[2]=='R' &&
                                    magic_number[3]=='T' &&
                                    magic_number[4]=='E' &&
                                    magic_number[5]=='X' )
                                {
                                    ret = true;
                                }
                                else
                                {
                                    printf("Binary header is missing the final magic number; we read %s, mean read len %d and total seq %'qd\n", magic_number, *mean_read_len, (long long int)*total_seq);
                                }
                                
                            }
                            else
                            {
                                printf("Binary header does not contain total seq\n");
                            }
                        }
                        else
                        {
                            printf("Binary header does not contain a mean read length\n");
                        }
                        
                    }
                    else
                    {
                        printf("You are loading  binary with %d colours into a graph with %d colours - incompatible\n", num_cols, NUMBER_OF_COLOURS);
                    }
                }
                else
                {
                    printf("Kmer of binary matches the current graph. However this binary was dumped with a different max_kmer size to that of the current graph\n");
                    printf("This binary uses %d bitfields, and current graph uses %d\n", num_bitfields, NUMBER_OF_BITFIELDS_IN_BINARY_KMER);
                }
            }
            else
            {
                printf("You are loading a binary with kmer=%d into a graph with kmer=%d - incompatible\n", kmer_size2, kmer_size);
            }
        }
        else
        {
            printf("Binary versions do not match.\n");
        }
    }
    else
    {
        printf("Binary does not have magic number in header. Corrupt, or not a Cortex binary\n");
    }
    
    return ret;
}


long long seq_kmer_iterator(FILE * fp, void * args, void  (*kmer_action) (BinaryKmer * k, void * arg) ,int (*file_reader) (FILE * fp, Sequence * seq,
                                                                                                                           int max_read_length, boolean new_entry, boolean * full_entry),
                            long long *bad_reads, int fastq_ascii_offset, char quality_cut_off,
                            int max_read_length, short kmer_size)

{
	//----------------------------------
	// preallocate the memory used to read the sequences
	//----------------------------------
	Sequence *seq = malloc(sizeof(Sequence));
	if (seq == NULL) {
		fputs("Out of memory trying to allocate Sequence\n", stderr);
		exit(1);
	}
    
	alloc_sequence(seq, max_read_length, LINE_MAX, fastq_ascii_offset);
    
	long long seq_length = 0;
	long long seq_count = 0;
    
	//max_read_length/(kmer_size+1) is the worst case for the number of sliding windows, ie a kmer follow by a low-quality/bad base
	int max_windows = max_read_length / (kmer_size + 1);
    
	//number of possible kmers in a 'perfect' read
	int max_kmers = max_read_length - kmer_size + 1;
    
	//----------------------------------
	//preallocate the space of memory used to keep the sliding_windows. NB: this space of memory is reused for every call -- with the view
	//to avoid memory fragmentation
	//NB: this space needs to preallocate memory for orthogonal situations:
	//    * a good read -> few windows, many kmers per window
	//    * a bad read  -> many windows, few kmers per window
	//----------------------------------
	KmerSlidingWindowSet *windows = malloc(sizeof(KmerSlidingWindowSet));
	if (windows == NULL) {
		fputs("Out of memory trying to allocate a KmerArraySet",
		      stderr);
		exit(1);
	}
	//allocate memory for the sliding windows
	binary_kmer_alloc_kmers_set(windows, max_windows, max_kmers);
    
	int entry_length;
    
	boolean full_entry = true;
	boolean prev_full_entry = true;
    
	while ((entry_length =
            file_reader(fp, seq, max_read_length, full_entry, &full_entry)))
	{
		seq_count++;
        
		//printf("entry %s %i %i %s\n",seq->name,seq->start,seq->end,prev_full_entry?"not connect":"connect");
        
        
        
		int i, j;
		seq_length += (long long)(entry_length - (prev_full_entry == false ? kmer_size : 0));
        
		//printf("Length %qd %i %s %i %i\n",seq_length,entry_length, seq->name, seq->start, seq->end);
        
		int nkmers = get_sliding_windows_from_sequence(seq->seq, seq->qual,
                                                       entry_length,
                                                       quality_cut_off,
                                                       kmer_size,
                                                       windows, max_windows,
                                                       max_kmers,false,0);
        
		//printf("number kmers:%i %i\n",nkmers,windows->nwindows);
        
		if (nkmers == 0) {
			(*bad_reads)++;
		} else {
            
            
			for (i = 0; i < windows->nwindows; i++) {	//for each window
				KmerSlidingWindow *current_window = &(windows->window[i]);
                
				for (j = 0; j < current_window->nkmers; j++) {	//for each kmer in window
					kmer_action(&(current_window->kmer[j]), args );
                    
				}
			}
            
		}
        
		if (full_entry == false) {
			shift_last_kmer_to_start_of_sequence(seq, entry_length, kmer_size);
		}
        
		prev_full_entry = full_entry;
        
	}
    
	free_sequence(&seq);
	binary_kmer_free_kmers_set(&windows);
	return seq_length;
}

