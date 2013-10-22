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
 
#ifndef FILE_READER_H_
#define FILE_READER_H_

#include <nucleotide.h>
#include <dB_graph.h>
#include <seq.h>
#include <file_format.h>
typedef struct{
    FILE * fp;
    Sequence * seq;
    int max_read_length;
    boolean new_entry;
    boolean full_entry;
    FileFormat format;
    short kmer_size;
} DBGraphFileReaderInnerArgs;

typedef struct {
    char * filename; 
    short colour;
    long long bad_reads;
    char quality_cut_off;
    int max_read_length;
    int fastq_ascii_offset;
    DBGraphFileReaderInnerArgs * inner_args;//This holds thread-specific information. 
    float maximum_ocupancy;
    boolean stop_on_full;
    boolean insert;
    dBGraph * db_graph;
}DBGraphFileReaderArgs;


//int file_reader_wrapper(KmerFileReaderInnerArgs * fria);

//long long load_seq_into_db_graph(DBGraphFileReaderArgs * fra);



DBGraphFileReaderInnerArgs * open_db_graph_file_reader_wrapper(char * filename,  short kmer_size,int max_read_length,  FileFormat format, char offset);

long long load_kmers_binary_from_filename_update_coverage(char *filename, dBGraph * db_graph, short colour);

//void close_kmer_file_reader_wrapper(DBGraphFileReaderInnerArgs ** fria );
/**
 * these routines return the length of the read sequence, for the binary file is all the kmers conctenated
 * The first part is available when compiling for read pairs. 
 * The second block of definitions is mostly the same, just taking in account the read pairs. 
 * This if block system may be used when adding extra functionality. The point about this is to don't repeat
 * code and give some kind of "polymorphism" depending on the program which is being compiled. 
 */

//for short fasta entries - reads or similar 
long long load_fasta_from_filename_into_graph(char* filename, short colour, long long * bad_reads, int max_chunk_length, dBGraph* db_graph);

//for fastq
long long load_fastq_from_filename_into_graph(char* filename, short colour, long long * bad_reads,  char quality_cut_off, int max_read_length, int fastq_ascii_offset, dBGraph* db_graph);



//Actually loading the sequence
long long load_seq_into_graph(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry), long long * bad_reads, int fastq_ascii_offset, char quality_cut_off, int max_read_length, short colour, dBGraph * db_graph);

long long load_seq_cov_into_graph(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry), long long * bad_reads, int fastq_ascii_offset, char quality_cut_off, int max_read_length, short colour, dBGraph * db_graph);

long long load_fastq_cov_from_filename_into_graph(char *filename, short colour, long long *bad_reads, char quality_cut_off, int max_read_length, int fastq_ascii_offset, dBGraph * db_graph);

//for binary
long long load_binary_from_filename_into_graph(char *filename,
                                               dBGraph * db_graph, short colour,
                                               boolean all_entries_are_unique);



//for reference
void read_ref_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry),
                                                                            long long * bad_reads, char quality_cut_off, int max_read_length, dBGraph * db_graph);

void read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference(char* f_name, dBGraph* db_graph);

void read_all_ref_chromosomes_and_mark_graph(dBGraph* db_graph);

// This would be in seq.c, but that has no knowledge of dBGraph
int get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(char * seq,  char * qualities, int length, char quality_cut_off,
                                                                                  KmerSlidingWindowSet * windows, int max_windows, int max_kmers, dBGraph* db_graph);

// for dumping clean fasta files from fastq - ie only holding reads that lie entirely in a (presumably cleaned) graph
void read_fastq_and_print_reads_that_lie_in_graph(FILE* fp, FILE* fout, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry), 
						  long long * bad_reads, int max_read_length, dBGraph * db_graph,
						  boolean is_for_testing, char** for_test_array_of_clean_reads, int* for_test_index);

//This are used to validate that the graph is correct.
void print_binary_signature(FILE * fp, uint32_t kmer_size, uint32_t num_cols, uint32_t mean_read_len, uint64_t total_seq);
boolean check_binary_signature(FILE * fp, uint32_t kmer_size, uint32_t binary_version, uint32_t* number_of_colours_in_binary, uint32_t* mean_read_len, uint64_t* total_seq);
void binary_kmers_sliding_window_reset_iterator(KmerSlidingWindowSet * ksws);
//To make some analysis, iterator over the kmers of a sequence. 
long long seq_kmer_iterator(FILE * fp, void * args, void  (*kmer_action) (BinaryKmer * k, void * arg) ,int (*file_reader) (FILE * fp, Sequence * seq,
					int max_read_length, boolean new_entry, boolean * full_entry),
		    		long long *bad_reads, int fastq_ascii_offset, char quality_cut_off,
					int max_read_length, short kmer_size);

#endif /* FILE_READER_H_ */
