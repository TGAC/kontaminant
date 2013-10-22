//
//  kmer_reader.h
//  Untitled
//
//  Created by Ricardo Ramirez-Gonzalez (TGAC) on 12/12/2011.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef kmer_reader_h
#define kmer_reader_h


typedef struct{
    FILE * fp;
    Sequence * seq;
    int max_read_length;
    boolean new_entry;
    boolean full_entry;
    FileFormat format;
    short kmer_size;
} KmerFileReaderInnerArgs;

typedef struct {
    char * filename; 
    short colour;
    long long bad_reads;
    char quality_cut_off;
    int max_read_length;
    int fastq_ascii_offset;
    KmerFileReaderInnerArgs * inner_args;//This holds thread-specific information. 
    float maximum_ocupancy;
    boolean stop_on_full;
    boolean insert;
    KmerHash * KmerHash;
}KmerFileReaderArgs;


//int file_reader_wrapper(KmerFileReaderInnerArgs * fria);

long long load_seq_into_kmers_hash(KmerFileReaderArgs * fra);

void kmer_hash_dump_binary(char *filename, boolean(*condition) (dBNode * node), KmerHash * db_graph);
long long load_binary_from_filename_into_kmers_hash(char *filename, dBGraph * db_graph, short colour, boolean all_entries_are_unique);

long long load_kmers_binary_from_filename_update_coverage(char *filename, dBGraph * db_graph, short colour);

KmerFileReaderInnerArgs * open_kmer_file_reader_wrapper(char * filename,  short kmer_size,int max_read_length,  FileFormat format, char offset);

int file_reader_wrapper(KmerFileReaderInnerArgs * fria);
void close_kmer_file_reader_wrapper(KmerFileReaderInnerArgs ** fria );
#endif
