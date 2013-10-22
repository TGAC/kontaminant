//
//  kmer_reader.c
//  Untitled
//
//  Created by Ricardo Ramirez-Gonzalez (TGAC) on 12/12/2011.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <locale.h>
#include <getopt.h>
#include <err.h>
#include <assert.h>
#include <unistd.h>

#include <open_hash/hash_table.h>
#include <dB_graph.h>
#include <logger.h>
#include <file_format.h>
#include <nucleotide.h>
#include <seq.h>

#include <kmer_hash.h>
#include <kmer_reader.h>


//Returns the count of how many kmers had the coverage updated
static int kmer_hash_load_sliding_windows(Element *previous_node,KmerHash *kmer_hash,short colour,boolean prev_full_entry,KmerFileReaderArgs *fra,short kmer_size, KmerSlidingWindowSet *windows) {
    Element *current_node = NULL;
    BinaryKmer tmp_kmer;
    int i, j, added =0 ;
    for (i = 0; i < windows->nwindows; i++) {	//for each window
        KmerSlidingWindow *current_window = &(windows->window[i]);
        for (j = 0; j < current_window->nkmers; j++) {	//for each kmer in window
            boolean found;
            Key key = element_get_key(&(current_window->kmer[j]),kmer_size, &tmp_kmer);
            if (fra->insert) {
                current_node = hash_table_find_or_insert(key, &found, kmer_hash);
            }else{
                current_node =	hash_table_find(key, kmer_hash);
            }
            if(current_node != NULL){
                if (!(i == 0 && j == 0 && prev_full_entry == false && current_node == previous_node)) {	//otherwise is the same old last entry
                    element_update_coverage(current_node, colour, 1);
                    added++;
                    
                }
            }
            previous_node = current_node;
            
        }
    }
    return added;
}

//Returns null if there is no file to read. 
KmerFileReaderInnerArgs * open_kmer_file_reader_wrapper(char * filename,  short kmer_size,int max_read_length,  FileFormat format, char offset){
    if (filename == NULL) {
        return NULL;
    }
    KmerFileReaderInnerArgs * fria = calloc(1, sizeof(KmerFileReaderInnerArgs));
    
    fria->format = format;
    fria->full_entry = false;
    fria->new_entry = true;
    fria->kmer_size = kmer_size;
    fria->max_read_length = max_read_length;
    fria->seq = sequence_new(max_read_length, max_read_length, offset);
    fria->fp = fopen(filename, "r");
    
    if (fria->seq == NULL) {
        fprintf(stderr, "Unable to allocate Sequence.\n");
        exit(-1);
    }
    
    if (fria->fp == NULL) {
        fprintf(stderr, "Unable to open file %s\n", filename);
        exit(-1);
    }
    
    
    return fria;

}

void close_kmer_file_reader_wrapper(KmerFileReaderInnerArgs ** fria ){
    
    fclose(fria[0]->fp);
    free(fria[0]->seq);
    fria[0] = NULL;
}

int file_reader_wrapper(KmerFileReaderInnerArgs * fria){
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

static void print_file_reader_args(KmerFileReaderArgs * fra){
#ifdef DEBUG
    log_and_screen_printf("FileReaderArgs: %p\n", fra);
    log_and_screen_printf("Filename: %s\n", fra->filename );
    log_and_screen_printf("Colour: %d\n", fra->colour);
    log_and_screen_printf("bad_reads: %lld\n", fra->bad_reads);
    log_and_screen_printf("quality_cut_off: %d\n", fra->quality_cut_off);
    log_and_screen_printf("max_read_length: %d\n", fra->max_read_length);
    log_and_screen_printf("ascii offset: %d\n", fra->fastq_ascii_offset);
    log_and_screen_printf("Inner ARGS: %p\n", fra->inner_args);
    
    log_and_screen_printf("\tFP: %p\n", fra->inner_args->fp);
    log_and_screen_printf("\tseq: %p\n", fra->inner_args->seq);
    log_and_screen_printf("\tmax_read_length: %d\n", fra->inner_args->max_read_length);
    log_and_screen_printf("\tnew_entry: %d\n", fra->inner_args->new_entry);
    log_and_screen_printf("\tfull_entr: %d\n", fra->inner_args->full_entry);
    log_and_screen_printf("\tformat: %s\n", file_format_to_string(fra->inner_args->format));
    log_and_screen_printf("\tkmer_size: %d\n", fra->inner_args->kmer_size);
    
    log_and_screen_printf("Maximum ocupancy: %f\n", fra->maximum_ocupancy);
    log_and_screen_printf("Stop on full: %d\n", fra->stop_on_full);
    log_and_screen_printf("Insert: %d\n", fra->insert);
    log_and_screen_printf("KmerHash: %p\n", fra->KmerHash);
    hash_table_print_stats(fra->KmerHash);
#endif
    
}


long long load_seq_into_kmers_hash(KmerFileReaderArgs * fra){
    
    assert(fra != NULL);
    assert(fra->inner_args != NULL);
    assert(fra->KmerHash != NULL);
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
    KmerHash * kmer_hash = fra->KmerHash;
    
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
    long long contaminated_reads =0;
	short kmer_size = kmer_hash->kmer_size;
    
//int max_windows;
//int max_kmers;
   KmerSlidingWindowSet * windows = binary_kmer_sliding_window_set_new_from_read_length(kmer_size,max_read_length);

    
	int entry_length;
    int kmers_loaded;
    
	boolean prev_full_entry = true;
    
	Element *previous_node = NULL;
    boolean keep_reading = true;
    
    print_file_reader_args(fra);
	while ((entry_length = file_reader_wrapper(fra->inner_args)) && keep_reading)
	{
		
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
            kmers_loaded += kmer_hash_load_sliding_windows(previous_node,kmer_hash,colour,prev_full_entry,fra,kmer_size, windows);
        }
        if (fra->inner_args->full_entry == false) {
            shift_last_kmer_to_start_of_sequence(seq, entry_length, kmer_size);
        }else{
            hash_table_add_number_of_reads(1, colour, kmer_hash);

            
            kmer_hash->contaminated_kmers_per_read[kmers_loaded < MAX_READ_LENGTH? kmers_loaded:MAX_READ_LENGTH]++; //We allow up to the maximum read lenght. the last element is the cumulative of the reads/contigs that have more than the space we have allocated
            if(kmers_loaded > 0){
                contaminated_reads++;
            }

            seq_count++;
            kmers_loaded=0;
        }
        prev_full_entry = fra->inner_args->full_entry ;
        if (seq_count % 10000) {
            if (hash_table_percentage_occupied(kmer_hash) > fra->maximum_ocupancy) {
                keep_reading = false;
            }
        } 
    }
    
    if(colour == KMER_HASH_SAMPLE_INDEX){
        kmer_hash_add_contaminated_reads(contaminated_reads, kmer_hash);
    }
    
    free_sequence(&seq);
    fra->inner_args->seq = NULL;
    binary_kmer_free_kmers_set(&windows);
    return seq_length;
}

static void kmer_print_binary_signature(FILE * fp, uint32_t kmer_size, uint32_t num_cols, uint32_t mean_read_len, uint64_t total_seq)
{
    char magic_number[6];
    uint32_t version = BINVERSION;
    uint32_t num_bitfields = NUMBER_OF_BITFIELDS_IN_BINARY_KMER;
    
    magic_number[0]='K';
    magic_number[1]='M';
    magic_number[2]='E';
    magic_number[3]='R';
    magic_number[4]='S';
    magic_number[5]=' ';
    
    fwrite(magic_number,sizeof(char),6,fp);
    fwrite(&version,sizeof(uint32_t),1,fp);
    fwrite(&kmer_size,sizeof(uint32_t),1,fp);
    fwrite(&num_bitfields, sizeof(uint32_t),1,fp);
    fwrite(&num_cols, sizeof(uint32_t), 1, fp);
    fwrite(&mean_read_len, sizeof(uint32_t), 1, fp);
    fwrite(&total_seq, sizeof(uint64_t), 1, fp);
    fwrite(magic_number,sizeof(char),6,fp);
}


struct print_binary_args{
    boolean(*condition) (dBNode * node);
    long long count;
    FILE * fout;
    short kmer_size;
};


static void print_node_binary(dBNode * node, void * args) {
    struct print_binary_args * arguments = (struct print_binary_args * ) args;
    
    if (arguments->condition(node)) {
        arguments->count++;
        db_node_print_binary(arguments->fout, node, arguments->kmer_size);
    }
}

// We want to make this for the new elemetns., at the moment is a copy of the standard cortex output routine
void kmer_hash_dump_binary(char *filename, boolean(*condition) (dBNode * node), KmerHash * db_graph)
{
	FILE *fout;		//binary output
	fout = fopen(filename, "wb");
	if (fout == NULL) {
		fprintf(stderr, "cannot open %s", filename);
		exit(1);
	}
    
	int mean_read_len=0;//Mario - you can plumb this in if you want. See cortex_var/core/graph_info.c, and how the GraphInfo object is used in my file_reader * cprtex_var.c, 
	// for how I did it.
	long long total_seq = hash_table_get_number_of_reads(0, db_graph);
	kmer_print_binary_signature(fout,db_graph->kmer_size,1, mean_read_len, total_seq);
    
	//routine to dump graph
	
    struct print_binary_args args;
    args.condition = condition;
    args.count = 0;
    args.fout = fout;
    args.kmer_size = db_graph->kmer_size;
    
    void *  pass[1];
    pass[0] = (void*)&args;
    
//	hash_table_traverse(&print_node_binary, db_graph);
 
    
    hash_table_traverse_with_args(&print_node_binary, (void **) &pass, db_graph);
  //  void hash_table_traverse_with_args(void (*f)(Element *, void *),void ** args, HashTable * hash_table){
    
    
	fclose(fout);
    
	log_and_screen_printf("%'lld kmers dumped\n", args.count);
}

//return yes if signature is consistent
//third argument is ignored for now.
static boolean check_binary_signature_kmers(FILE * fp, uint32_t kmer_size, uint32_t binary_version, uint32_t* number_of_colours_in_binary, uint32_t* mean_read_len, uint64_t* total_seq)
{
    size_t read;
    char magic_number[6];
    boolean ret = false;
    
    read = fread(magic_number,sizeof(char),6,fp);
    if (read>0 &&
        magic_number[0]=='K' &&
        magic_number[1]=='M' &&
        magic_number[2]=='E' &&
        magic_number[3]=='R' &&
        magic_number[4]=='S' &&
        magic_number[5]==' ' )
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
                                    magic_number[0]=='K' &&
                                    magic_number[1]=='M' &&
                                    magic_number[2]=='E' &&
                                    magic_number[3]=='R' &&
                                    magic_number[4]=='S' &&
                                    magic_number[5]==' ' )
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
        printf("Binary does not have magic number in header. Corrupt, or not a KMERS binary\n");
    }
    
    return ret;
}

long long load_binary_from_filename_into_kmers_hash(char *filename, dBGraph * db_graph, short colour, boolean all_entries_are_unique)

{
    
	FILE *fp_bin = fopen(filename, "r");
    
	if (fp_bin == NULL) {
		printf("cannot open file:%s\n", filename);
		exit(1);	//TODO: prefer to print warning and skip file and return an error code?
	}
    
	uint32_t num_colours_in_binary=1;
	uint32_t mean_read_len;
	uint64_t total_seq;
	if ( !check_binary_signature_kmers(fp_bin, db_graph->kmer_size, BINVERSION, &num_colours_in_binary, &mean_read_len, &total_seq)){
        errx(1,"binary version or kmer_size are inconsistent");
	}
    hash_table_add_number_of_reads(total_seq,colour, db_graph);
    
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
        
		db_node_set_edges(current_node, colour, db_node_get_edges(&node_from_file));
        
		element_update_coverage(current_node, colour,
                                element_get_coverage_all_colours
                                (&node_from_file));
	}
    
	fclose(fp_bin);
    
	return count;
    
}

long long load_kmers_binary_from_filename_update_coverage(char *filename, dBGraph * db_graph, short colour){
    
	FILE *fp_bin = fopen(filename, "r");
    
	if (fp_bin == NULL) {
		printf("cannot open file:%s\n", filename);
		exit(1);	//TODO: prefer to print warning and skip file and return an error code?
	}
    
	uint32_t num_colours_in_binary=1;
	uint32_t mean_read_len;
	uint64_t total_seq;
	if ( !check_binary_signature_kmers(fp_bin, db_graph->kmer_size, BINVERSION, &num_colours_in_binary, &mean_read_len, &total_seq)){
        errx(1,"binary version or kmer_size are inconsistent");
	}
    hash_table_add_number_of_reads(total_seq,colour, db_graph);
    
	dBNode node_from_file;
	long long count = 0;
	BinaryKmer tmp_kmer;
    
	//Go through all the entries in the binary file
	while (db_node_read_binary
	       (fp_bin, db_graph->kmer_size, &node_from_file)) {
		count++;
        
		dBNode *current_node = NULL;		
        current_node = hash_table_find(element_get_key(element_get_kmer(&node_from_file), db_graph->kmer_size, &tmp_kmer), db_graph);
        if (current_node) {
            element_update_coverage(current_node, colour,
                                    element_get_coverage_all_colours
                                    (&node_from_file));
        }
	}
    
	fclose(fp_bin);
    
	return  count;
    
}

