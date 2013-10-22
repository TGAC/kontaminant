//
//  kmer_hash.h
//  Untitled
//
//  Created by Ricardo Ramirez-Gonzalez (TGAC) on 08/12/2011.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef kmer_hash_h
#define kmer_hash_h

#define KMER_HASH_SAMPLE_INDEX 1
#define KMER_HASH_REFERENCE_INDEX 0
typedef  dBGraph KmerHash; //TODO make this more independent, and dont depend from the graph
void kmer_hash_reset_coverage( KmerHash * kmer_hash );
void kmer_hash_calculate_stats(KmerHash * kmer_hash);
int kmer_hash_present_kmers_in_sliding_window_set(KmerSlidingWindowSet * ksws, KmerHash * kmer_hash);
void kmer_hash_add_contaminated_reads(long long count, KmerHash * kmer_hash);
void kmer_hash_print_kmer_stats(char * output_filename,char * reference_kmers_file,char * input_filename, KmerHash * kmer_hash);
void kmer_hash_print_contaminated_kmers_histogram(char * filename, KmerHash * kmer_hash);
#endif
