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
#include <nucleotide.h>
#include <binary_kmer.h>
#include <kmer_hash.h>
#include <stdio.h>


static void kmer_reset_coverage(Element * e){
    e->coverage[KMER_HASH_SAMPLE_INDEX] = 0;
}

void kmer_hash_reset_coverage( KmerHash * kmer_hash )
{
	printf("Reseting coverage flags...\n");
#ifdef THREADS
	hash_table_threaded_traverse(&kmer_reset_coverage, kmer_hash);
#else
	hash_table_traverse(&kmer_reset_coverage, kmer_hash);
#endif
	printf("\n");
	fflush(stdout);
}

void kmer_hash_calculate_stats(KmerHash * kmer_hash){
    db_graph_calculate_stats(kmer_hash);
}

int kmer_hash_present_kmers_in_sliding_window_set(KmerSlidingWindowSet * windows, KmerHash * kmer_hash){
    int count = 0;
    int i, j;
    BinaryKmer tmp_kmer;
    Element * e;
    int s = 0;
    short kmer_size = kmer_hash->kmer_size;
    for (i = 0; i < windows->nwindows; i++) {	//for each window
        KmerSlidingWindow * current_window = &(windows->window[i]);
        for (j = 0; j < current_window->nkmers; j++) {	//for each kmer in window
            Key key = element_get_key(&(current_window->kmer[j]),kmer_size, &tmp_kmer);
            e = hash_table_find(key, kmer_hash);   
            s++;
            if (e != NULL) {
                count++;
            }
        }
    }
    
    return count;
}

void kmer_hash_add_contaminated_reads(long long count, KmerHash * kmer_hash){
#ifdef KMER_TOOLS 
    kmer_hash->contaminated_reads += count;
#else 
    fprintf(stderr, "Unsupported function kkmer_hash_add_contaminated_reads");
    exit(-1);
#endif
}

void kmer_hash_print_kmer_stats(char * output_filename,char * reference_kmers_file,char * input_filename, KmerHash * kmer_hash){
    
    kmer_hash_calculate_stats( kmer_hash);
    
    long long reference_kmers =  kmer_hash->colour_kmers[KMER_HASH_REFERENCE_INDEX];
    double reference_coverage =  kmer_hash->average_coverage[KMER_HASH_REFERENCE_INDEX];
    long long sample_kmers = kmer_hash->colour_kmers[KMER_HASH_SAMPLE_INDEX];
    long long common_kmers =  kmer_hash->common_kmers_in_all_colours; 
    double sample_coverage = kmer_hash->average_coverage[KMER_HASH_SAMPLE_INDEX];
    long long sample_reads = kmer_hash->number_of_reads[KMER_HASH_SAMPLE_INDEX];
    long long contaminated_reads = kmer_hash->contaminated_reads;
    log_and_screen_printf("\nReference kmers:  %'lld", reference_kmers);	
    log_and_screen_printf("\nReference coverage:  %f", reference_coverage);
    log_and_screen_printf("\nSample kmers:  %'lld",  sample_kmers);	
    log_and_screen_printf("\nSample coverage:  %f",  sample_coverage);	
    log_and_screen_printf("\nCommon kmers:  %'lld",  common_kmers);	
    
    
    
    
    float percentage = 100 * ( (float) common_kmers / (float)  kmer_hash->colour_kmers[KMER_HASH_REFERENCE_INDEX] ) ;
    float percentage_reads = 100 * ( (float) contaminated_reads / (float)  kmer_hash->number_of_reads[KMER_HASH_SAMPLE_INDEX] ) ;
    
    log_and_screen_printf("\nPercemtage of kmers sample in reference:  %f",  percentage);
    log_and_screen_printf("\nContaminated reads:  %'lld",  contaminated_reads);	
    log_and_screen_printf("\nSample reads count:  %'lld",   kmer_hash->number_of_reads[KMER_HASH_SAMPLE_INDEX] );	
    log_and_screen_printf("\nPercentage of contaminated reads:  %f",  percentage_reads);	
    
    FILE * output = NULL;
    struct flock fl = {F_WRLCK, SEEK_SET,   0,      0,     0 };
    fl.l_pid = getpid();
    short attemps = 5;
    if(output_filename){
        if (access(output_filename, F_OK) == -1){
            output = fopen(output_filename, "a");
            //We print the header here
            if(output == NULL){
                log_and_screen_printf("Unable to open file: %s\n");
                return;
            }
            while (fcntl(fileno(output), F_SETLKW, &fl) == -1 && attemps-- > 0) {
                sleep(rand()%10);//Sleep from 0 to 10 seconds...
            }
            if (attemps == 0) {
                log_and_screen_printf("Unable to lock file: %s\n");
            }
            
            fprintf(output, "REFERENCE\tSAMPLE\tREFERENCE_READS\tSAMPLE_READS\tKMER_SIZE\tREFERENCE_KMERS\tREFERENCE_COVERAGE\tSAMPLE_KMERS\tSAMPLE_COVERAGE\tPERCENTAGE\tCOMMON\tSAMPLE_READS\tCONTAMINTED_READS\tPERCENTAGE_CONTAMINATED_READS\n");
        }else{
            output = fopen(output_filename, "a");
            if(output == NULL){
                log_and_screen_printf("Unable to open file: %s\n");
            }
            
            while (fcntl(fileno(output), F_SETLKW, &fl) == -1 && attemps-- > 0) {
                sleep(rand()%10);//Sleep from 0 to 10 seconds...
            }
            if (attemps == 0) {
                log_and_screen_printf("Unable to lock file: %s\n");
            }
            
        }
        
        fprintf(output, "%s\t%s\t%lld\t%lld\t%d\t%lld\t%f\t%lld\t%f\t%f\t%lld\t%lld\t%lld\t%f\n", reference_kmers_file, input_filename, hash_table_get_number_of_reads(KMER_HASH_REFERENCE_INDEX, kmer_hash),hash_table_get_number_of_reads(KMER_HASH_SAMPLE_INDEX, kmer_hash), kmer_hash->kmer_size, reference_kmers, reference_coverage, sample_kmers, sample_coverage, percentage, common_kmers, sample_reads, contaminated_reads, percentage_reads);
        
	}
    if (output!=NULL) {
        fflush(output);
        fl.l_type=F_UNLCK;
        
        if (output != NULL && fcntl(fileno(output), F_SETLKW, &fl) == -1 ) {
            log_and_screen_printf("ERROR: Unable to unlock!");
        }
        fclose(output);

    }
   }

void kmer_hash_print_contaminated_kmers_histogram(char * filename, KmerHash * kmer_hash){

    if (filename == NULL) {
        return;
    }
    
    FILE * output = fopen(filename, "w");
    if(output == NULL){
        log_and_screen_printf("Unable to open file: %s\n", filename);
        return;
    }
    
    int max_present = 0, i=0;
    
    for (i =0; i  < MAX_READ_LENGTH; i++) {
        if (kmer_hash->contaminated_kmers_per_read[i]) {
            max_present = i;
        }
    }
    fprintf(output, "NO_OF_KMERS\tCOUNT\n");
    for (i=0; i<=max_present; i++) {
        fprintf(output, "%d\t%lld\n",i,kmer_hash->contaminated_kmers_per_read[i]);
    }
    fclose(output);
    
}

