/*
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
 * This file is to declare the analysis functions. Stuff that is not necesarily part of the main assembly, but that can be useful to analyze the graph. 
 */

#import <stdio.h>

#include <global.h>
#include <flags.h>
#include <binary_kmer.h>
#include <element.h>
#include <dB_graph.h>
#include <file_reader.h>
#include <path.h>
#include <perfect_path.h>
#include <branches.h>
#include <y_walk.h>
#include <analysis.h>
#include <logger.h>

void analysis_print_reference_coverage(char * reference_file, char * output_file, dBGraph * db_graph){
	log_write_timestamp(true);
	log_and_screen_printf("\nPrinting the coverage from a reference.\n");
	fflush(stdout);
	
	//long long kmers = seq_kmer_iterator(FILE * fp, void * args, void  (*kmer_action) (BinaryKmer * k, void * arg) ,int (*file_reader) (FILE * fp, Sequence * seq,
	//					int max_read_length, boolean new_entry, boolean * full_entry),
	//		    		long long *bad_reads, int fastq_ascii_offset, char quality_cut_off,
//						int max_read_length, short kmer_size);

	//TODO: Implement this!
}

void analysis_init_cortex_environment(short kmer_size){
    path_array_initialise_buffers(kmer_size);
}



