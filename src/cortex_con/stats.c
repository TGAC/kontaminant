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
#include <flags.h>
#include <limits.h>
#include <flags.h>
#include <nucleotide.h>
#include <binary_kmer.h>
#include <element.h>
#include <logger.h>
#include <seq.h>

#include <cmd_line.h>
int main(int argc, char **argv)
{
	CmdLine cmd_line = parse_cmdline(argc, argv, sizeof(Element));

	int max_read_length = MAX_READ_LENGTH;

	FILE *fp = fopen(cmd_line.input_filename, "r");
	if (fp == NULL) {
		fprintf(stderr, "cannot open file:%s\n",
			cmd_line.input_filename);
		exit(1);	//TODO - prefer to print warning and skip file and return an error code?
	}
	
	fseek(fp, 0L, SEEK_END);
	long long sz = ftell(fp);
	fseek(fp, 0L, SEEK_SET);
	
	//----------------------------------
	// preallocate the memory used to read the sequences
	//----------------------------------
	Sequence *seq = malloc(sizeof(Sequence));
	if (seq == NULL) {
		fputs("Out of memory trying to allocate Sequence\n", stderr);
		exit(1);
	}
	alloc_sequence(seq, max_read_length, LINE_MAX, cmd_line.quality_offset);

	int entry_length = 0;
	boolean new_entry = true;
	boolean full_entry;
	SequenceStats totalStats;
	clean_stats(&totalStats);
	FILE *fq;
	long long base_count = 0;
	// printf("Fileformat %d",cmd_line.input_file_format );
	printf("\n");
	log_progress_bar(0);
	long long seq_count = 0;
	
	do {
		if(cmd_line.verbose){
			if(seq_count++ % 10000 == 0){
				log_progress_bar((ftell(fp)*100)/sz );
			}
		}
		switch (cmd_line.input_file_format) {
		case FASTQ:
			entry_length =
			    read_sequence_from_fastq(fp, seq, max_read_length);
			break;

		case FASTA:

			entry_length =
			    read_sequence_from_fasta(fp, seq, max_read_length,
						     new_entry, &full_entry, 0);
			new_entry = full_entry;
			break;
		case ROCHE:
			// printf("WE ARE IN THE CASE!");
			fq = fopen(cmd_line.qual_filename, "r");
			if (fq == NULL) {
				fprintf(stderr, "cannot open file:%s\n",
					cmd_line.qual_filename);
				exit(1);	//TODO - prefer to print warning and skip file and return an error code?
			}
			
			entry_length =
			    read_sequence_from_fasta_and_qual(fp, fq, seq,
							      max_read_length);
			new_entry = full_entry;
			break;
		default:
			fprintf(stderr, "Unsupported format for statistics");
			exit(-1);
			break;
		}

		
		if(cmd_line.low_coverage_node_clip){
			clean_stats(&totalStats);
			sequence_stats(&totalStats, seq);
			
			double th = (double)cmd_line.node_coverage_threshold / 100;
			
			double base_c = base_content(Adenine, &totalStats)  ; 
			
			printf("%f < %f \n", base_c, th);
			
			if(base_c < th){
				FILE *fp = fopen("less.fa", "a");
				sequence_print_fasta(fp, seq);
				fclose(fp);
			}else{
				FILE *fp = fopen("more.fa", "a");
				sequence_print_fasta(fp, seq);
				fclose(fp);
			}
			
		}else if(cmd_line.remove_bubbles ){
			base_count += seq->length;
		}else{
			sequence_stats(&totalStats, seq);
		}
		
		/*int i;

		 */
	}
	while (entry_length > 0);
	log_progress_bar(100);
	printf("\n");
	fclose(fp);
	fp = fopen(cmd_line.output_fasta_filename, "w");
	if (fp == NULL) {
		fprintf(stderr, "Unable to open output file");
		exit(-1);
	}
	
	 if(cmd_line.remove_bubbles ){		
		fprintf(fp, "Bases count\t %lld\n", base_count);
	}else{
		print_stats(fp, &totalStats);
	}
	fclose(fp);

	return 0;
}
