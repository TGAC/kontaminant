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

/*
  cmd_line.h - manipulation of command line
  
*/

#ifndef CMD_LINE_H_
#define CMD_LINE_H_

#include <stdio.h>
#include <global.h>
#include <file_format.h>

typedef enum{
	PERFECT_PATH = 0,
	BRANCHES = 1,
	Y_WALK = 2,
	BUBBLES = 3,
	READ_PAIR = 5,
    METACORTEX = 6
}TraverseAlgorithm;

typedef struct
{
    //core parameters 
    boolean verbose;
    int kmer_size;
    int bucket_size;
    int number_of_buckets_bits;
    
    //----------------
    //actions on/off
    //----------------
    //input
    boolean input_file; //if it is present
    boolean input_file_format_known; 
    boolean health_check_binary; 
	boolean input_reference_known;
    
    //cleaning
    boolean tip_clip;
    boolean remove_bubbles;
    boolean low_coverage_path_clip;
    boolean low_coverage_node_clip;
    boolean remove_low_coverage_supernodes; 
    boolean remove_spurious_links;
    
    //output 
    boolean dump_binary;
    boolean output_fasta;
    boolean dump_hash;
    boolean detect_bubbles; 
    boolean output_coverages;
    boolean graphviz;
    boolean print_uncertain_as_n;
    boolean output_log;
	boolean output_kmer_coverage_known;
	boolean output_reference_coverage_file_known;
	boolean output_kmer_coverage_know ;
    //read pairs
#if defined (ENABLE_MARK_PAIR) || defined ( ENABLE_READ_PAIR)
    boolean read_pair_enabled;
#endif

    //colour space
#ifdef SOLID
    boolean output_base_space;
#endif
    
    //-----------
    //parameters 
    //-----------
    //core parameters
    int threads;
    short max_double_y_complexity; 
    
    //input
    FileFormat input_file_format;
    int quality_score_threshold;
    char input_filename[LENGTH_FILENAME];
    char qual_filename[LENGTH_FILENAME];
	char input_reference[LENGTH_FILENAME];
	
    int quality_score_offset;
    int max_read_len; 
    int binary_version;
    
    //cleaning
    int node_coverage_threshold;
    int tip_length;
    int remove_low_coverage_supernodes_threshold;
    int bubble_max_depth;
    int bubble_max_length;
    int remove_spurious_links_min_coverage;
    int remove_spurious_links_max_difference;
    int path_coverage_threshold;
    int tip_clip_iterations;
    
    //output
    char output_ctx_filename[LENGTH_FILENAME];
    char output_fasta_filename[LENGTH_FILENAME];
    char output_graphviz_filename[LENGTH_FILENAME];
    char output_hash_filename[LENGTH_FILENAME];
    char log_filename[LENGTH_FILENAME];
    char output_kmer_coverage[LENGTH_FILENAME];
	char output_reference_coverage_file[LENGTH_FILENAME];
	int max_length;
    int singleton_length;
    int min_subgraph_size;
    TraverseAlgorithm algorithm;
    
    //read pairs
#if defined (ENABLE_MARK_PAIR) || defined ( ENABLE_READ_PAIR)
    int read_pair_distance;
    int read_pair_coverage;
    int read_pair_tolerance;
    int read_pair_min_bits;
    int read_pair_start_length;
    char read_pair_filename[LENGTH_FILENAME];
    int read_pair_max_paths;
    int read_pair_min_kmers;
#endif
    
#ifdef ENABLE_MARK_PAIR
    char mark_pair_filename[LENGTH_FILENAME];
#endif
    
    

} CmdLine;

CmdLine parse_cmdline( int argc, char* argv[],int unit_size); 
int default_opts(CmdLine *);

#endif /* CMD_LINE_H_ */
