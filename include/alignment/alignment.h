/**    aligner.h
 * 
 *  	A simple Smith-Watersman algorithm implementation to be used under the 
 * 		a similar way than cortex. 
 *      
 *      Copyright 2010 by M. Caccamo (mario.caccamo@bbsrc.ac.uk),  Z. Iqbal (zam@well.ox.ac.uk) and R. Ramirez-Gonzalez(Ricardo.Ramirez-Gonzalez@bbsrc.ac.uk)
 *      
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; either version 2 of the License, or
 *      (at your option) any later version.
 *      
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *      
 *      You should have received a copy of the GNU General Public License
 *      along with this program; if not, write to the Free Software
 *      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *      MA 02110-1301, USA.
 *
 */
#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#ifdef BLAST_MATRIX
#ifndef GAP_PENALTY
#define GAP_PENALTY -10
#endif

#ifndef TRANSITION_SCORE
#define TRANSITION_SCORE -4
#endif

#ifndef TRANSVERSION_SCORE
#define TRANSVERSION_SCORE -4
#endif

#ifndef MATCH_SCORE
#define MATCH_SCORE 5
#endif

#ifndef ALIGNMENT_SCORE_DATATYPE
#define ALIGNMENT_SCORE_DATATYPE int
#endif

#ifndef GAP
#define GAP '-'
#endif
#endif


#ifndef GAP_PENALTY
#define GAP_PENALTY -10
#endif

#ifndef TRANSITION_SCORE
#define TRANSITION_SCORE -5
#endif

#ifndef TRANSVERSION_SCORE
#define TRANSVERSION_SCORE -7
#endif

#ifndef MATCH_SCORE
#define MATCH_SCORE 2
#endif

#ifndef ALIGNMENT_SCORE_DATATYPE
#define ALIGNMENT_SCORE_DATATYPE int
#endif

#ifndef GAP
#define GAP '-'
#endif

typedef enum {NONE, DIAGONAL, LEFT, UP} Direction;
/**
 * Definition of the struct. To Be able to align in paralell
 * and to reduce the copy times, the sequence a and b shall be 
 * immutable. 
 * 
 */ 

typedef struct {
	Sequence * a;
	Sequence * b;
	Sequence * a_aligned;
	Sequence * b_aligned;
	//Sequence * consensus;
	int * a_mapping;
	int * b_mapping;
	int max_length;
	ALIGNMENT_SCORE_DATATYPE * * scores;
	int start_a, start_b, end_a, end_b, score,mismatches ,matches, length;
	int start_longest_a,start_longest_b,end_longest_a,end_longest_b, longest_lenght;
	boolean aligned;
} Alignment;

typedef struct{
	Sequence * reference;
	SequenceArray * sequences;
	Sequence * concensus;
	Alignment ** alignments;
	PeptideSequence * concensus_protein;
	PeptideSequence * reference_protein;
	int capacity;
	char qual_threashold;
	int min_alignment_score;
	short qual_window;
	int min_correct_bases;
	int max_name_length;
	int max_lenght;
}AlignmentToReference;

typedef struct{
	Sequence * reference;
	Sequence * sequence;
	Sequence * unmatched;
	int matches;
	int mismatches;
	int max_mismatches;
	int max_distance_between_mismatches;
	char search_window;
	int longest;
	int fix_count;
	int interations;
}MaskAlignment;

typedef struct{
	Sequence * sequence;
	BinaryTree * numps;
	BinaryTree * chimeras;
	
}RelatedSequences;


Alignment * alignment_new(int max_length, int max_name_lenght);

Alignment * alignment_set_sequences(Sequence * a, Sequence * b, Alignment * al);

Alignment * alignment_clean(Alignment * al);

void alignment_copy(Alignment * tgt, Alignment * src);

void alignment_destroy(Alignment **);

Alignment * alignment_align(Alignment * al, ALIGNMENT_SCORE_DATATYPE(* score_rule)(char a, char b ));

Alignment * alignment_consense(Alignment * al);

ALIGNMENT_SCORE_DATATYPE default_nucleotide_score(char a, char b);

ALIGNMENT_SCORE_DATATYPE half_gap_nucleotide_score(char a, char b);

Alignment * alignment_fix_homopolymers(Alignment * al);

void alignment_print_formated(FILE * f, Alignment * al);

Alignment * alignment_align_nucleotides(Alignment * al);

ALIGNMENT_SCORE_DATATYPE extreme_nucleotide_score(char a, char b);

//Functions for the aligment to reference
AlignmentToReference * alignment_to_reference_new(int alignments, int max_length, int max_name_lenght);

void alignment_to_reference_clean( AlignmentToReference * atr);

void alingment_to_reference_destroy(AlignmentToReference ** atr);

Sequence * alignment_to_reference_get_reference(AlignmentToReference * atr);

Sequence * alignment_to_reference_get_sequence(int pos, AlignmentToReference * atr);

void alignment_to_reference_align_all(AlignmentToReference *atr, ALIGNMENT_SCORE_DATATYPE(* score_rule)(char a, char b ));

void alignment_to_reference_print(FILE * f, AlignmentToReference * atr);

void alignment_to_reference_set_threadshodlds(int qual, int score, AlignmentToReference * atr);

void alignment_to_reference_consense_by_quality(AlignmentToReference * atr);

void alignment_to_reference_consence(AlignmentToReference *atr);

void alignment_to_reference_quality_filter_window(AlignmentToReference * atr);

void alignment_to_reference_print_stats_header(FILE * f);

boolean alignment_to_reference_print_stats(FILE * f,  AlignmentToReference * atr);

int alignment_to_reference_count_correct_bases(AlignmentToReference * atr);

void alignment_to_reference_mask_primer(AlignmentToReference *atr);

void alignment_to_reference_clean_reference( AlignmentToReference * atr);
//Functions for the masked alginment. 

MaskAlignment * mask_alignment_new(Sequence * reference);

void mask_alignment_destroy(MaskAlignment ** ma);

void mask_alignment_align(MaskAlignment * ma);

int mask_alignment_longest_match(MaskAlignment * ma);

void mask_alignment_set_sequence(Sequence * seq, MaskAlignment * ma);

void mask_alignment_print_alignment_fasta(FILE * f, MaskAlignment * ma);

int mask_alignment_compare( void * a,  void * b);

int mask_alignment_simple_compare( void * a,  void * b);

void mask_alignment_print_stats_header(FILE * f);

void mask_alignment_print_print_stats(FILE * f,  MaskAlignment * ma);

void mask_alignment_fix_homopolymers(MaskAlignment * ma);

void mask_alignment_fix_insertions(MaskAlignment * ma);

void mask_alignment_fix_deletions(MaskAlignment * ma);

void mask_alignment_fix_local_del_in(MaskAlignment * ma);

void mask_alignment_fix_local_in_del(MaskAlignment * ma);



//Functions to find related sequences. 

RelatedSequences * related_sequences_new(Sequence * seq);

void related_sequence_destroy(RelatedSequences ** rs);
double nump_score(Alignment * al, double identity_threshold, double  lenght_window);
double chimera_score(Alignment * al, double identity_threshold, double  min_lenght, double max_lenght);
void related_sequences_to_fasta(FILE * f, RelatedSequences * rs);
#endif
