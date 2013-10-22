/**    aligner.c
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


#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>

#include <nucleotide.h>
#include <binary_kmer.h>
#include <seq.h>
#include <peptide.h>
#include <binary_tree.h>
#include <alignment.h>

#include <string.h>
/**
 * Preallocates the space required for the alignment, and the sequences. 
 * It assumes that the sequences to align will be allocated somewhere else.
 * This is to avoid the need of coping the sequences back and forth many times.   
 * This also implies that the sequences have to be allocated once an not modified
 * by any thread. 
 * The size of the buffers are twice the max_lenght since it represents
 * the max_lenght of each individual sequence. 
 * 
 * WARNING. The max_lenght has to be consisten with the lenght of the sequences being read. 
 * WARNING. The design implies that the Alighment objects will be recycled when possible. 
 */

Alignment * alignment_new(int max_length, int max_name_length){
	Alignment * al = calloc(1,sizeof(Alignment));
	al->max_length = max_length;
	//al->max_name_length = max_name_length;
	
	al->a = NULL;
	al->b = NULL;
	
	al->a_aligned = calloc(1, sizeof(Sequence));
	
	al->b_aligned = calloc(1, sizeof(Sequence));
	//al->consensus = calloc(1, sizeof(Sequence));
	if(al->a_aligned == NULL || al->b_aligned == NULL ){
		fprintf(stderr, "[alignment_new] Unable to allocate memory for the alignments\n");
		exit(-1);
	} 
	alloc_sequence(al->a_aligned,  max_length, max_name_length, 0);
	alloc_sequence(al->b_aligned, max_length, max_name_length, 0);
	//	alloc_sequence(al->consensus, max_length, max_name_length);
	
	al->scores = calloc(max_length, sizeof(ALIGNMENT_SCORE_DATATYPE *));
	
	if(al->scores == NULL){
		fprintf(stderr, "al->scores[i]");
		exit(-1);
	}
	
	/*al->directions = calloc(max_length, sizeof(Direction *));
	 
	 if(al->directions == NULL){
	 fprintf(stderr, "al->directions[i]");
	 exit(-1);
	 }*/ 
	
	al->a_mapping = calloc(max_length, sizeof(int));
	al->b_mapping = calloc(max_length, sizeof(int));
	
	int i;
	for(i=0;i<max_length;i++){
		al->scores[i] = calloc(max_length, sizeof(ALIGNMENT_SCORE_DATATYPE));
		//al->directions[i] = calloc(max_length, sizeof(Direction));
		if(al->scores[i]==NULL){
			fprintf(stderr, "[alignment_new] Unable to allocate memory for al->scores[%d]\n",i);
			exit(-1);
		}
	}
	al->max_length = max_length;
	return al;
}

Alignment * alignment_set_sequences(Sequence * a, Sequence * b, Alignment * al){
	al->a = a;
	al->b = b;
	return al;
}


double  alignment_identity(Alignment * al){
	assert(al != NULL);
	
	return (double) al->matches/(double)(al->length);
}

Alignment * alignment_clean(Alignment * al){
	/*al->a = NULL;
	 al->b = NULL;
	 */ 
	sequence_clean(al->a_aligned);
	sequence_clean(al->b_aligned);
	//  sequence_clean(al->consensus);
	int i, j, max_length = al->max_length;
	
	for(i=0;i<max_length;i++){
	    al->a_mapping[i] = -1;
	    al->b_mapping[i] = -1;
		for(j=0;j<max_length;j++){
			al->scores[i][j]=0;
			//al->directions[i][j] = NONE;
		}
	}
	return al;
}

/**
 * Copy an alignment src to tgt. It requires that all the pointers
 * of the tgt are already initalized. Exceptions are the original
 * sequences that will still ve pointer to the original sequences.			
 *
 */
void alignment_copy(Alignment * tgt, Alignment * src){
	
	int i, j;
	
	assert(src != NULL);
	assert(tgt != NULL);
	assert(src->max_length == tgt->max_length);
	
	for (i = 0; i < src->max_length; i++) {
		tgt->a_mapping[i] = src->a_mapping[i];
		tgt->b_mapping[i] = src->b_mapping[i];
	}
	
	for (i=0;  i < src->max_length; i++) {
		for (j=0; j < src->max_length; j++) {

			tgt->scores[i][j] = src->scores[i][j];
		}
	}
	
	tgt->a = src->a;
	tgt->b = src->b;
	sequence_copy(tgt->b_aligned, src->b_aligned);
	sequence_copy(tgt->a_aligned, src->a_aligned);
	
	tgt->max_length = src->max_length;
	
	tgt->start_a = src->start_a;
	tgt->start_b = src->start_b;  
	tgt->end_a = src->end_a;
	tgt->end_b = src->end_b;
	tgt->score = src->score;
	tgt->mismatches = src->mismatches;
	tgt->matches = src->matches;
	tgt->length = src->length;
	tgt->start_longest_a = src->start_longest_a;
	tgt->start_longest_b = src->start_longest_b;
	tgt->end_longest_a = src->end_longest_a;
	tgt->end_longest_b = src ->end_longest_b;
	tgt->longest_lenght = src->longest_lenght;
	tgt->aligned= src->aligned;
	
}

void alignment_destroy(Alignment ** al){
	Sequence * tmp =(*al)->a_aligned ;
	
	free_sequence(&tmp);
	tmp =(*al)->b_aligned ;
	free_sequence(&tmp);
	
	
	int i;
    for(i = 0; i < (*al)->max_length; i++)
	{
		free((*al)->scores[i]);
	}
	free((*al)->scores);
	free((*al)->a_mapping);
	free((*al)->b_mapping);
	free(*al);
	*al = NULL;
}

Alignment * alignment_align_nucleotides(Alignment * al){
	alignment_align(al, &default_nucleotide_score);
	return al;
	
}

static ALIGNMENT_SCORE_DATATYPE simple_score(char a, char b){
	
	if(a == GAP ||  b == GAP){
		return GAP_PENALTY ;
	}
	if(a == b){
		return MATCH_SCORE;
	}else{//MISMATCH
		return TRANSITION_SCORE;
	}
}

static boolean is_transversion(char a, char b){
	if(a == 'C' || a == 'T')
		if(b == 'A' || b == 'G')
			return true;
	
	if(a == 'A' || a == 'G')
		if(b == 'C' || b == 'T')
			return true;	
	
	return false; 
}

static boolean is_transition(char a, char b){
	if(a == 'A' && b == 'G' )
		return true;
	
	if(a == 'G' && b == 'A' )
		return true;
	
	if(a == 'C' && b == 'T' )
		return true;
	
	if(a == 'T' && b == 'C' )
		return true;
	return false;
}

ALIGNMENT_SCORE_DATATYPE default_nucleotide_score(char a, char b){
	
	if(a == GAP ||  b == GAP){
		return GAP_PENALTY ;
	}else
		if(a == b){
			return MATCH_SCORE;
		}else
			if(is_transition(a, b))
				return TRANSITION_SCORE;
			else if(is_transversion(a,b))
				return TRANSVERSION_SCORE;
			else 
				return MATCH_SCORE;
}

ALIGNMENT_SCORE_DATATYPE extreme_nucleotide_score(char a, char b){
	
	if(a == GAP ||  b == GAP){
		return GAP_PENALTY * 2 ;
	}else
		if(a == b){
			return MATCH_SCORE / 2;
		}else
			if(is_transition(a, b))
				return TRANSITION_SCORE;
			else if(is_transversion(a,b))
				return TRANSVERSION_SCORE;
			else 
				return MATCH_SCORE;
}

ALIGNMENT_SCORE_DATATYPE half_gap_nucleotide_score(char a, char b){
	
	if(a == GAP ||  b == GAP){
		return GAP_PENALTY /2;
	}else
		if(a == b){
			return MATCH_SCORE;
		}else
			if(is_transition(a, b))
				return TRANSITION_SCORE;
			else if(is_transversion(a,b))
				return TRANSVERSION_SCORE;
			else 
				return GAP_PENALTY;
}


void score_matrix_iterator(Alignment * al, void (* function)(ALIGNMENT_SCORE_DATATYPE)){
	int i, j;
	for(i = 0; i < al->a->length; i++){
		for(j = 0; j < al->b->length; j++){
			function(al->scores[i][j]);
		}
	}
}

char alignment_direction_to_char(Direction d){
	switch(d){
		case UP: return 'u';
		case DIAGONAL: return 'd';
		case LEFT: return 'l';
		case NONE: return '-';
		default: return ' ';
	}
} 

void print_score_matrix(FILE * f, Alignment * al){
	
#ifndef NO_BOUNDS_CHECK
	if(al == NULL){
		fprintf(stderr, "[print_score_matrix] Alignment cant be null!!!!");
		exit(-1);
	}
#endif
	int i, j;
	
	fprintf(f, "Score matrix %s-%s\n-\t", al->a->name, al->b->name);
	
	for(j = 0; j < al->b->length; j++){
		fprintf(f, "%c\t", al->b->seq[j]);
	}
	fprintf(f, "\n");
	for(i = 0; i < al->a->length; i++){
		fprintf(f,  "%c\t", al->a->seq[i]);
		for(j = 0; j < al->b->length; j++){
			fprintf(f, "%d\t", al->scores[i][j]);
		}
		fprintf(f, "\n");
	}
}

Alignment * alignment_align(Alignment * al, ALIGNMENT_SCORE_DATATYPE(* score_rule)(char a, char b )){
#ifndef NO_BOUNDS_CHECK	
	if(al== NULL){
		fprintf(stderr, "[alignment_align] The alignment can't be NULL");
		exit(-1);
	}
	
	if(al->a == NULL){
		fprintf(stderr, "[alignment_align] The sequence a for the alignment can't be NULL");
		exit(-1);
	}
	
	if(al->b == NULL){
		fprintf(stderr, "[alignment_align] The sequence b for the alignment can't be NULL");
		exit(-1);
	}
#endif
	if(!al->a->upper_case)
		sequence_to_upper_case(al->a);
	if(!al->b->upper_case)
		sequence_to_upper_case(al->b);
	
	if(score_rule== NULL){
		score_rule = &simple_score;
	}
	
	int i, j;
	int current_a;
	int current_b;
	ALIGNMENT_SCORE_DATATYPE max_score = 0, max_i = 0, max_j = 0, left, up, left_up, tmp;
	Direction d;
	//Generate the scores matrix
	for(i = 0; i < al->a->length; i++){
		for(j = 0; j < al->b->length; j++){
			//Initializing the values to 0, in case we don't have them quite yet
			left_up =  0;
			left = GAP_PENALTY; 
			up =  GAP_PENALTY;
			tmp = -2147483648;
		    d = NONE;
		    
		    if(i > 0 && j >  0){//We are not in a border
				left_up =  al->scores[i-1][j-1];
			}
			if(i > 0){//We are not in the top
				up = al->scores[i-1][j];
			}
			if(j > 0){//We are in the left side
				left = al->scores[i][j-1];
			}
			
			left_up += score_rule(al->a->seq[i], al->b->seq[j]);
			up += score_rule(al->a->seq[i], GAP);
			left += score_rule(GAP, al->b->seq[j]);
			
			
			
			if(left_up > tmp){
				tmp = left_up;
				d = DIAGONAL;
			}
			if(left > tmp){
				tmp = left;
				d = LEFT;
			}
			if(up > tmp){
				tmp = up;
				d = UP;
			}
			
			if(tmp < 0){
				tmp = 0;
				d = NONE;
			} 
			
			
			
			al->scores[i][j] = tmp;
			//al->directions[i][j] = d;
			if (tmp >= max_score){		
				max_score = tmp;
				max_j = j;
				max_i = i;
			}
		}
	}
	//print_score_matrix(stdout, al);
	
	//Trackback from the highest score
	i = max_i;
	j = max_j;
	al->end_a = i ;
	al->end_b = j ;
	tmp = max_score;
	al->score = max_score;
	sequence_add_base(al->a->seq[i], al->a->qual[i], al->a_aligned);
	sequence_add_base(al->b->seq[j], al->b->qual[j], al->b_aligned);
	
	current_a = 0;
	current_b = 0;
	al->matches = 0; 
	al->length = 0;
	al->start_longest_a = 0;
	al->start_longest_b = 0;
	al->end_longest_a   = 0;
	al->end_longest_b   = 0;
	al->longest_lenght  = 0;
	int tmp_len = 0;
	int tmp_first_a = 0;
	int tmp_first_b = 0;
	
	do{
		
		left_up = 0;
		up = 0;
		left =   0;
		//printf("%d ", max_score);
		d = NONE;
		//d = al->directions[i][j];	
		if(i > 0 && j >  0){//We are not in a border
			left_up =  al->scores[i-1][j-1];
			
		}
		if(i > 0){//We are not in the top
			up = al->scores[i-1][j];
		}
		if(j > 0){//We are in the left side
			left = al->scores[i][j-1];
		}	
		
		if(left_up >= up && left_up >= left){//it is a match
			//if(d==DIAGONAL){
			
			
			
			i--;
			j--;
			al->a_mapping[current_a++] = i; //The length grows in each iteration. to retrive back when the
			al->b_mapping[current_b++] = j; //sequence is reversed, use the complement of this numbers
			
			//fprintf(stderr, "%d ", current_b);
			sequence_add_base(al->a->seq[i], al->a->qual[i], al->a_aligned);
			sequence_add_base(al->b->seq[j], al->b->qual[j], al->b_aligned);
			//	current_a++;
			//	current_b++;
			max_score = left_up;
			if (al->a->seq[i] == al->b->seq[j]) {
				al->matches++;
				
				
				if (tmp_len == 0) {
					tmp_first_a = i;
					tmp_first_b = j;
				}
				tmp_len++;
				
			}else {
				al->mismatches++;
				if(tmp_len > al->longest_lenght){
					al->start_longest_a = tmp_first_a- tmp_len;
					al->start_longest_b = tmp_first_b- tmp_len;
					al->longest_lenght  = tmp_len;
					al->end_longest_a = tmp_first_a ;
					al->end_longest_b = tmp_first_b ;
					
					tmp_len = 0;
					tmp_first_a = 0;
					tmp_first_b = 0;
				}
			}
			
			
		}else if(up > left){//we go up
			//}else if(d == UP){
			i--;
			//al->a_mapping[current_a++] =i;
			//current_a++;
			sequence_add_base(al->a->seq[i], al->a->qual[i], al->a_aligned);
			sequence_add_base(GAP, '-', al->b_aligned);
			
			if(tmp_len > al->longest_lenght){
				al->start_longest_a = tmp_first_a;
				al->start_longest_b = tmp_first_b;
				al->longest_lenght  = tmp_len;
				al->end_longest_a = tmp_first_a + tmp_len;
				al->end_longest_b = tmp_first_b + tmp_len;
				
				tmp_len = 0;
				tmp_first_a = 0;
				tmp_first_b = 0;
			}
			
			max_score = up ;
		}else /*if(d == LEFT)*/{//we go left
			j--;
			sequence_add_base(GAP, '-', al->a_aligned);
			sequence_add_base(al->b->seq[j], al->b->qual[j], al->b_aligned);
			max_score = left;
			if(tmp_len > al->longest_lenght){
				al->start_longest_a = tmp_first_a;
				al->start_longest_b = tmp_first_b;
				al->longest_lenght  = tmp_len;
				al->end_longest_a = tmp_first_a + tmp_len;
				al->end_longest_b = tmp_first_b + tmp_len;
				
				tmp_len = 0;
				tmp_first_a = 0;
				tmp_first_b = 0;
			}
		}
		al->length++;
	}while(i > 0 && j > 0 && max_score > 0/* && d != NONE*/); 
	
	al->start_a = i ;
	al->start_b = j ;
	al->aligned = true;
	sequence_reverse(al->a_aligned);
	sequence_reverse(al->b_aligned);
	
	//printf("\nAligned: \n%s\n%s\n", al->a_aligned->seq, al->b_aligned->seq);
	//printf("\nQualities: \n%s\n%s\n", al->a_aligned->qual, al->b_aligned->qual);
	//printf("aligned:");
	//sequence_print_fastq(stdout, al->a_aligned);
	//sequence_print_fastq(stdout, al->b_aligned);
	
	return al;
}

void print_mappings(Alignment * al){
	int i = 0;
	for(i = 0; i < al->b->length; i++){
		fprintf(stderr, "%d ", al->b_mapping[i]);
	}
}


Alignment * alignment_consense(Alignment * al){
	return al;
}


void alignment_print_formated(FILE * f, Alignment * al){
#ifndef NO_BOUNDS_CHECK	//For debugging purposes, in "production", if coded properly, we can skip the validation
	if(al == NULL){
		fprintf(stderr, "[alignment_print_formated] aligment can't be null");
		exit(-1);
	}
	if(f == NULL){
		fprintf(stderr, "[alignment_print_formated] File can't be null");
		exit(-1);
	}
	if(!al->aligned){
		if(f == NULL){
			fprintf(stderr, "[alignment_print_formated] You haven't aligned yet");
			exit(-1);
		}
	}
#endif 
    int length = sequence_get_length(al->a_aligned);
	int i, block_end = 80, block_start = 0;
	char current_a, current_b;
	
	
	int offset_a, offset_b;
	
	if(length > 1)
		do{
			offset_a = 0;
			offset_b = 0;
			fprintf(f, "%s (score: %d)\n %d \t", al->a->name,al->score ,al->start_a + block_start + 1);
			for(i = block_start; i < length && i < block_end; i++){
				current_a = sequence_get_base(i, al->a_aligned);
				if(current_a != '-'){
					offset_a++;
				}
				fprintf(f,"%c", current_a);
			}
			fprintf(f, "\t%d \n \t",al->start_a + block_start+offset_a );
			
			if (block_start > 10) {
				//fprintf(f, "\t",al->start_a + block_start+offset_a );
				fprintf(f, "\t");
			}
			
			
			for(i = block_start; i < length && i < block_end; i++){
				current_a = sequence_get_base(i, al->a_aligned);
				current_b = sequence_get_base(i, al->b_aligned);
				fprintf(f,"%c", current_a == current_b?'|': current_a == '-' || current_b=='-'? ' ' : '*');
			}
			fprintf(f, "\n %d \t", al->start_b + block_start + 1);
			for(i = block_start; i < length && i < block_end ; i++){
				current_b = sequence_get_base(i, al->b_aligned);
				if(current_b != '-'){
					offset_b++;
				}
				fprintf(f,"%c", current_b);
			}
			fprintf(f, "\t%d\n%s\n\n", al->start_b + block_start + offset_b , al->b->name);
			block_start += 80;
			block_end += 80;
			if(block_end > length){
				block_end = length;
			}
		}while(block_start < length);
	
}

AlignmentToReference * alignment_to_reference_new(int alignments, int max_length, int max_name_length){
	AlignmentToReference * tmp = calloc(1, sizeof(AlignmentToReference));
	tmp->capacity = alignments;
	tmp->reference = sequence_new(max_length, max_name_length,0);
	tmp->concensus = sequence_new(max_length, max_name_length,0);
	tmp->sequences = sequence_array_new(alignments);
	tmp->alignments = calloc(alignments, sizeof(Alignment *));
	tmp->concensus_protein = peptide_sequence_new(1+(max_length/3) , max_name_length,0);
	tmp->reference_protein = peptide_sequence_new(1+(max_length/3) , max_name_length,0);
	tmp->max_name_length = max_name_length;
	tmp->max_lenght = max_length;
	int i;
	for(i = 0; i < alignments; i++){
		sequence_array_add(max_length, max_name_length, tmp->sequences);
		tmp->alignments[i] = alignment_new(max_length, max_name_length);
		alignment_set_sequences(tmp->reference, sequence_array_get_sequence(i, tmp->sequences),tmp->alignments[i]);
	}
	
	return tmp;
}

void alingment_to_reference_destroy(AlignmentToReference ** atr){
	AlignmentToReference * ar = (*atr);
	int i;
	for(i = 0; i < ar->capacity; i++){
		alignment_destroy(&(ar->alignments[i]));
	}
	//peptide_sequence_destroy(&(atr->concensus_protein)); TODO FIX THIS LEAK
	//peptide_sequence_destroy(&(atr->reference_protein));
	free(ar->alignments);
	sequence_array_destroy(&(ar->sequences));
	free_sequence(&(ar->reference));
	free(ar);
	(*atr) = NULL;
}

Sequence * alignment_to_reference_get_reference(AlignmentToReference * atr){
#ifndef NO_BOUNDS_CHECK	//For debugging purposes, in "production", if coded properly, we can skip the validation
	if(atr == NULL){
		fprintf(stderr, "[alignment_to_reference_get_reference] The aligment to reference cant be NULL!");
		exit(-1);
	}
#endif 
	return atr->reference;
}


Sequence * alignment_to_reference_get_sequence(int pos, AlignmentToReference * atr){
#ifndef NO_BOUNDS_CHECK	//For debugging purposes, in "production", if coded properly, we can skip the validation
	if(atr == NULL){
		fprintf(stderr, "[alignment_to_reference_get_reference] The aligment to reference cant be NULL!");
		exit(-1);
	}
	
	if(pos >= atr->capacity){
		fprintf(stderr, "[alignment_to_reference_get_reference] Can't get a sequence from  possition %d on aligment to refernce with capacity %d!", pos, atr->capacity);
		exit(-1);
	}
#endif 
	return sequence_array_get_sequence(pos, atr->sequences);
}

void alignment_to_reference_clean( AlignmentToReference * atr){
	
#ifndef NO_BOUNDS_CHECK	//For debugging purposes, in "production", if coded properly, we can skip the validation
	if(atr == NULL){
		fprintf(stderr, "[alignment_to_reference_clean] The aligment to reference cant be NULL!");
		exit(-1);
	}
#endif		
	sequence_array_clean(atr->sequences);
	sequence_clean(atr->concensus);
	peptide_sequence_clean(atr->concensus_protein);
	int i;
	for(i = 0; i < atr->capacity; i++){
		alignment_clean(atr->alignments[i]);
	}
}

void alignment_to_reference_clean_reference( AlignmentToReference * atr){
	
#ifndef NO_BOUNDS_CHECK	//For debugging purposes, in "production", if coded properly, we can skip the validation
	if(atr == NULL){
		fprintf(stderr, "[alignment_to_reference_clean] The aligment to reference cant be NULL!");
		exit(-1);
	}
#endif		
	sequence_clean(atr->reference);
	sequence_clean(atr->concensus);
	peptide_sequence_clean(atr->concensus_protein);
	int i;
	for(i = 0; i < atr->capacity; i++){
		alignment_clean(atr->alignments[i]);
	}
}

void alignment_to_reference_mask_primer(AlignmentToReference *atr){
	int i;
	sequence_copy(atr->concensus, atr->reference);
	
	for(i = 0; i < atr->capacity; i++){
		if(atr->alignments[i]->score > atr->min_alignment_score){
			if(atr->alignments[i]->start_a >= atr->reference->length - atr->sequences->sequences[i]->length){
				sequence_mask(atr->alignments[i]->start_a,atr->reference->length,atr->concensus); 
			}
			
			if(atr->alignments[i]->end_a <= atr->sequences->sequences[i]->length){
				sequence_mask(0,atr->alignments[i]->end_a +1,atr->concensus); 
			}
		}
		
	}
}

void alignment_to_reference_align_all(AlignmentToReference *atr, ALIGNMENT_SCORE_DATATYPE(* score_rule)(char a, char b )){
	
	assert(atr != NULL);
	
	int i, j;
	Alignment * als[4];
	SequenceArray * sa = sequence_array_new(4);
	for(j = 0; j < 4; j++){
		als[j] = alignment_new(atr->max_lenght,atr->max_name_length); 
		sequence_array_add(atr->max_lenght,atr->max_name_length, sa);
	}
	
	ALIGNMENT_SCORE_DATATYPE max_score = 0;
	int max = 0;
	
	
	
	
	for(i = 0; i < atr->capacity; i++){
		max_score = 0;
		
		for(j = 0; j < 4; j++){
			
			sequence_copy(sequence_array_get_sequence(j, sa), alignment_to_reference_get_sequence(i, atr));
			alignment_set_sequences(alignment_to_reference_get_reference(atr),sequence_array_get_sequence(j, sa),  als[j]);
			switch(j){
				case 0:
					break;
					
				case 1:
					sequence_reverse_complement(sequence_array_get_sequence(j, sa));
					break;
					
				case 2:
					sequence_reverse(sequence_array_get_sequence(j, sa));
					break;
					
				case 3:
					sequence_complement(sequence_array_get_sequence(j, sa));
					break;
			};
			
			alignment_align(als[j], score_rule);
			if(als[j]->score > max_score){
				max_score = als[j]->score;
				max = j;
			}
		}
	//	printf("Max: %d ",als[max]->score);
	//	printf("\n%s··············", atr->alignments[i]->b_aligned->qual);
		alignment_copy(atr->alignments[i], als[max]);
	//	printf("%s\n", atr->alignments[i]->b_aligned->qual);
//		j= max;
		sequence_copy(sequence_array_get_sequence(max, sa), alignment_to_reference_get_sequence(i, atr));
		alignment_set_sequences(sequence_array_get_sequence(max, sa), alignment_to_reference_get_reference(atr), als[max]);
		
		switch(max){
		 
		 case 1:
				
		 break;
		 
		 case 2:
		 sequence_reverse_complement(alignment_to_reference_get_sequence(i, atr));
		 break;
		 
		 case 3:
		 sequence_reverse(alignment_to_reference_get_sequence(i, atr));
		 break;
		 
		 case 4:
		 sequence_complement(alignment_to_reference_get_sequence(i, atr));
		 break;
		 };
		
		
	/*	alignment_clean(atr->alignments[i]);
		alignment_align(atr->alignments[i], score_rule);
		
		 alignment_clean(atr->alignments[i]);
		 alignment_align(atr->alignments[i], score_rule);
		 alignment_align(atr->alignments[i], score_rule);
		 
		 if(atr->alignments[i]->score < atr->min_alignment_score){
		 sequence_reverse_complement(alignment_to_reference_get_sequence(i, atr));
		 alignment_clean(atr->alignments[i]);
		 alignment_align(atr->alignments[i], score_rule);
		 }
		 
		 if(atr->alignments[i]->score < atr->min_alignment_score){
		 sequence_reverse(alignment_to_reference_get_sequence(i, atr));
		 alignment_clean(atr->alignments[i]);
		 alignment_align(atr->alignments[i], score_rule);
		 }
		 
		 if(atr->alignments[i]->score < atr->min_alignment_score){
		 sequence_complement(alignment_to_reference_get_sequence(i, atr));
		 alignment_clean(atr->alignments[i]);
		 alignment_align(atr->alignments[i], score_rule);
		 }*/
		
		
	}
	
	
	for(j = 0; j < 4; j++){
		alignment_destroy(&als[j]);
	}
	
}

void alignment_to_reference_print(FILE * f, AlignmentToReference * atr){
#ifndef NO_BOUNDS_CHECK	//For debugging purposes, in "production", if coded properly, we can skip the validation
	if(atr == NULL){
		fprintf(stderr, "[alignment_to_reference_print] The aligment to reference cant be NULL!");
		exit(-1);
	}
#endif
	int i;
	for(i = 0; i < atr->capacity; i++){
		alignment_print_formated(f, atr->alignments[i]);
	}
}

void alignment_to_reference_print_aligned(FILE * f, AlignmentToReference * atr){
#ifndef NO_BOUNDS_CHECK	//For debugging purposes, in "production", if coded properly, we can skip the validation	
	if(atr == NULL){
		fprintf(stderr, "[alignment_to_reference_print_aligned] The aligment to reference cant be NULL!");
		exit(-1);
	}
#endif	
	int i;
	for(i = 0; i < atr->capacity; i++){
		alignment_print_formated(f, atr->alignments[i]);
	}
}



void alignment_to_reference_set_threadshodlds(int qual, int score, AlignmentToReference * atr){
#ifndef NO_BOUNDS_CHECK	//For debugging purposes, in "production", if coded properly, we can skip the validation
	if(atr == NULL){
		fprintf(stderr, "[alignment_to_reference_set_threadshodlds] The aligment to reference cant be NULL!");
		exit(-1);
	}
#endif
	atr->qual_threashold = qual;
	atr->min_alignment_score = score;
}


void alignment_to_reference_consense_by_quality(AlignmentToReference * atr){
#ifndef NO_BOUNDS_CHECK	//For debugging purposes, in "production", if coded properly, we can skip the validation
	if(atr == NULL){
		fprintf(stderr, "[alignment_to_reference_consense_by_quality] The aligment to reference cant be NULL!");
		exit(-1);
	}
	if(atr->reference == NULL){
		fprintf(stderr, "[alignment_to_reference_consense_by_quality] The reference cant be NULL!");
		exit(-1);
	}
#endif
	
	int i;//iterator for the lenght of the sequence
	int j;//iterator for each alignment
	int  local_index;
	int total_length = sequence_get_length(atr->reference);
	char best_base;
	char best_qual;
	char tmp_base;
	char tmp_qual, tmp_ref;
	
	
	int * offset = calloc(atr->capacity, sizeof(int)); //This fills the array with 0s, so it shall be ok to use the counters straight away. 
	int gaps = 0;
	for(i = 0; i < total_length; i++){
		best_base = sequence_get_base(i, atr->reference);
		best_qual = 0;
		//printf("|");
		for(j = 0; j < atr->capacity; j ++) {
			gaps = 0;
		//	printf(".");
			//TODO the actual code which decides the base to use. 
			if( atr->alignments[j]->score > atr->min_alignment_score && i >= atr->alignments[j]->start_a  && i < atr->alignments[j]->end_a ){
				local_index = i - atr->alignments[j]->start_a+offset[j];
				tmp_ref  = sequence_get_base(local_index, atr->alignments[j]->a_aligned);
				
				while(tmp_ref == '-'){
					offset[j]++;
					gaps++;
					local_index = i - atr->alignments[j]->start_a+offset[j];
					tmp_ref  = sequence_get_base(local_index, atr->alignments[j]->a_aligned);
			    }
				tmp_base = sequence_get_base(local_index, atr->alignments[j]->b_aligned);
				tmp_qual = sequence_get_qual(local_index, atr->alignments[j]->b_aligned) ;
				
				//printf("%d < %d\t",tmp_qual , best_qual);
				if( tmp_base != '-' &&  tmp_base != 'N' &&  tmp_qual > best_qual && tmp_qual > atr->qual_threashold ){
					best_base = tmp_base;
					best_qual = tmp_qual;
				}
				
			}
			
		}
		if(gaps > 0){
			best_qual = 0;
			
		}
		if(best_qual == 0){
			best_base = tolower(sequence_get_base(i, atr->reference));
			//best_base = '-';
		}
		sequence_add_base(best_base, best_qual, atr->concensus);	
	}
	free(offset);
}

int alignment_to_reference_count_correct_bases(AlignmentToReference * atr){
	int i, len = atr->concensus->length, count = 0;
	char * seq = atr->concensus->seq;
	for(i = 0; i < len; i++){
		if(isupper(seq[i]))
			count++;
	}
	return count;
}

void alignment_to_reference_quality_filter_window(AlignmentToReference * atr){
	
	Sequence * reference = atr->reference;
	Sequence * concensus = atr->concensus;
	
	int i, j, min_qual, length = sequence_get_length(reference), low, up;
	
	for(i = 0; i < length; i ++){
		min_qual = 100;
		low = i - atr->qual_window/2;
		up = i + atr->qual_window/2 ;
		for(j = low ; j < up && j < length && min_qual >= atr->qual_threashold; j++){
			if(j > 0 && concensus->qual[j] < min_qual){
				min_qual = concensus->qual[j];
			}
		}
		if(min_qual < atr->qual_threashold){
			concensus->seq[i]  = tolower(reference->seq[i]);
			//concensus->qual[i] = atr->qual_threashold;
		}
		
		/**
		 * The following block is to have a quick and dirty way to print to stderr the counts of the variations...
		 * we want to put it in it's own function.  
		 */
		if( toupper(concensus->seq[i]) != toupper(reference->seq[i])){
			fprintf(stderr, "%s: %c%d->%c\n", concensus->name, reference->seq[i], i+1, concensus->seq[i] ); 
		}	
		
	}
}



void alignment_to_reference_print_stats_header(FILE * f){
	
	fprintf(f, "REF\tID\tFROM\tPOSITION\tTO\tREF_LENGTH\tTYPE\n");	
}
/**
 * It requires the header!
 * For convenience, it returns true if it found a mutation
 * WARNING! it also prints if the sequences are synonimous
 */ 
//TODO get the string to print it on the header too. Use a single line for each read, with all the information. 
//TODO mark sequence. 
boolean alignment_to_reference_print_stats(FILE * f,  AlignmentToReference * atr){
	Sequence * reference = atr->reference;
	Sequence * concensus = atr->concensus;
	int i, length = sequence_get_length(reference);
	Peptide *r = atr->reference_protein->peptide;
	Peptide *c = atr->concensus_protein->peptide;
	int p_len = atr->reference_protein->length;
	boolean found = false;
	char name[100];
	sscanf(concensus->name,"%s.", &name[0]);
	int correct = alignment_to_reference_count_correct_bases(atr);
	
	if(correct > atr->min_correct_bases){
		for(i = 0; i < length; i ++){
			
			if( toupper(concensus->seq[i]) != toupper(reference->seq[i])){
				fprintf(f,"'%s'\t'%s'\t%c\t%d\t%c\t%d\tBASE\n",reference->name ,&name[0], reference->seq[i], i+1, concensus->seq[i], reference->length );
			//	fprintf(stdout,"'%s'\t'%s'\t%c\t%d\t%c\t%d\tBASE\n",reference->name ,&name[0], reference->seq[i], i+1, concensus->seq[i], reference->length );
				//TODO Review the consistency of the columns (reference-clone)
				found = true; 
				if(r[i	/3] != c[i/3]){
					fprintf(f,"'%s'\t'%s'\t%c\t%d\t%c\t%d\tNO_SYN\n",reference->name ,&name[0], peptide_to_char(r[i/3]), (i/3)+1, peptide_to_char(c[i/3]), p_len );
									}else{	
					fprintf(f,"'%s'\t'%s'\t%c\t%d\t%c\t%d\tSYN\n",reference->name ,&name[0], peptide_to_char(r[i/3]), (i/3)+1, peptide_to_char(c[i/3]), p_len );
				}
			}
		}
	}
	fprintf(f,"'%s'\t'%s'\t'-'\t%d\t'-'\t%d\t%s\n",reference->name ,&name[0], correct, atr->min_correct_bases, correct > atr->min_correct_bases?"VALID":"INVALID");
	
	return found;
}

//Functions for the masked alginment. 

MaskAlignment * mask_alignment_new(Sequence * reference){
	
	MaskAlignment * ma = calloc(1, sizeof(MaskAlignment));
	ma->unmatched = sequence_new(reference->max_length, reference->max_name_length, 0);
	ma->reference = reference;
	return ma;
}

void mask_alignment_set_sequence(Sequence * seq, MaskAlignment * ma){
	ma->sequence = seq;
	sequence_clean(ma->unmatched);
	
}

void mask_alignment_destroy(MaskAlignment ** ma){
	free_sequence(&(*ma)->unmatched);
	free(*ma);
	*ma = NULL;
}

int mask_alignment_longest_match(MaskAlignment * ma){
#ifndef NO_BOUNDS_CHECK	//For debugging purposes, in "production", if coded properly, we can skip the validation
	if(ma == NULL){
		fprintf(stderr, "[mask_alignment_longest_match] The MaskAlignment cant be NULL!");
		exit(-1);
	}
	if(ma->reference == NULL){
		fprintf(stderr, "[mask_alignment_longest_match] The reference cant be NULL!");
		exit(-1);
	}
	
	if(ma->sequence == NULL){
		fprintf(stderr, "[mask_alignment_longest_match] The sequence cant be NULL!");
		exit(-1);
	}
#endif
	
	/*if(!ma->sequence->upper_case)
	 sequence_to_upper_case(ma->sequence);
	 if(!ma->reference->upper_case)
	 sequence_to_upper_case(ma->sequence);*/
	int i;
	char * r  = ma->reference->seq;
	char * s  = ma->sequence->seq;
	char * al = ma->unmatched->seq;
	int len   =  ma->reference->length;
	int matches = 0;
	int mismatches = 0;
	int longest = 0;
	int constant_match=0;
	int local_mismatch = 0;
	for(i = 0; i <len; i++){
		if(r[i] != 'N'){
			if(base_is_valid(s[i], r[i] )){
				al[i] = '*';
				matches++;
				constant_match++;
				if(constant_match > ma->search_window)
					longest = i;
				local_mismatch = 0;
			}else{
				al[i] = s[i];
				local_mismatch++;
				mismatches++;
				if(local_mismatch > ma->max_mismatches){
					constant_match = 0;
				}
			}
		}else{
			//matches++;
			constant_match++;
			if(constant_match > ma->search_window)
				longest = i;
			al[i] = '-';
		}
	}
	al[i] = '\0';
	ma->matches = matches;
	ma->mismatches = mismatches;
	ma->longest = longest;
	return longest;
	
}

int mask_alignment_score(MaskAlignment * ma){
	return  ma->matches - (ma->mismatches )- 2*sequence_count_gaps(ma->sequence, ma->sequence->length) - (ma->fix_count);
	
}

int mask_alignment_compare( void * a,  void * b){
	MaskAlignment * ma_a = (MaskAlignment*)a;
	MaskAlignment * ma_b = (MaskAlignment*)b;
	
	int diff = ma_a->sequence->length - ma_b->sequence->length;
	/*if (diff != 0) {
	 return diff;
	 }*/
	
	int score_a = mask_alignment_score(ma_a);
	int score_b = mask_alignment_score(ma_b);
	
	diff = score_b - score_a;
	//fprintf(stderr, "Compared(%d): %s, %s\n", diff,ma_a->sequence->name, ma_b->sequence->name);
	if(diff == 0){
		diff = strcasecmp(ma_a->sequence->seq, ma_b->sequence->seq);
		
	}
	if(diff == 0){
		diff = strcmp(ma_a->sequence->name, ma_b->sequence->name);
		
	}
	
	
	return diff;
}

int mask_alignment_simple_compare( void * a,  void * b){
	MaskAlignment * ma_a = (MaskAlignment*)a;
	MaskAlignment * ma_b = (MaskAlignment*)b;
	
	
	
	int diff = strcmp(ma_a->sequence->name, ma_b->sequence->name);
	
	
	
	
	return diff;
}



void mask_alignment_print_alignment_fasta(FILE * f, MaskAlignment * ma){
	//sequence_set_name(ma->sequence->name, ma->aligned);
#ifndef NO_BOUNDS_CHECK	//For debugging purposes, in "production", if coded properly, we can skip the validation
	if(ma == NULL){
		fprintf(stderr, "[mask_alignment_print_alignment_fasta] ma is NULL\n");
		exit(-1);
	}
	if(ma->unmatched == NULL){
		fprintf(stderr, "[mask_alignment_print_alignment_fasta] ma->unmatched is NULL\n");
		exit(-1);
	}
	if(ma->unmatched->name == NULL){
		fprintf(stderr, "[mask_alignment_print_alignment_fasta] ma->unmatched->name is NULL\n");
		exit(-1);
	}
	if(f == NULL){
		fprintf(stderr, "[mask_alignment_print_alignment_fasta] f is NULL\n");
		exit(-1);
	}
#endif
	fprintf(f, ">%s\tmatches=%d\tmismatches=%d\tlongest=%d\n", ma->unmatched->name, ma->matches, ma->mismatches, ma->longest);
	fprintf(f, "%s\n", ma->unmatched->seq);
	
}

void mask_alignment_print_stats_header(FILE * f){
	
	fprintf(f, "REF\tID\tLONGEST\tMATCHES\tMISMATCHES\tGAPS\tEDITS\tSCORE\n");	
}
/**
 * It requires the header!
 * 
 * 
 */ 
void mask_alignment_print_print_stats(FILE * f,  MaskAlignment * ma){
	fprintf(f, "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", ma->reference->name, ma->sequence->name, ma->longest, ma->matches, ma->mismatches, sequence_count_gaps(ma->sequence, ma->reference->length), ma->fix_count, mask_alignment_score(ma));	
}

void mask_alignment_align(MaskAlignment * ma){
	
	assert(ma!=NULL);
	assert(ma->reference != NULL);
	assert(ma->sequence != NULL);
	
	int seq_len = ma->reference->max_length;
	int name_len = ma->reference->max_name_length;
	
	Alignment * al = alignment_new(seq_len, name_len);
	alignment_set_sequences(ma->reference, ma->sequence, al);
	
	alignment_align(al, &default_nucleotide_score);
	
	sequence_set_name(ma->sequence->name, al->b_aligned);
	sequence_copy(ma->sequence,al->b_aligned);
	
	alignment_destroy(&al);
	
}

void mask_alignment_fix_homopolymers(MaskAlignment * ma){
	Sequence * seq = ma->sequence;
	int i, hl, hr, edits;
	edits = 0;
	for(i = 0; i<seq->length; i++){
		if(seq->seq[i]=='-'){
			hl = 0;
			hr = 0;
			if(i > 0){
				hl = sequence_count_homopolymer(false, i-1, seq);
			}
			
			hr = sequence_count_homopolymer(true, i+1, seq);
			
			if(hl > hr){
				seq->seq[i] = tolower(seq->seq[i-1]);
				edits++;
			}else if(hl < hr){
				seq->seq[i] = tolower(seq->seq[i+1]);
				edits ++;
			}else if(hl > 1){
				if(seq->seq[i+1] == seq->seq[i-1]){
					seq->seq[i] = tolower(seq->seq[i+1]);
					edits ++;
				}
			}else if(seq->seq[i+1] == seq->seq[i-1]){
				seq->seq[i] = tolower(seq->seq[i+1]);
				edits ++;
			}/*else{
			  seq->seq[i] = tolower(seq->seq[i+1]); 
			  edits++;
			  }*/
			
		}
	}
	ma->fix_count += edits;
}

static int shift_from_pos(int initial_pos, int last_to_search, int max_mismatches ,Sequence * seq, Sequence * ref){
	int i,j;
	int pos = initial_pos - 1;
	int missed_bases;
	char ref_char;
	char seq_char;
	do{
		pos++;//we are detecting insertions, so we shift  to the left the sequence.
		missed_bases = 0;
		
		for(i = pos, j = initial_pos ; i < seq->length && j < ref->length && i < last_to_search;i++, j++){ 
			
			seq_char = sequence_get_base(i, seq);
			ref_char = sequence_get_base(j, ref);
			
			
			if(ref_char != 'N' && seq_char != '-'){//Dont count gaps, other parts of the algorithm take care of them... 
				if(!base_is_valid( seq_char, ref_char)){
					missed_bases++;
				}
			}
		}
		
	}while( pos > 0 &&  missed_bases > max_mismatches );
	
	return pos - initial_pos;
}

static int back_shift_from_pos(int initial_pos, int last_to_search, int max_mismatches ,Sequence * seq, Sequence * ref){
	int i,j;
	int pos = initial_pos + 1;
	int missed_bases;
	char ref_char;
	char seq_char;
	do{
		pos--;//we are detecting deletions, so we shift  to the left the sequence.
		missed_bases = 0;
		
		for(i = pos, j = initial_pos ; i < seq->length && j < ref->length && i < last_to_search;i++, j++){ 
			
			seq_char = sequence_get_base(i, seq);
			ref_char = sequence_get_base(j, ref);
			
			
			if(ref_char != 'N' && seq_char != '-'){//Dont count gaps, other parts of the algorithm shall take care of them... we are
				if(!base_is_valid( seq_char, ref_char)){
					missed_bases++;
					
				}
			}
		}
		
	}while( pos > 0 &&  missed_bases > max_mismatches );
	
	return initial_pos - pos;
}

static int first_correct_base_from_pos(int initial_pos, int last_to_search, int max_mismatches ,Sequence * seq, Sequence * ref){
	int i = initial_pos + 1;
	int shift;
	do{
		shift = shift_from_pos( i,  last_to_search,  max_mismatches , seq, ref);
	}while(i++ < last_to_search && shift != 0);
	
	return i == seq->length?-1:i-1;
}

void mask_alignment_fix_local_del_in(MaskAlignment * ma){
	Sequence * seq = ma->sequence;
	Sequence * ref = ma->reference;
	int i, edits;
	edits = 0;
	int max_shifts = 3;
	int correct_from = 0;
	int shift = 0;

	int prev_homo, prev_anch;
	int next_homo;
	for(i = 0; i<seq->length && i < ref->length && edits < max_shifts; i++){
		
		//if(seq->seq[i] !='-' && ref->seq[i] !='N' && seq->seq[i] != ref->seq[i]){
		if(seq->seq[i] !='-' && ref->seq[i] !='N' && !base_is_valid( seq->seq[i], ref->seq[i])){
			
			correct_from = first_correct_base_from_pos(i, seq->length, max_shifts ,seq, ref);
			
			if( correct_from > i){
				
				//	fprintf(stdout, "seq:%s(pos:%d, edits:%d):correct_from %d\n",seq->name, i,edits,  correct_from);
				
				shift = shift_from_pos( i,  correct_from,  max_shifts , seq, ref);
				prev_anch = sequence_prev_anchor_base(i-1, seq);
				prev_homo = sequence_prev_hompoplymer(i, prev_anch, seq);
				if(prev_homo == -1){
					printf("Not prev homo\n");
					prev_homo = i - 1;
				}
				//		fprintf(stdout, "seq:%s(pos:%d, edits:%d):shift %d\n",seq->name, i,edits,  shift);
				if(shift == 1){	
					next_homo = correct_from ;
					while(seq->seq[next_homo] == seq->seq[next_homo+1]){
						next_homo++;
					}
					
					sequence_remove_base_up_to_limit(prev_homo, next_homo, seq);
					i = 0;
				}
			}
		}
	}
}

void mask_alignment_fix_local_in_del(MaskAlignment * ma){
	Sequence * seq = ma->sequence;
	Sequence * ref = ma->reference;
	int i, edits;
	edits = 0;
	int max_shifts = 3;
	int correct_from = 0;
	
	int shift = 0;
	
	int prev_homo, prev_anch;
	int next_homo;
	for(i = 0; i<seq->length && i < ref->length && edits < max_shifts; i++){
		
		if(seq->seq[i] !='-' && ref->seq[i] !='N' && !base_is_valid( seq->seq[i], ref->seq[i])){
			
			correct_from = first_correct_base_from_pos(i, seq->length, max_shifts ,seq, ref);
			
			if( correct_from > i){
				
				//	fprintf(stdout, "[mask_alignment_fix_local_in_del] seq:%s(pos:%d, edits:%d):correct_from %d\n",seq->name, i,edits,  correct_from);
				
				shift = back_shift_from_pos( i,  correct_from,  max_shifts , seq, ref);
				prev_anch = sequence_prev_anchor_base(i-1, seq);
				prev_homo = sequence_prev_hompoplymer(i, prev_anch, seq);
				if(prev_homo == -1){
					printf("Not prev homo\n");
					prev_homo = i - 1;
				}
				//	fprintf(stdout, "seq:%s(pos:%d, edits:%d):shift %d\n",seq->name, i,edits,  shift);
				if(shift == 1){	
					next_homo = correct_from ;
					while(seq->seq[next_homo] == seq->seq[next_homo+1]){
						next_homo++;
					}
					
					sequence_insert_base_up_to_limit('-',prev_homo, next_homo, seq);
					i = 0;
				}
			}
		}
	}
}

void mask_alignment_fix_deletions(MaskAlignment * ma){
	Sequence * seq = ma->sequence;
	Sequence * ref = ma->reference;
	int i, edits;
	edits = 0;
	int max_shifts = 100;
	int prev_homo;
	
	for(i = 0; i<seq->length && i < ref->length && edits < max_shifts; i++){
		if (seq->seq[i]=='-') {
			
			prev_homo = sequence_prev_hompoplymer(i, i-10, seq);
			sequence_remove_base(i, seq);
			//printf("inserting-del!!!\n");
			sequence_insert_base_up_to_limit(seq->seq[prev_homo], prev_homo, seq->length, seq);
			edits++;
			i = 0;
		}/*else if (ref->seq[i] != 'N' && seq->seq[i] != ref->seq[i]) {
		  prev_homo = sequence_prev_hompoplymer(i, i-10, seq);
		  //printf("inserting!!!\n");
		  sequence_insert_base_up_to_limit(seq->seq[prev_homo], prev_homo, seq->length, seq);
		  edits++;
		  i = 0;
		  }*/
	}
	
	ma->fix_count += edits;
}

void mask_alignment_fix_insertions(MaskAlignment * ma){
	Sequence * seq = ma->sequence;
	Sequence * ref = ma->reference;
	int i, edits;
	edits = 0;
	int max_shifts = 3;
	int prev_homo;
	int prev_anch;
	int shift, shift2;
	
	boolean shifted = false;
	for(i = 0; i<seq->length && i < ref->length && edits < max_shifts; i++){
		if(seq->seq[i] !='-' && ref->seq[i] !='N' && !base_is_valid( seq->seq[i], ref->seq[i])){
			shift = shift_from_pos(i, seq->length, max_shifts ,seq, ref);
			//fprintf(stdout, "seq:%s(pos:%d, edits:%d):shift %d\n",seq->name, i,edits,  shift);
			if(!shifted && sequence_count_homopolymer(false, i, seq)>1 && sequence_count_homopolymer(true, i+1, seq) >1 && shift == 0  ){
				if(base_is_valid( seq->seq[i+1], ref->seq[i])){//ONLY IF THE HOMOPOLYMER IS ON THE REFERENCE!
			 //	if(seq->seq[i+1]==ref->seq[i]){//ONLY IF THE HOMOPOLYMER IS ON THE REFERENCE!
					//	fprintf(stdout, "_____TOREF_1!seq:%s(pos:%d, edits:%d):shift %d\n",seq->name, i,edits,  shift);
					seq->seq[i] = ref->seq[i]; //TODO check what happens now that we have the IUAPC codes...
					//shifted = true;
					i = 0;
				}
			}else if(!shifted && sequence_count_homopolymer(false, i-1, seq)>1 && sequence_count_homopolymer(true, i, seq) >1 && shift == 0  ){
			if(base_is_valid(seq->seq[i-1],ref->seq[i])){
			//	if(seq->seq[i-1]==ref->seq[i]){
					//					fprintf(stdout, "_____TOREF_2!seq:%s(pos:%d, edits:%d):shift %d\n",seq->name, i,edits,  shift);
					seq->seq[i]=ref->seq[i];
					//	shifted = true;
					i = 0;
				}
			}else if(shift < max_shifts && shift > 0){
				if(sequence_count_homopolymer(false, i, seq) || sequence_count_homopolymer(false, i-1, seq)){
					sequence_remove_base(i-1, seq);
					edits++;
					i = 0;
				}else if(sequence_count_homopolymer(true, i, seq) || sequence_count_homopolymer(false, i+1, seq)){
					sequence_remove_base(i+1, seq);
					edits++;
					i = 0;
				}else {
					prev_anch = sequence_prev_anchor_base(i, seq);
					prev_anch = sequence_prev_anchor_base(prev_anch, seq);
					prev_homo = sequence_prev_hompoplymer(i, prev_anch, seq);	
					//	fprintf(stdout, "_____prevhomo******??seq:%s(pos:%d, edits:%d):shift %d\n",seq->name, i,edits,  shift);
					if(prev_homo > 0){
						shift2 = shift_from_pos(prev_homo, 90, max_shifts ,seq, ref);
						//		fprintf(stdout, "_____prevhomo_______seq:%s(pos:%d, edits:%d):shift %d\n",seq->name, i,edits,  shift);
						if(shift2 < shift){
							
							sequence_remove_base(prev_homo, seq);
						}
					}
				}
			}
		}
	}
	ma->fix_count += edits;
}

//Functions to find related sequences. 

/**
 * Calculates an score for a nump, a good score needs to have a high identy and a similar lenght to the parental sequence
 * The score is between 0 and 1, where 0 is not likely to be a nump and 1 is the sequence itself. Basically, it checks on how close the sequences are.  
 * 
 * Alingment * al the alignment
 * double identity_threshold how close the sequences have to be. 1 represents 100%
 * lenght_window, as multiplying factor, where .1 represents +- 10% of the sequence lenght
 * 
 */
double nump_score(Alignment * al, double identity_threshold, double  lenght_window){
	
	assert(al != NULL);
	assert(al->a != NULL);
	assert(al->b != NULL);
	double min = (double) sequence_get_length(al->a)  * (1-lenght_window);
	double max = (double)al->a->length * (1+lenght_window);
	if(al->length > max || al->length < min) {//if the sequence is not in the range, score as 0
		return 0;
	}
	
	double identity = alignment_identity(al);
	
	if (identity < identity_threshold) {//We want a high identity, to know that the sequence is parental. 
		return 0;
	}
	
	return identity;
}



/**
 * Used to search for numps, it is better scored if the sequences
 * are very close and the alignment is nearly total. 
 *
 */
int compare_numps(void *p1, void *p2){
	
	Alignment * al1 = (Alignment*) p1;
	Alignment * al2 = (Alignment*) p2;
	
	//	assert(al1->a == al2->a);//Validate that both "parent" sequences are the same. 
	/*if (al1->a == al2->a) {
	 return 0;
	 }*/
	
	if (al1->b == al2->b) {
		return 0;
	}
	
	double score1 = nump_score(al1, .8, .1);
	double score2 = nump_score(al2, 0.8, .1);
	
	/*	if(score1 == 0 || score2 == 0)
	 return 0;//This prevents the non nump to be added
	 */
	if (score1 > score2) {
		return -1;
	}
	else if (score1 < score2) {
		return 1;
	}else {
		return strcmp(al1->b->name, al2->b->name);
	}
	
}


double chimera_score(Alignment * al, double identity_threshold, double  min_lenght, double max_lenght){
	
	assert(al != NULL);
	assert(al->a != NULL);
	assert(al->b != NULL);
	double min = (double) sequence_get_length(al->a)  * min_lenght;
	double max = (double)al->a->length * max_lenght;
	if(al->longest_lenght > max || al->longest_lenght < min) {//if the sequence is not in the range, score as 0
		return 0;
	}
	/*
	 double identity = alignment_identity(al);
	 
	 if (identity < identity_threshold) {//We want a high identity, to know that the sequence is parental. 
	 return 0;
	 }*/
	
	return (double)al->longest_lenght/(double)al->a->length;
}

int compare_chimeras(void *p1, void *p2){
	
	Alignment * al1 = (Alignment*) p1;
	Alignment * al2 = (Alignment*) p2;
	
	//	assert(al1->a == al2->a);//Validate that both "parent" sequences are the same. 
	/*if (al1->a == al2->a) {
	 return 0;
	 }*/
	
	if (al1->b == al2->b) {
		return 0;
	}
	
	double score1 = chimera_score(al1, 0.95, 0.20, 0.90);
	double score2 = chimera_score(al2, 0.95, 0.20, 0.90);
	
	/*	if(score1 == 0 || score2 == 0)
	 return 0;//This prevents the non nump to be added
	 */
	if (score1 > score2) {
		return -1;
	}
	else if (score1 < score2) {
		return 1;
	}else {
		return strcmp(al1->b->name, al2->b->name);
	}
	
}


RelatedSequences * related_sequences_new(Sequence * seq){
	RelatedSequences * rs = calloc(1, sizeof(RelatedSequences));
	rs->sequence = seq;
	rs->numps = binary_tree_new(100,&compare_numps);
	rs->chimeras = binary_tree_new(100, &compare_chimeras);
	return rs;
}

void related_sequence_destroy(RelatedSequences ** rs){
	
}

void related_sequences_to_fasta(FILE * f, RelatedSequences * rs){
	
	sequence_print_fasta(f, rs->sequence);
	fprintf(f, ">NUMPS\nNUMPS_NUMPS_NUMPS_NUMPS_NUMPS_NUMPS_NUMPS_NUMPS_NUMPS_\n");
	
	void p(void * v){
		Alignment  * al = (Alignment *) v;
		Sequence * s = al->b;
		sequence_print_fasta(f, s);
	}
	
	binary_tree_sorted_walk(&p, rs->numps);
	
	
	
	fprintf(f, ">CHIMERAS\nCHIMERAS_CHIMERAS_CHIMERAS_CHIMERAS_CHIMERAS_CHIMERAS_CHIMERAS_CHIMERAS\n");
	void p2(void * v2){
		
		Alignment  * al2 = (Alignment *) v2;
		Sequence * s = al2->b;
		sequence_print_fasta_subseq(f, al2->start_longest_b, al2->end_longest_b, s);
		//			&&sequence_print_fasta(f, s);
		
		
	}
	
	binary_tree_sorted_walk(&p2, rs->chimeras);
	
	
	fprintf(f, ">NEXT\n____________________________________________________________________\n");
	
}
