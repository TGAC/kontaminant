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
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <nucleotide.h>
#include <binary_kmer.h>
#include <seq.h>
#include <peptide.h>

/*
typedef enum {
	Adenine = 0, Cytosine = 1, Guanine = 2, Thymine = 3, Undefined = 4,
//Never used for the binary representation of a kmer!.
} Nucleotide;
*/
static Peptide genetic_code[4][4][4];
static boolean generated = false;

static void print_genetic_code()
{
	int i, j, k;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			
			for (k = 0; k < 4; k++) {
				printf("%c%c%c->%c\t ",
				       binary_nucleotide_to_char(i),
				       binary_nucleotide_to_char(j),
				       binary_nucleotide_to_char(k),
				       peptide_to_char(genetic_code[i][j][k]));
				fflush(stdout);
			}
			printf("\t");
		}
		printf("\n");
	}

}

//Adenine = 0, Cytosine = 1, Guanine = 2, Thymine = 3, Undefined = 4,
static void init_code()
{
	generated = true;
	//TT
	genetic_code[Thymine][Thymine][Thymine] = Phenylalanine;
	genetic_code[Thymine][Thymine][Cytosine] = Phenylalanine;
	genetic_code[Thymine][Thymine][Adenine] = Leucine;
	genetic_code[Thymine][Thymine][Guanine] = Leucine;

	//CT
	genetic_code[Cytosine][Thymine][Thymine] = Leucine;
	genetic_code[Cytosine][Thymine][Cytosine] = Leucine;
	genetic_code[Cytosine][Thymine][Adenine] = Leucine;
	genetic_code[Cytosine][Thymine][Guanine] = Leucine;

	//AT
	genetic_code[Adenine][Thymine][Thymine] = Isoleucine;
	genetic_code[Adenine][Thymine][Cytosine] = Isoleucine;
	genetic_code[Adenine][Thymine][Adenine] = Isoleucine;
	genetic_code[Adenine][Thymine][Guanine] = Methionine;	//Also start codon

	//GT
	genetic_code[Guanine][Thymine][Thymine] = Valine;
	genetic_code[Guanine][Thymine][Cytosine] = Valine;
	genetic_code[Guanine][Thymine][Adenine] = Valine;
	genetic_code[Guanine][Thymine][Guanine] = Valine;

	//TC
	genetic_code[Thymine][Cytosine][Thymine] = Serine;
	genetic_code[Thymine][Cytosine][Cytosine] = Serine;
	genetic_code[Thymine][Cytosine][Adenine] = Serine;
	genetic_code[Thymine][Cytosine][Guanine] = Serine;

	//CC
	genetic_code[Cytosine][Cytosine][Thymine] = Proline;
	genetic_code[Cytosine][Cytosine][Cytosine] = Proline;
	genetic_code[Cytosine][Cytosine][Adenine] = Proline;
	genetic_code[Cytosine][Cytosine][Guanine] = Proline;

	//AC
	genetic_code[Adenine][Cytosine][Thymine] = Threonine;
	genetic_code[Adenine][Cytosine][Cytosine] = Threonine;
	genetic_code[Adenine][Cytosine][Adenine] = Threonine;
	genetic_code[Adenine][Cytosine][Guanine] = Threonine;

	//GC
	genetic_code[Guanine][Cytosine][Thymine] = Alanine;
	genetic_code[Guanine][Cytosine][Cytosine] = Alanine;
	genetic_code[Guanine][Cytosine][Adenine] = Alanine;
	genetic_code[Guanine][Cytosine][Guanine] = Alanine;

	//TA
	genetic_code[Thymine][Adenine][Thymine] = Tyrosine;
	genetic_code[Thymine][Adenine][Cytosine] = Tyrosine;
	genetic_code[Thymine][Adenine][Adenine] = Stop;
	genetic_code[Thymine][Adenine][Guanine] = Stop;

	//CA
	genetic_code[Cytosine][Adenine][Thymine] = Histidine;
	genetic_code[Cytosine][Adenine][Cytosine] = Histidine;
	genetic_code[Cytosine][Adenine][Adenine] = Glutamine;
	genetic_code[Cytosine][Adenine][Guanine] = Glutamine;

	//AA
	genetic_code[Adenine][Adenine][Thymine] = Asparagine;
	genetic_code[Adenine][Adenine][Cytosine] = Asparagine;
	genetic_code[Adenine][Adenine][Adenine] = Lysine;
	genetic_code[Adenine][Adenine][Guanine] = Lysine;

	//GA
	genetic_code[Guanine][Adenine][Thymine] = Aspartic_acid;
	genetic_code[Guanine][Adenine][Cytosine] = Aspartic_acid;
	genetic_code[Guanine][Adenine][Adenine] = Glutamic_acid;
	genetic_code[Guanine][Adenine][Guanine] = Glutamic_acid;

	//TG
	genetic_code[Thymine][Guanine][Thymine] = Cysteine;
	genetic_code[Thymine][Guanine][Cytosine] = Cysteine;
	genetic_code[Thymine][Guanine][Adenine] = Stop;
	genetic_code[Thymine][Guanine][Guanine] = Tryptophan;

	//CG
	genetic_code[Cytosine][Guanine][Thymine] = Arginine;
	genetic_code[Cytosine][Guanine][Cytosine] = Arginine;
	genetic_code[Cytosine][Guanine][Adenine] = Arginine;
	genetic_code[Cytosine][Guanine][Guanine] = Arginine;

	//AG
	genetic_code[Adenine][Guanine][Thymine] = Serine;
	genetic_code[Adenine][Guanine][Cytosine] = Serine;
	genetic_code[Adenine][Guanine][Adenine] = Arginine;
	genetic_code[Adenine][Guanine][Guanine] = Arginine;

	//GG
	genetic_code[Guanine][Guanine][Thymine] = Glycine;
	genetic_code[Guanine][Guanine][Cytosine] = Glycine;
	genetic_code[Guanine][Guanine][Adenine] = Glycine;
	genetic_code[Guanine][Guanine][Guanine] = Glycine;

	//print_genetic_code();

}

Peptide peptide_from_char_array(char a, char b, char c)
{
/*#ifndef NO_BOUNDS_CHECK
	if(triplet == NULL){
			fprintf(stderr, "[peptide_from_char_array] The triplet cant be NULL!");
			exit(-1);
	}
#*/
	if (!generated)
		init_code();
//      printf("|%c%c%c|",a,b,c);
	Nucleotide na = char_to_binary_nucleotide(a);
	Nucleotide nb = char_to_binary_nucleotide(b);
	Nucleotide nc = char_to_binary_nucleotide(c);
	if (na > 3 || nb > 3 || nc > 3) {
		printf("Invalid triplet!!!!!%c%c%c", a, b, c);
	}
//      fflush(stdout);
	return genetic_code[na][nb][nc];

}

/**
 * Translates the Sequnce seq into a PeptideSequence. The offset is set to start in certain possition the translation. 
 */

PeptideSequence *peptide_sequence_translate(Sequence * seq, int offset,
					    PeptideSequence * protein)
{

#ifndef NO_BOUNDS_CHECK
	if (seq == NULL) {
		fprintf(stderr,
			"[peptide_sequence_translate] The sequence cant be NULL!");
		exit(-1);
	}
	if (protein == NULL) {
		fprintf(stderr,
			"[peptide_sequence_translate] The protein cant be NULL (%s)!",
			seq->name);
		exit(-1);
	}

	if (seq->length > protein->max_length * 3) {
		fprintf(stderr,
			"[peptide_sequence_translate] The protein doesn't have enought space  ( capacity: %d)for the translation (required: %d) NULL!",
			protein->max_length, 1 + (seq->length / 3));
		exit(-1);
	}
#endif

	int i, len = sequence_get_length(seq), index;

	for (i = offset, index = 0; i < len; i += 3, index++) {
		protein->peptide[index] =
		    peptide_from_char_array(seq->seq[i], seq->seq[i + 1],
					    seq->seq[i + 2]);
		//fflush(stdout);
		//printf("%c%c%c(%d)->%c\n", seq->seq[i],seq->seq[i+1],seq->seq[i+2], i, peptide_to_char( protein->peptide[index]));
	}
	protein->length = index;
	return protein;
}

char peptide_to_char(Peptide p)
{
	switch (p) {
	case Phenylalanine:
		return 'F';
		break;
	case Leucine:
		return 'L';
		break;
	case Isoleucine:
		return 'I';
		break;
	case Methionine:
		return 'M';
		break;
	case Valine:
		return 'V';
		break;
	case Serine:
		return 'S';
		break;
	case Proline:
		return 'P';
		break;
	case Threonine:
		return 'T';
		break;
	case Alanine:
		return 'A';
		break;
	case Tyrosine:
		return 'Y';
		break;
	case Stop:
		return '*';
		break;
	case Histidine:
		return 'H';
		break;
	case Glutamine:
		return 'Q';
		break;
	case Asparagine:
		return 'N';
		break;
	case Lysine:
		return 'K';
		break;
	case Aspartic_acid:
		return 'D';
		break;
	case Glutamic_acid:
		return 'E';
		break;
	case Cysteine:
		return 'C';
		break;
	case Tryptophan:
		return 'W';
		break;
	case Arginine:
		return 'R';
		break;
	case Glycine:
		return 'G';
		break;
	default:
		fprintf(stderr,
			"[peptide_to_char] %d is not a valid aminoacid!\n", p);
		exit(-1);
		break;

	};

}

PeptideSequence *peptide_sequence_new(int max_length, int max_name_length)
{
	PeptideSequence *ps = calloc(1, sizeof(PeptideSequence));
	ps->name = calloc(max_name_length, sizeof(char));
	ps->peptide = calloc(max_length, sizeof(Peptide));
	ps->max_length = max_length;
	ps->max_name_length = max_name_length;
	return ps;
}

void peptide_sequence_clean(PeptideSequence * ps)
{
	ps->length = 0;

}

void peptide_sequence_destroy(PeptideSequence ** ps)
{
	free((*ps)->name);
	free((*ps)->peptide);
	free((*ps));
	(*ps) = NULL;
}

void peptide_sequence_print_fasta(FILE * f, PeptideSequence * protein)
{
	int i;
	fprintf(f, ">%s\n", protein->name);
	for (i = 0; i < protein->length; i++) {
		fprintf(f, "%c", peptide_to_char(protein->peptide[i]));

	}

	fprintf(f, "\n");
}

void peptide_sequence_iterator(void (*f) (Peptide, int),
			       PeptideSequence * protein)
{
	int i;
	for (i = 0; i < protein->length; i++) {
		f(protein->peptide[i], i);
	}
}

/**
 * 
 * 
 * 
 */
int peptide_sequence_count_differences(PeptideSequence * a, PeptideSequence * b)
{

#ifndef NO_BOUNDS_CHECK
	if (a == NULL) {
		fprintf(stderr,
			"[peptide_sequence_count_differences] Peptide sequence a is NULL");
	}

	if (b == NULL) {
		fprintf(stderr,
			"[peptide_sequence_count_differences] Peptide sequence b is NULL");
	}

	if (a->length != b->length) {
		fprintf(stderr,
			"[peptide_sequence_count_differences] a->length (%d) != b->length (%d) ",
			a->length, b->length);
	}
#endif

	int diff = 0;
	int i;
	int len = a->length;
	for (i = 0; i < len; i++) {
		if (a->peptide[i] == b->peptide[i]) {
			diff++;
		}
	}
	return diff;
}
