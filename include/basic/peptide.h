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
 
#ifndef PEPTIDE_H_
#define PEPTIDE_H_

typedef enum {
	Phenylalanine=1, 
	Leucine=2, 
	Isoleucine=3, 
	Methionine=4,
	Valine=5, 
	Serine=6, 
	Proline=7,
	Threonine=8, 
	Alanine=9, 
	Tyrosine=10, 
	Stop=0, 
	Histidine=11, 
	Glutamine=12, 
	Asparagine=13, 
	Lysine=14, 
	Aspartic_acid=15, 
	Glutamic_acid=16, 
	Cysteine=17, 
	Tryptophan=18, 
	Arginine=19, 
	Glycine=20

} Peptide;

typedef struct {
	Peptide * peptide;
	char * name; 
	int max_name_length;
	int length;
	int max_length;
}PeptideSequence;

Peptide peptide_from_char_array(char a, char b, char c);

char peptide_to_char(Peptide p);

PeptideSequence * peptide_sequence_new();

void peptide_sequence_clean(PeptideSequence * ps);

void  peptide_sequence_destroy(PeptideSequence ** ps);

PeptideSequence * peptide_sequence_translate(Sequence * seq, int offset, PeptideSequence * protein);

void peptide_sequence_print_fasta(FILE * f, PeptideSequence * protein);

void peptide_sequence_iterator(void(*f)(Peptide, int), PeptideSequence * protein);

int peptide_sequence_count_differences(PeptideSequence * a, PeptideSequence *b);


#endif
