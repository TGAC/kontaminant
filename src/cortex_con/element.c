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
/*
 element.c -- implements the nodes of the dBruijn graph
 */

#include <stdlib.h>
#include <global.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#ifdef THREADS
#include <pthread.h>
#endif
#include <string.h>

#include <global.h>
#include <flags.h>
#include <nucleotide.h>
#include <binary_kmer.h>
#include <element.h>
//If this is just for the flags, maybe we want to have the Path-specific node actions in the path
#include <path.h> 
#ifdef ENABLE_READ_PAIR
#include <binary_tree.h>
#include <dB_graph.h>
#include <read_pair.h>
#endif

long long int visited_count = 0;

void element_assign(Element * e1, Element * e2)
{
	int i;

	binary_kmer_assignment_operator((*e1).kmer, (*e2).kmer);
	e1->flags = e2->flags;

	for (i = 0; i < NUMBER_OF_COLOURS; i++) {
		e1->edges[i] = e2->edges[i];
		e1->coverage[i] = e2->coverage[i];
#ifdef INCLUDE_QUALITY_SCORES
		e1->quality_string_arrays[i].number_of_strings =
		    e2->quality_string_arrays[i].number_of_strings;
		e1->quality_string_arrays[i].limit =
		    e2->quality_string_arrays[i].limit;
		// WARNING: Probably need to change this to copy the quality strings array
		e1->quality_string_arrays[i].quality_strings =
		    e2->quality_string_arrays[i].quality_strings;
#endif
	}
}

boolean element_is_key(Key key, Element e, short kmer_size)
{
	if (key == NULL) {
		printf
		    ("Do not call element_is_key wth a NULL pointer. Exiting\n");
		exit(1);
	}

	return binary_kmer_comparison_operator(*key, e.kmer);
}

#ifdef INCLUDE_QUALITY_SCORES
void element_preallocate_quality_strings(Element * e, int c, int n)
{
	QualityStringArray *qa = &e->quality_string_arrays[c];
    
    if (n > 0) {
        qa->quality_strings = calloc(n, sizeof(QualityString));
        if (qa->quality_strings == NULL) {
            printf("Can't allocate memory for quality strings.\n");
            exit(1);
        }
    }
    
	qa->number_of_strings = 0;
	qa->limit = n;
}

void element_add_quality_string(Element * e, short c, char *q)
{
	QualityStringArray *qa = &e->quality_string_arrays[c];

    // Sanity checking
    assert(e != NULL);
    assert(q != NULL);
    assert(c >= 0);
    assert(qa->number_of_strings >= 0);    
    
	if (qa->number_of_strings == qa->limit) {
		qa->limit += 1;
		if (qa->limit == 1) {
			qa->quality_strings = malloc(qa->limit * sizeof(QualityString));
		} else {
			qa->quality_strings = realloc(qa->quality_strings, qa->limit * sizeof(QualityString));
		}
		if (qa->quality_strings == NULL) {
			printf("Can't allocate memory for quality strings.\n");
			exit(1);
		}
	}
    
	qa->quality_strings[qa->number_of_strings].quality = malloc(strlen(q) + 1);
	if (qa->quality_strings[qa->number_of_strings].quality == NULL) {
		printf("Can't allocate memory for quality strings.\n");
		exit(1);
	}
	strcpy(qa->quality_strings[qa->number_of_strings].quality, q);
	qa->number_of_strings++;
}
#endif

/* No longer required?
boolean element_smaller(Element e1, Element e2) {
	//return e1.edges < e2.edges;
	return db_node_get_edges(e1) < db_node_get_edges(e2);
}
*/

//TODO - make API safer - this gets contents of hash table, not  copy
BinaryKmer *element_get_kmer(Element * e)
{
	return &(e->kmer);
}

//int element_get_coverage(Element * e) {
//      return e->coverage;
//}

int element_get_coverage_all_colours(Element * e)
{
	int coverage = 0;
#ifndef ENABLE_BUBBLEPARSE
	coverage = e->coverage[0];
#else
	int c;
	for (c = 0; c < NUMBER_OF_COLOURS; c++) {
		coverage = coverage + e->coverage[c];
	}
#endif
	return coverage;
}

int element_get_coverage_by_colour(Element * e, short colour)
{
	if ((colour >= 0) && (colour < NUMBER_OF_COLOURS))
		return e->coverage[colour];
	else {
		fprintf(stderr,
			"[element_get_coverage_by_colour] Colour is out of allowable range.");
		exit(-1);
	}
}

int element_update_coverage(Element * e, short colour, int update)
{
	if ((colour >= 0) && (colour < NUMBER_OF_COLOURS)) {
		//printf("[element_update_coverage] Update element %x colour %d by %d Total %d\n", e, colour, update, element_get_coverage_all_colours(e));
		e->coverage[colour] += update;
		return e->coverage[colour];
	} else {
		fprintf(stderr,
			"[element_update_coverage] Colour is out of allowable range.");
		exit(-1);
	}
}

Key element_get_key(BinaryKmer * kmer, short kmer_size, Key preallocated_key)
{

	BinaryKmer local_rev_kmer;
	binary_kmer_initialise_to_zero(&local_rev_kmer);

	binary_kmer_reverse_complement(kmer, kmer_size, &local_rev_kmer);

	if (binary_kmer_less_than(local_rev_kmer, *kmer, kmer_size)) {
		binary_kmer_assignment_operator(*((BinaryKmer *)
						  preallocated_key),
						local_rev_kmer);
	} else {
		binary_kmer_assignment_operator(*((BinaryKmer *)
						  preallocated_key), *kmer);
	}

	return preallocated_key;

}

void element_initialise(Element * e, BinaryKmer * kmer, short kmer_size)
{

	//TODO - add check that the kmer passed in really is consistent with kmer_size
	int i;
	BinaryKmer tmp_kmer;
	binary_kmer_initialise_to_zero(&tmp_kmer);
	binary_kmer_assignment_operator(e->kmer,
					*(element_get_key
					  (kmer, kmer_size, &tmp_kmer)));

	e->flags = ASSIGNED;

	for (i = 0; i < NUMBER_OF_COLOURS; i++) {
		e->edges[i] = 0;
		e->coverage[i] = 0;
#ifdef INCLUDE_QUALITY_SCORES
		e->quality_string_arrays[i].number_of_strings = 0;
		e->quality_string_arrays[i].limit = 0;
		e->quality_string_arrays[i].quality_strings = NULL;
#endif
	}
}

Orientation db_node_get_orientation(BinaryKmer * k, dBNode * e, short kmer_size)
{
	if (binary_kmer_comparison_operator(e->kmer, *k) == true) {
		return forward;
	}

	BinaryKmer tmp_kmer;

	if (binary_kmer_comparison_operator(e->kmer,
					    *(binary_kmer_reverse_complement
					      (k, kmer_size,
					       &tmp_kmer))) == true) {
		return reverse;
	}

	printf
	    ("programming error - you have called  db_node_get_orientation with a kmer that is neither equal to the kmer in this node, nor its rev comp\n");
	char tmpseq1[kmer_size];
	char tmpseq2[kmer_size];
	printf("Arg 1 Kmer is %s and Arg 2 node kmer is %s\n",
	       binary_kmer_to_seq(k, kmer_size, tmpseq1),
	       binary_kmer_to_seq(&(e->kmer), kmer_size, tmpseq2));
	exit(1);

}

char binary_nucleotide_to_edge(Nucleotide base){
    return 1 << base;
}
//add one edge to element -- basically sets a bit in the edges char
void
db_node_add_labeled_edge(dBNode * e, Orientation o, short colour,
			 Nucleotide base)
{

	//set edge
	char edge = binary_nucleotide_to_edge(base);	// A (0) -> 0001, C (1) -> 0010, G (2) -> 0100, T (3) -> 1000

	if (o == reverse) {
		edge <<= 4;	//move to next nibble
	}
	//update edges
	e->edges[colour] |= edge;
}

/**
 * Given a pair of nodes and the node of the first one, get the orientation inw
 */

/*Orientation db_node_get_reverse_edge_between(dBNode * src, dBNode * tgt,
 Orientation src_o, short kmer_size) {
 BinaryKmer src_k, tgt_k, tmp_k;
 Nucleotide last_s, first_s;

 binary_kmer_assignment_operator(src_k, src->kmer);
 binary_kmer_assignment_operator(tgt_k, tgt->kmer);

 if (DEBUG) {
 char tmp1[kmer_size];
 char tmp2[kmer_size];
 //fprintf(stderr, "Comparing %s and %s \n", binary_kmer_to_seq(src_k,
 //	kmer_size, tmp1), binary_kmer_to_seq(tgt_k, kmer_size, tmp2));
 fprintf(stdout, "Comparing %s and %s \n", binary_kmer_to_seq(src_k,
 kmer_size, tmp1), binary_kmer_to_seq(tgt_k, kmer_size, tmp2));
 }

 if (src_o == reverse) {
 binary_kmer_assignment_operator(src_k,
 *(binary_kmer_reverse_complement(&src_k, kmer_size, &tmp_k)));
 }
 if (DEBUG) {
 char tmp1[kmer_size];
 char tmp2[kmer_size];
 //fprintf(stderr, "Oriented %s : %s \n", binary_kmer_to_seq(src_k,
 //kmer_size, tmp1), binary_kmer_to_seq(tgt_k, kmer_size, tmp2));
 fprintf(stdout, "Oriented %s : %s \n", binary_kmer_to_seq(src_k,
 kmer_size, tmp1), binary_kmer_to_seq(tgt_k, kmer_size, tmp2));
 }
 //last_s = binary_kmer_get_last_nucleotide(src_k);
 last_s = binary_kmer_get_last_nucleotide(src_k);
 first_s = binary_kmer_get_first_nucleotide(src_k, kmer_size);

 binary_kmer_assignment_operator(tmp_k, tgt->kmer);
 boolean opp = true;
 if (src_o == reverse) {
 binary_kmer_left_shift(tmp_k, 2, kmer_size);
 binary_kmer_modify_base(tmp_k, last_s, kmer_size, 0);
 if (DEBUG) {
 //	fprintf(stderr, "1");
 fprintf(stdout, "1");
 }
 opp = false;

 } else {
 if (DEBUG) {
 printf("2");
 //	fprintf(stderr, "2");
 }

 binary_kmer_right_shift(tmp_k, 2);
 binary_kmer_modify_base(tmp_k, first_s, kmer_size, kmer_size - 1);
 opp = true;
 }

 if (DEBUG) {
 char tmp1[kmer_size];
 char tmp2[kmer_size];
 printf("testing the forward! %s and %s \n", binary_kmer_to_seq(src_k,
 kmer_size, tmp1), binary_kmer_to_seq(tmp_k, kmer_size, tmp2));
 //fprintf(stderr, "testing the forward! %s and %s \n",
 //		binary_kmer_to_seq(src_k, kmer_size, tmp1), binary_kmer_to_seq(
 //			tmp_k, kmer_size, tmp2));
 }

 if (binary_kmer_comparison_operator(src_k, tmp_k)) {
 //return db_node_get_orientation(tmp_k, tgt, kmer_size);
 return forward;
 }

 binary_kmer_assignment_operator(tgt_k, tgt->kmer);
 binary_kmer_assignment_operator(tgt_k, *(binary_kmer_reverse_complement(
 &tgt_k, kmer_size, &tmp_k)));
 binary_kmer_assignment_operator(tmp_k, tgt_k);
 if (src_o == reverse) {
 binary_kmer_right_shift(tmp_k, 2);
 binary_kmer_modify_base(tmp_k, first_s, kmer_size, kmer_size - 1);
 if (DEBUG)
 printf("3");
 //	fprintf(stderr, "3");
 opp = false;
 } else {
 opp = true;
 if (DEBUG) {
 printf("4 : %d ", first_s);
 //		fprintf(stderr, "4 : %d ", first_s);
 }
 binary_kmer_left_shift(tmp_k, 2, kmer_size);
 binary_kmer_modify_base(tmp_k, last_s, kmer_size, 0);

 }

 if (DEBUG) {
 char tmp1[kmer_size];
 char tmp2[kmer_size];
 printf("testing the reverse! %s and %s \n", binary_kmer_to_seq(src_k,
 kmer_size, tmp1), binary_kmer_to_seq(tmp_k, kmer_size, tmp2));
 //fprintf(stderr, "testing the reverse! %s and %s \n",
 //			binary_kmer_to_seq(src_k, kmer_size, tmp1), binary_kmer_to_seq(
 //					tmp_k, kmer_size, tmp2));
 }
 if ((binary_kmer_comparison_operator(src_k, tmp_k)))
 //	return db_node_get_orientation(tmp_k, tgt, kmer_size);
 return  reverse;

 printf(
 "Querying for the connection between two kmers which can not be connected!\n");
 //exit - 1;
 return forward;
 }*/

//adding an edge between two nodes implies adding two labeled edges (one in each direction)
//be aware that in the case of self-loops in palindromes the two labeled edges collapse in one
//the orientation refers which k-mer should be linked (it could the key or its reverse)


boolean db_node_add_edge(dBNode * src_e, dBNode * tgt_e, Orientation src_o, Orientation tgt_o, short kmer_size, short colour){

	BinaryKmer src_k, tgt_k, tmp_kmer;

	/*if(DEBUG){
	   printf("Trying to assign %p and %p ", src_e, tgt_e);
	   } */
	binary_kmer_assignment_operator(src_k, src_e->kmer);
	binary_kmer_assignment_operator(tgt_k, tgt_e->kmer);

	if (src_o == reverse) {
		binary_kmer_assignment_operator(src_k,
						*(binary_kmer_reverse_complement
						  (&src_k, kmer_size,
						   &tmp_kmer)));
	}

	if (tgt_o == reverse) {
		binary_kmer_assignment_operator(tgt_k,
						*(binary_kmer_reverse_complement
						  (&tgt_k, kmer_size,
						   &tmp_kmer)));
	}

	db_node_add_labeled_edge(src_e, src_o, colour,
				 binary_kmer_get_last_nucleotide(&tgt_k));

	db_node_add_labeled_edge(tgt_e, opposite_orientation(tgt_o), colour,
				 binary_kmer_get_last_nucleotide
				 (binary_kmer_reverse_complement
				  (&src_k, kmer_size, &tmp_kmer)));

	return true;
}

// Check if edge exists (defaults to colour 0)
boolean
db_node_edge_exist(dBNode * element, Nucleotide base, Orientation orientation)
{
	char edge = db_node_get_edges(element);

	if (orientation == reverse) {
		edge >>= 4;
	}

	edge >>= base;

	edge &= 1;

	if (edge == 1) {
		return true;
	} else {
		return false;
	}
}

// Check if edge exists, regardless of colour
boolean
db_node_edge_exist_any_colour(dBNode * element, Nucleotide base,
			      Orientation orientation)
{
	char edge = db_node_get_edges_all_colours(element);

	if (orientation == reverse) {
		edge >>= 4;
	}

	edge >>= base;

	edge &= 1;

	if (edge == 1) {
		return true;
	} else {
		return false;
	}
}

boolean db_node_is_supernode_end(dBNode * element, Orientation orientation)
{
	char edges = db_node_get_edges(element);

	if (orientation == reverse) {
		//shift along so the 4 most significant bits become the 4 least - we've passed out argument by copy so not altering the original
		edges >>= 4;
	}

	edges &= 15;		// AND with 00001111 so that we only look at the 4 least significant bits

	//is supernode end EITHER if it has 0 edges out, or >1.
	return ((edges == 0) || ((edges != 1) && (edges != 2) && (edges != 4)
				 && (edges != 8)));
}

Orientation opposite_orientation(Orientation o)
{
   // assert(o==forward || o == reverse); //TODO: temporary assert, to expensive to have it all the time. 
    
	return o ^ 1;

}

// Reset an edge in a single colour
void db_node_reset_edge(dBNode * node, Orientation orientation, short colour, Nucleotide nucleotide)
{

	char edge = 1 << nucleotide;

	if (orientation == reverse) {
		edge <<= 4;
	}
	//toggle 1->0 0->1

	edge ^= (unsigned char)0xFF;	//xor with all 1's, ie 00010000 -> 11101111

	node->edges[colour] &= edge;	//reset one edge
}

// Reset an edge in all colours
void db_node_reset_edge_all_colours(dBNode * node, Orientation orientation, Nucleotide nucleotide)
{
	int c;
	for (c = 0; c < NUMBER_OF_COLOURS; c++) {
		db_node_reset_edge(node, orientation, c, nucleotide);
	}
}

//set all edges
void db_node_set_edges(dBNode * node, short colour, Edges edges)
{
	node->edges[colour] |= edges;
}

/* No longer required?
boolean db_node_edges_reset(dBNode * node) {
	db_node_reset_edges(node);
	return true;
}
*/

// Reset all edges for given colour
void db_node_reset_edges(dBNode * node, short colour)
{
	node->edges[colour] = 0;
}

// Reset all edges for all colours
void db_node_reset_edges_all_colours(dBNode * node)
{
	int c;
	for (c = 0; c < NUMBER_OF_COLOURS; c++) {
		db_node_reset_edges(node, c);
	}
}

// Defaults to colour 0
boolean db_node_has_precisely_one_edge(dBNode * node, Orientation orientation, Nucleotide * nucleotide)
{	
	assert(orientation != undefined);
    assert(node != NULL);
	Nucleotide n = Undefined;
	Edges edges = db_node_get_edges(node);
	short edges_count = 0;

	if (orientation == reverse) {
		edges >>= 4;
	}

	for (n = 0; n < 4; n++) {

		if ((edges & 1) == 1) {
			*nucleotide = n;
			edges_count++;

		}

		edges >>= 1;
	}

	return (edges_count == 1);
}

// Colourblind version
boolean db_node_has_precisely_one_edge_all_colours(dBNode * node, Orientation orientation, Nucleotide * nucleotide)
{

	Nucleotide n;
	Edges edges = db_node_get_edges_all_colours(node);
	short edges_count = 0;

	if (orientation == reverse) {
		edges >>= 4;
	}

	for (n = 0; n < 4; n++) {

		if ((edges & 1) == 1) {
			*nucleotide = n;
			edges_count++;
		}

		edges >>= 1;
	}

	return (edges_count == 1);
}

boolean db_node_has_unvisited_edge(dBNode * node, Orientation orientation, Nucleotide * nucleotide)
{
	short edges_count = 0;
#ifndef SHORT_FLAGS
	Nucleotide n;
	Edges edges = db_node_get_edges(node);
	Flags visited = node->flags & EDGES_FLAGS;

	

	if (orientation == reverse) {
		edges >>= 4;
		visited = (node->flags & EDGES_REV_FLAGS) >> 4;
	}
	printf("Checking visited: %x > %x\n", edges, visited);

	for (n = 0; n < 4 && edges_count == 0; n++) {
		if (((edges & 1) == 1 && (visited & 1) == 0)	/*|| ((edges & 1) == 1
								   && db_node_check_for_flag(node, PLAIN_NODE)) */ ) {
			//Probably make a function for this cryptic line.
			*nucleotide = n;
			edges_count++;
		}
		edges >>= 1;
		visited >>= 1;
	}
#else
    // TODO: db_node_has_unvisited_edge for short flags needs to be written.
    printf("ERROR: db_node_has_unvisited_edge for SHORT_FLAGS not yet written!");
    exit(1);
    //#warning Dont db_node_has_unvisited_edge, the method not relaying on the flags is not written yet 
#endif
	return (edges_count > 1);

}

//a conflict - bifurcation
boolean
db_node_has_precisely_two_edges(dBNode * node, Orientation orientation,
				Nucleotide * nucleotide1,
				Nucleotide * nucleotide2)
{

	Nucleotide n;
	Edges edges = db_node_get_edges(node);
	short edges_count = 0;

	if (orientation == reverse) {
		edges >>= 4;
	}

	for (n = 0; n < 4; n++) {

		if ((edges & 1) == 1) {
			if (edges_count == 0) {
				*nucleotide1 = n;
			}

			if (edges_count == 1) {
				*nucleotide2 = n;
			}
			edges_count++;
		}

		edges >>= 1;
	}

	return (edges_count == 2);
}

// The 'colourblind' version defaults to colour 0
int db_node_edges_count(dBNode * node, Orientation orientation)
{
	return db_node_edges_count_by_colour(node, orientation, 0);
}

// Edge counting by colour
int
db_node_edges_count_by_colour(dBNode * node, Orientation orientation,
			      short colour)
{
	Edges edges = db_node_get_edges_by_colour(node, colour);
	int count = 0;

	if (orientation == reverse) {
		edges >>= 4;
	}

	int n;
	for (n = 0; n < 4; n++) {
		if ((edges & 1) == 1) {
			count++;
		}
		edges >>= 1;
	}

	return count;
}

// Count all edges for all colours
int db_node_edges_count_all_colours(dBNode * node, Orientation orientation)
{
	Edges edges = 0;
	int count = 0;
	int i;

	// OR together edges for different colours
	for (i = 0; i < NUMBER_OF_COLOURS; i++)
		edges |= db_node_get_edges_by_colour(node, i);

	if (orientation == reverse) {
		edges >>= 4;
	}

	int n;
	for (n = 0; n < 4; n++) {
		if ((edges & 1) == 1) {
			count++;
		}
		edges >>= 1;
	}

	return count;
}

// Count number of different colours leaving node in given orientation
int db_node_count_number_of_colours_out(dBNode * node, Orientation orientation)
{
	int i;
	int number_of_colours_out = 0;

	for (i = 0; i < NUMBER_OF_COLOURS; i++) {
		if (db_node_edges_count_by_colour(node, orientation, i) > 0)
			number_of_colours_out++;
	}

	return number_of_colours_out;
}

int db_node_count_number_of_colours_out_any_orientation(dBNode * node)
{
	int i;
	int number_of_colours_out = 0;

	for (i = 0; i < NUMBER_OF_COLOURS; i++) {
		if ((db_node_edges_count_by_colour(node, reverse, i) > 0) ||
		    (db_node_edges_count_by_colour(node, forward, i) > 0))
			number_of_colours_out++;
	}

	return number_of_colours_out;
}

boolean db_node_is_blunt_end_all_colours(dBNode * node, Orientation orientation)
{

	Edges edges = db_node_get_edges_all_colours(node);

	if (orientation == reverse) {
		edges >>= 4;
	}

	edges &= 15;		// AND with 00001111 so that we only look at the 4 least significant bits

	return edges == 0;
}

boolean db_node_is_blunt_end(dBNode * node, Orientation orientation)
{

	Edges edges = db_node_get_edges(node);

	if (orientation == reverse) {
		edges >>= 4;
	}

	edges &= 15;		// AND with 00001111 so that we only look at the 4 least significant bits

	return edges == 0;
}

// Gets edges by orientation for ALL colours (they are OR'd together)
Edges
db_node_get_edges_for_orientation_all_colours(dBNode * node,
					      Orientation orientation)
{
	Edges edges = 0;
	short i;

	for (i = 0; i < NUMBER_OF_COLOURS; i++) {
		edges |=
		    db_node_get_edges_for_orientation_by_colour(node,
								orientation, i);
	}
	return edges;
}

// Gets edges by orientation for a single speciified colour
Edges
db_node_get_edges_for_orientation_by_colour(dBNode * node,
					    Orientation orientation,
					    short colour)
{

	Edges edges = db_node_get_edges_by_colour(node, colour);
	if (orientation == reverse) {
		edges >>= 4;
	}
	edges &= 15;		// AND with 00001111 so that we only look at the 4 least significant bits
	return edges;
}

// Defaults to colour 0
Edges db_node_get_edges_for_orientation(dBNode * node, Orientation orientation)
{
	return db_node_get_edges_for_orientation_by_colour(node, orientation,
							   0);
}

boolean
db_node_is_visited_on_all_the_paths(dBNode * node, Orientation orientation)
{
	fprintf(stderr, "[db_node_is_visited_on_all_the_paths] Not yet implemented, \
			if implemented it should iterate to the next nodes and check if any of \
			the next nodes is free\n");
	exit(-1);
	//return db_node_get_edges_for_orientation(node,
	//					 orientation ==
	//					 (node->flags | EDGES_FLAGS));
}

/** Note: reads single instance of coverage and edges - only store one
*       colour per binary file. 
* 
*/

boolean db_node_read_binary(FILE * fp, short kmer_size, dBNode * node)
{
	BinaryKmer kmer;
	Edges edges;
	uint32_t coverage;
	int read;

	//int number_of_bitfields = ((kmer_size * 2) / (sizeof(bitfield_of_64bits)*8))+1;
	//read = fread(&kmer,number_of_bitfields*sizeof(bitfield_of_64bits), 1,fp);
    
	read = fread(&kmer, sizeof(bitfield_of_64bits), NUMBER_OF_BITFIELDS_IN_BINARY_KMER, fp);

	if (read > 0) {
		element_initialise(node, &kmer, kmer_size);	//Moved the initialization at the begining to be able to keep the read pair in a single blok of code

		read = fread(&coverage, sizeof(uint32_t), 1, fp);
		if (read == 0) {
			puts("error with input file\n");
			exit(1);
		}

		read = fread(&edges, sizeof(Edges), 1, fp);
		if (read == 0) {
			puts("error with input file\n");
			exit(1);
		}
	} else {
		return false;
	}

	node->edges[0] = edges;
	node->coverage[0] = coverage;

#ifdef SOLID
#warning solid0
	Edges first_base;
	read = fread(&first_base, sizeof(Edges), 1, fp);
	if (read == 0) {
		puts("error with input file\n");
		exit(1);
	}

//	db_node_set_all_starting_bases(starting_base, node);
	node->first_base = first_base;
#endif 

	return true;
}

/** Note: writes single instance of coverage and edges - only store one
*         colour per binary file. 
*   Note2: Only writes one pair. each CTX shall be generated from each pair.
*  WARN: We are not checking that the code is compiled with solid
*/
void db_node_print_binary(FILE * fp, dBNode * node, int kmer_size)
{
	BinaryKmer kmer;
	binary_kmer_assignment_operator(kmer, *element_get_kmer(node));
	Edges edges = db_node_get_edges(node);
	//int number_of_bitfields = ((kmer_size * 2) / (sizeof(bitfield_of_64bits)*8))+1;
	uint32_t total_coverage = element_get_coverage_all_colours(node);
	uint32_t coverage = node->coverage[0];
	
	// Can only write this node if all colour 0
	if (coverage != total_coverage) {
		fclose(fp);
		fprintf(stderr,
			"[db_node_print_binary] All colours should be 0. File incomplete.");
		exit(-1);
	}

	//fwrite(&kmer,  number_of_bitfields*sizeof(bitfield_of_64bits), 1, fp);
    
	fwrite(&kmer, sizeof(bitfield_of_64bits), NUMBER_OF_BITFIELDS_IN_BINARY_KMER, fp);
	fwrite(&coverage, sizeof(uint32_t), 1, fp);
	fwrite(&edges, sizeof(Edges), 1, fp);
#ifdef SOLID
#warning solid1
	//Edges first_base = db_node_get_all_starging_bases(node);
	Edges first_base = node->first_base;
	fwrite(&first_base, sizeof(Edges), 1, fp );
#endif
}

/** Note: writes single instance of coverage and edges - only store one
 *         colour per binary file. 
 *   Note2: Only writes one pair. each CTX shall be generated from each pair.
 */
void
db_node_print_binary_by_colour(FILE * fp, dBNode * node, short colour,
			       int kmer_size)
{
	BinaryKmer kmer;
	binary_kmer_assignment_operator(kmer, *element_get_kmer(node));
	Edges edges = db_node_get_edges_by_colour(node, colour);
	//int number_of_bitfields = ((kmer_size * 2) / (sizeof(bitfield_of_64bits)*8))+1;
	uint32_t coverage = node->coverage[colour];

	//fwrite(&kmer,  number_of_bitfields*sizeof(bitfield_of_64bits), 1, fp);

	fwrite(&kmer, sizeof(bitfield_of_64bits), NUMBER_OF_BITFIELDS_IN_BINARY_KMER, fp);
	fwrite(&coverage, sizeof(uint32_t), 1, fp);
	fwrite(&edges, sizeof(Edges), 1, fp);

#ifdef ENABLE_READ_PAIR
	//ReadPairSignature signature = db_node_get_signature(0, 0, node);
	//fwrite(&signature, sizeof(ReadPairSignature), 1, fp);
#endif

#ifdef SOLID
#warning solid2
	//Edges first_base = db_node_get_all_starging_bases(node);
	Edges first_base = node->first_base;
	fwrite(&first_base, sizeof(Edges), 1, fp );
#endif
}

// Get edges for a given colour
Edges db_node_get_edges_by_colour(dBNode * node, short colour)
{
	return node->edges[colour];
}

// db_node_get_edges is 'colourblid' - gets colour 0
Edges db_node_get_edges(dBNode * node)
{
	return db_node_get_edges_by_colour(node, 0);
}

// db_node_get_edges is 'colourblid' - gets colour 0
Edges db_node_get_edges_all_colours(dBNode * node)
{
	int i;
	Edges e = 0;

	for (i = 0; i < NUMBER_OF_COLOURS; i++)
		e |= db_node_get_edges_by_colour(node, i);

	return e;
}

void db_node_action_set_current_path(dBNode * node, Orientation o)
{
	if (o == forward)
		db_node_action_set_flag(node, CURRENT_PATH_FORWARD);
	else
		db_node_action_set_flag(node, CURRENT_PATH_REVERSE);
}

void db_node_action_unset_current_path(dBNode * node, Orientation o)
{
	if (o == forward)
		db_node_action_unset_flag(node, CURRENT_PATH_FORWARD);
	else
		db_node_action_unset_flag(node, CURRENT_PATH_REVERSE);
}

boolean db_node_action_is_in_current_path(dBNode * node, Orientation o)
{
	if (o == forward)
		return db_node_check_for_flag(node, CURRENT_PATH_FORWARD);
	else
		return db_node_check_for_flag(node, CURRENT_PATH_REVERSE);

}

void db_node_action_do_nothing(dBNode * node)
{

}

boolean db_node_condition_always_true(dBNode * node)
{
	return true;
}

void db_node_action_clear_flags(dBNode * node)
{
	//node->flags = ALL_OFF;
	Flags f = (ASSIGNED | PRUNED) & node->flags;
	flags_action_clear_flags(&(node->flags));
	flags_action_set_flag(f, &(node->flags));
}

void db_node_action_set_flag(dBNode * node, Flags f)
{
	// if (DEBUG)
//printf("SET: %x | %x = ", node->flags, f);
	//node->flags = node->flags | f;

	flags_action_set_flag(f, &(node->flags));
	//if (DEBUG)
	//printf("%x\n", node->flags);
}

void db_node_action_unset_flag_current_path(dBNode * node)
{
	db_node_action_unset_flag(node,
				  CURRENT_PATH_FORWARD | CURRENT_PATH_REVERSE);
}

void db_node_action_unset_flag(dBNode * node, Flags f)
{
	//if(DEBUG) printf("UNSET: %x & ~%x = %x\n", node->flags, f, (node->flags & ~f));
	//node->flags = node->flags & ~f;
	flags_action_unset_flag(f, &(node->flags));
}

Flags db_node_get_flags(dBNode * node, Flags f)
{
	//if (DEBUG)
	//printf("GETING: %x & %x = %x\n", node->flags, f, (node->flags & f));
	return node->flags & f;
}
//The next functions asume we have a big set of flags. We are trying to reduce the footprint, so if
//they provide to be absolutly necesary, we should bring them back
#ifndef SHORT_FLAGS
void db_node_action_set_visited(dBNode * node, Orientation o, Nucleotide n,
			   Nucleotide n_r)
{
	//TODO: check if it really uses the path from reverse properly...

	Flags f = Undefined == n ? 0 : (1 << n);
	Flags f_r = Undefined == n_r ? 0 : (1 << n_r);
	if (o == forward) {
		db_node_action_set_flag(node, VISITED_FORWARD);
		f_r <<= 4;
	} else {
		db_node_action_set_flag(node, VISITED_REVERSE);
		f <<= 4;
	}

	if (DEBUG)
		printf("Setting flag: %x, for nucleotide %d\n", f | f_r, n);
	db_node_action_set_flag(node, f | f_r);

}

void db_node_action_unset_visited(dBNode * node, Orientation o)
{
	o == forward ? db_node_action_unset_flag(node, VISITED_FORWARD)
	    : db_node_action_unset_flag(node, VISITED_REVERSE);
}

boolean db_node_check_visited(dBNode * node, Orientation o)
{
	return (o == forward) ? db_node_check_for_flag(node, VISITED_FORWARD)
	    : db_node_check_for_flag(node, VISITED_REVERSE);

}
#endif
boolean db_node_check_for_flag(dBNode * node, Flags flag)
{
	//if(DEBUG) printf("CHECK: %x & %x = %x\n", node->flags, flag, (node->flags & flag));

	return flags_check_for_flag(flag, &(node->flags));

}

boolean db_node_check_for_any_flag(dBNode * node, Flags flag)
{
	//if (DEBUG)
	//      printf("CHECK FOR ANY: %x & %x = %x\n", node->flags, flag, (node->flags
	//                      & flag));
	return flags_check_for_any_flag(flag, &(node->flags));
}
#ifndef SHORT_FLAGS
Flags
db_node_set_print_orientation(Orientation current_orientation,
			      dBNode * current_node)
{
	Flags f = ALL_OFF;

	   f = db_node_get_flags(current_node, PRINT_FORWARD | PRINT_REVERSE);

	if (f == ALL_OFF) {

		if (current_orientation == forward) {

			db_node_action_set_flag(current_node, PRINT_FORWARD);

			f = PRINT_FORWARD;
		} else {
		
			db_node_action_set_flag(current_node, PRINT_REVERSE);

			f = PRINT_REVERSE;
		}
	}
	return f;
}
#endif
void db_node_action_set_flag_pruned(dBNode * node)
{
	db_node_action_set_flag(node, PRUNED);
}

void db_node_action_set_flag_visited(dBNode * node)
{
    if (!(node->flags & VISITED)) {
        visited_count++;
        db_node_action_set_flag(node, VISITED);
    }
}

long long int get_visited_count(void) {
    return visited_count;
}

void clear_visited_count(void) {
    visited_count = 0;
}

boolean db_node_check_flag_not_pruned(dBNode * node)
{
	return !db_node_check_for_flag(node, PRUNED);
}

boolean db_node_check_flag_visited(dBNode * node)
{
	return db_node_check_for_flag(node, VISITED);
}

boolean db_node_check_for_flag_ALL_OFF(dBNode * node)
{
	return node->flags == ALL_OFF;
}

boolean element_check_for_flag_ALL_OFF(Element * node)
{
	return node->flags == ALL_OFF;
}

boolean element_check_for_flag(Element * node, Flags flag)
{
	return db_node_check_for_flag((dBNode *) node, flag);
}

#ifdef ENABLE_READ_PAIR_OLD

boolean db_node_action_add_read_pair(long long count, short pair, short colour, dBNode * node)
{
	//ReadPairSignature old = node->signature[colour][pair];
	ReadPairSignature old = db_node_get_signature( pair,  colour, node);
	READ_PAIR_DATATYPE one = 1;
	
	node->signature[colour][pair] |= one << ((count - 1) % READ_PAIR_LENGTH); 
	
//	printf("count: %lld, pair: %i, old %d, new: %d\n", count, pair, old,  node->signature[colour][pair]);
	
	return node->signature[colour][pair] != old;	//returns true if the signature changed the node or not.
}

int db_node_count_signature_matches(ReadPairSignature signature, short pair, short colour, dBNode * node)
{
	ReadPairSignature matches = signature & db_node_get_signature(pair, colour, node);
	ReadPairSignature tmp;
	int i;
	int count = 0;
	for (i = 0; i < READ_PAIR_LENGTH; i++) {
		tmp = 1 << i;
		if (tmp & matches) {
			count++;
		}
	}
	return count;
}

ReadPairSignature db_node_get_signature(short pair, short colour, dBNode * node)
{
	return node->signature[colour][pair];
}

int db_node_set_signature(ReadPairSignature signature, short pair, short colour, dBNode * node)
{
	int matches = db_node_count_signature_matches(signature, pair, colour, node);
	node->signature[colour][pair] |= signature;
	return matches;
}
#elif defined ENABLE_MARK_PAIR

uint32_t db_node_get_supernode(dBNode * node){
    return node->supernode;
}

void db_node_set_supernode_id(uint32_t id, dBNode * node){
    node->supernode = id;
}



#elif defined ENABLE_READ_PAIR

ReadPairSignature db_node_get_signature(short pair, ReadPairBitfield bitfield, dBNode * node)
{
//    short signature_index = (pair * 4) + bitfield;

    short signature_index = pair;
    if(bitfield == Bitfield_PerfectPath){
        signature_index = NUMBER_OF_SIGNATURES - 1;
    }
    assert(pair >= 0);
    assert(bitfield >= 0);
    assert(bitfield <= 4);
    
    return node->signatures[signature_index];
}

int db_node_count_signature_matches(ReadPairSignature signature, short pair, ReadPairBitfield bitfield, dBNode * node)
{
    assert(pair >= 0);
    assert(bitfield >= 0);
    assert(bitfield <= 4);
    
	ReadPairSignature matches = signature & db_node_get_signature(pair, bitfield, node);
	ReadPairSignature tmp;
	int i;
	int count = 0;

    

	for (i = 0; i < READ_PAIR_LENGTH; i++) {
		tmp = 1 << i;
		if (tmp & matches) {
			count++;
		}
	}
    
	return count;
}

int db_node_set_signature(ReadPairSignature signature, short pair, ReadPairBitfield bitfield, dBNode * node)
{
    assert(pair >= 0);
    assert(bitfield >= 0);
    assert(bitfield <= 4);
    
//    short signature_index = (pair * 4) + bitfield;
     short signature_index = pair;
    if (bitfield == Bitfield_PerfectPath) {
        signature_index = NUMBER_OF_SIGNATURES - 1;
    }
 //   int matches = db_node_count_signature_matches(signature, pair, bitfield, node);
//TODO: validate there are no collisions. 
    
    
    node->signatures[signature_index] |= signature;
    	    
    return 0;
}

boolean db_node_action_add_read_pair(long long count, short pair, ReadPairBitfield bitfield, dBNode * node)
{
    assert(pair >= 0);
    assert(bitfield >= 0);
    assert(bitfield <= 4);
    
    //    short signature_index = (pair * 4) + bitfield;

    short signature_index = pair;
    if (pair == Bitfield_PerfectPath) {
        signature_index = NUMBER_OF_SIGNATURES - 1;
    }
	ReadPairSignature old = db_node_get_signature(pair, bitfield, node);
	READ_PAIR_DATATYPE one = 1;

   

//    node->signatures[signature_index] |= one << ((count - 1) % READ_PAIR_LENGTH);
    node->signatures[signature_index] |= one << ((count - 1) % READ_PAIR_LENGTH);
	return node->signatures[signature_index] != old;	//returns true if the signature changed the node or not.
}

#endif

//Code to deal with the SOLiD specific stuff. 
#ifdef SOLID

Edges db_node_get_all_starging_bases(dBNode * node){
	return node->first_base;
}

void db_node_set_all_starting_bases(Edges e, dBNode * node){
	node->first_base = e;
}

void db_node_add_starting_base(NucleotideBaseSpace base, Orientation o, dBNode * node){
    //set edge
    if(base != Undef){ //Do nothing if the base is not defineds. 
        char e = 1 << base;	// A (0) -> 0001, C (1) -> 0010, G (2) -> 0100, T (3) -> 1000
        if(o == reverse){
            e <<= 4 ;
        }
        node->first_base |= e;
    }
}

NucleotideBaseSpace db_node_get_starting_base(Orientation orientation, dBNode * node){
    Edges edges = node->first_base;
    NucleotideBaseSpace base = Undef;
    boolean found, valid;
    Edges mask = 1;
    Edges current;
    NucleotideBaseSpace i;
    
    if (orientation == reverse) {
		edges >>= 4;
	}
	edges &= 15;
    
    
    valid = true;
    found = false;
    for(i = 0; i < 4; i++){
        current = edges & mask;
        if(current && found){
            valid = false;
        }
        if(current){
            base = i;
            found = true;
        }
        mask <<= 1;
    }
    
    if(!valid){
        base = Undef;
    }
    return base;
    
}

NucleotideBaseSpace db_node_get_starting_base_search_reverse(Orientation orientation, dBNode * node,short kmer_size){
    NucleotideBaseSpace base = Undef;
    
    base = db_node_get_starting_base(orientation, node);
    
    if(base == Undef){
        base = db_node_get_starting_base(opposite_orientation(orientation), node);
    }
    
    base = binary_kmer_get_last_nucleotide_in_base_space(&node->kmer, orientation, base, kmer_size );
    
    return base;
    
}


#endif
