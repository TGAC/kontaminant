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
 element.h defines the interface for the de Bruijn graph node. The implementation is complemented by
 a hash table that stores every node indexed by kmers (BinaryKmers).

 The element routines, ie the one required by hash_table/priority queue, are prefixed with element_
 The de Bruijn based routines are prefixed with db_node
 */

#include <global.h>
#include <nucleotide.h>
#include <seq.h>
#include <binary_kmer.h>
#include <flags.h>

#ifndef ELEMENT_H_
#define ELEMENT_H_

#ifndef NUMBER_OF_COLOURS
#ifdef ENABLE_BUBBLEPARSE
#define NUMBER_OF_COLOURS 2
#else
#define NUMBER_OF_COLOURS 1
#endif
#endif


#ifdef ENABLE_READ_PAIR

#endif
#ifdef ENABLE_READ_PAIR

#ifndef READ_PAIR_DATATYPE
#define READ_PAIR_DATATYPE unsigned long long
//#define READ_PAIR_DATATYPE unsigned int
#endif

#ifndef READ_PAIR_LENGTH
#define READ_PAIR_LENGTH   (sizeof(READ_PAIR_DATATYPE)*8)

#endif

#ifndef READ_PAIR_SEGMENTS
#define READ_PAIR_SEGMENTS 1
#endif

// Before marking by perfect paths
//#define NUMBER_OF_SIGNATURES (4 * READ_PAIR_SEGMENTS) 
// To mark the perfect path and different segments
#define NUMBER_OF_SIGNATURES READ_PAIR_SEGMENTS + 1 

typedef READ_PAIR_DATATYPE ReadPairSignature;

typedef enum {
	Bitfield_First=0, Bitfield_Second=1, Bitfield_First_Inverse=2, Bitfield_Second_Inverse=3, Bitfield_PerfectPath=4
} ReadPairBitfield;

#endif

#ifdef ENABLE_READ_TRACK
#ifndef READ_TRACK_DATATYPE
#define READ_TRACK_DATATYPE unsigned long long
#endif

#ifndef READ_PAIR_LENGTH
#define READ_PAIR_LENGTH   (sizeof(READ_TRACK_DATATYPE)*8)
#endif

#define NUMBER_OF_SIGNATURES 2

#endif

//type definitions
typedef char Edges;



// We provide a smaller version of element for use in count_kmers
typedef struct {
	BinaryKmer kmer;
	uint32_t coverage[NUMBER_OF_COLOURS];
	// less significant nibble forward
	Edges edges[NUMBER_OF_COLOURS];
#ifdef SOLID
    Edges first_base;
#endif
	Flags flags;
#ifdef ENABLE_READ_PAIR_OLD
	ReadPairSignature signature[NUMBER_OF_COLOURS][READ_PAIR_SEGMENTS];
#elif defined ENABLE_READ_PAIR
    ReadPairSignature signatures[NUMBER_OF_SIGNATURES];
#endif
#ifdef INCLUDE_QUALITY_SCORES
	QualityStringArray quality_string_arrays[NUMBER_OF_COLOURS];
#endif
    
#ifdef ENABLE_MARK_PAIR
    uint32_t supernode;
#endif


} Element;


typedef struct {
	BinaryKmer kmer;
	Edges edges; // less significant nibble forward
} Element2;

typedef Element dBNode;

typedef BinaryKmer* Key;

typedef Element GraphNode;

void element_assign(Element* e1, Element* e2);

//reverse orientation
Orientation opposite_orientation(Orientation);

boolean element_is_key(Key, Element, short kmer_size);

Key element_get_key(BinaryKmer*, short kmer_size, Key preallocated_key);

//boolean element_smaller(Element, Element);

void element_initialise(Element *, BinaryKmer* kmer, short kmer_size);

#ifdef INCLUDE_QUALITY_SCORES
void element_add_quality_string(Element *e, short c, char *q);
void element_preallocate_quality_strings(Element *e, int c, int n);
#endif

BinaryKmer* element_get_kmer(Element *);

//int element_get_coverage(Element *);
int element_get_coverage_all_colours(Element *);
int element_get_coverage_by_colour(Element *, short);
int element_update_coverage(Element *, short, int);
boolean db_node_check_for_flag_ALL_OFF(dBNode * node);

Orientation db_node_get_orientation(BinaryKmer*, dBNode *, short kmer_size);

//add an edge between nodes -- NB: it adds both edges: forward and reverse
boolean db_node_add_edge(dBNode *, dBNode *, Orientation, Orientation,
		short kmer_size, short);

//returns true if the node side defined by the orientation is a conflict 
//or doesn't have any outgoing edge
boolean db_node_is_supernode_end(dBNode *, Orientation);

//returns yes if the label defined by the nucleotide coresponds to an 
//outgoing edge in the side defined by the orientation.   
boolean db_node_edge_exist(dBNode *, Nucleotide, Orientation);
boolean db_node_edge_exist_any_colour(dBNode *, Nucleotide, Orientation);

//returns the label of the first outgoing edge -- leaving from the side 
//defined by orientation. 
boolean db_node_has_precisely_one_edge(dBNode *, Orientation, Nucleotide *);

boolean db_node_has_precisely_one_edge_all_colours(dBNode *, Orientation, Nucleotide *);

//returns the label of the "two edges"
//defined by orientation. 
boolean db_node_has_precisely_two_edges(dBNode *, Orientation,
		Nucleotide *, Nucleotide *);

boolean db_node_has_unvisited_edge(dBNode * node, Orientation orientation,
		Nucleotide * nucleotide);

Edges db_node_get_edges(dBNode * node);

Edges db_node_get_edges_all_colours(dBNode * node);

Edges db_node_get_edges_by_colour(dBNode * node, short colour);

char db_node_get_edges_coverage(dBNode * node);

//forgets about the edges
void db_node_reset_edges(dBNode *, short);

void db_node_reset_edges_all_colours(dBNode * node);

void db_node_reset_edge(dBNode *, Orientation, short, Nucleotide);

void db_node_reset_edge_all_colours(dBNode *, Orientation, Nucleotide);

//check that the edges are 0's
//boolean db_node_edges_reset(dBNode *);

//set every edge in 'edges' 
void db_node_set_edges(dBNode * node, short colour, Edges edges);

boolean db_node_edge_has_single_coverage(dBNode * element, Nucleotide base,
		Orientation orientation);

int db_node_count_number_of_colours_out(dBNode *node, Orientation orientation);

int db_node_count_number_of_colours_out_any_orientation(dBNode *node);

//check if node doesn't have any edges in a given orientation
boolean db_node_is_blunt_end(dBNode * node, Orientation orientation);

boolean db_node_is_blunt_end_all_colours(dBNode * node, Orientation orientation);

void db_node_print_binary(FILE * fp, dBNode * node,int kmer_size);

void db_node_print_binary_by_colour(FILE * fp, dBNode * node, short colour, int kmer_size);

boolean db_node_read_binary(FILE * fp, short kmer_size, dBNode * node);

//actions and conditions 

void db_node_action_do_nothing(dBNode * node);

void db_node_action_clear_flags(dBNode * node);

void db_node_action_set_flag(dBNode * node, Flags f);

void db_node_action_unset_flag(dBNode * node, Flags f);

void db_node_action_set_flag_none(dBNode * node);

void db_node_action_set_flag_pruned(dBNode * node);

void db_node_action_set_flag_visited(dBNode * node);

long long int get_visited_count(void);

void clear_visited_count(void);

Flags db_node_get_flags(dBNode * node, Flags f);
#ifndef SHORT_FLAGS
void db_node_action_set_visited(dBNode * node, Orientation o, Nucleotide n,
		Nucleotide n_r);

void db_node_action_unset_visited(dBNode * node, Orientation o);

boolean db_node_check_visited(dBNode * node, Orientation o);
#endif
void db_node_action_set_current_path(dBNode * node, Orientation o);

void db_node_action_unset_current_path(dBNode * node, Orientation o);

boolean db_node_action_is_in_current_path(dBNode * node, Orientation o);

void db_node_action_unset_flag_current_path(dBNode * no);

boolean db_node_check_flag_visited(dBNode * node) ;

boolean db_node_check_for_flag(dBNode * node, Flags flag);

boolean db_node_check_for_any_flag(dBNode * node, Flags flag);

boolean db_node_check_nothing(dBNode * node);

int db_node_edges_count(dBNode * node, Orientation orientation);

char binary_nucleotide_to_edge(Nucleotide base);

int db_node_edges_count_by_colour(dBNode * node, Orientation orientation, short colour);

int db_node_edges_count_all_colours(dBNode * node, Orientation orientation);

boolean db_node_condition_always_true(dBNode* node);

boolean db_node_is_visited_on_all_the_paths(dBNode * node,
		Orientation orientation);

Edges db_node_get_edges_for_orientation_by_colour(dBNode * node, Orientation orientation, short colour);

Edges db_node_get_edges_for_orientation_all_colours(dBNode * node, Orientation orientation);

Edges db_node_get_edges_for_orientation(dBNode * node, Orientation orientation);
/**
 * Method to set the orientation in which the node will be printed. It is
 * important to notice that if the node already has an orientation, the method
 * doesn't have any effect and the returned value is the stored orientation.
 */
#ifndef SHORT_FLAGS
Flags db_node_set_print_orientation(Orientation current_orientation,
		dBNode * current_node);
#endif
boolean element_check_for_flag_ALL_OFF(Element * node);

boolean db_node_check_flag_not_pruned(dBNode * node);

void db_node_action_set_flag_none(dBNode * node);

boolean elemet_is_assigned(Element * node);

boolean element_check_for_flag(Element * node, Flags flag);

/**
 * 
 *Functions related to the read pair. They are in an ifdef block
 * to keep compatibility with programs that do not require the read pair
 * functionality  
 * 
 */
#ifdef ENABLE_READ_PAIR_OLD

/**
 * Function that sets the corresponding bit in the signature for the sequence count. 
 * We are assuming that the read pair files come in the same order and have the same number
 * of entries.The return value helps in counting the "collisions".
 * 
 * Returns true when setting the flag has an effect and false when the flag was already set. 
 */ 
boolean db_node_action_add_read_pair(long long count, short pair, short colour, dBNode *  node);

/**
 * Counts how many of the bits are present in the node from a given signature
 */
int db_node_count_signature_matches(ReadPairSignature signature, short pair, short colour, dBNode * node); 


/**
 * Gets the signature of a given pair in a given colour. 
 */ 
ReadPairSignature db_node_get_signature(short pair, short colour, dBNode * node);

/**
 * Function that sets the signature for a given node. It does and OR operation to set the field
 * It returns the number of bits that where already set up in the corresponing colour and pair. 
 * The return value is helpful to know how many of the merging signatures collide. 
 * 
 */ 
int db_node_set_signature(ReadPairSignature signature,  short pair, short colour, dBNode * node);

#elif defined ENABLE_READ_PAIR

ReadPairSignature db_node_get_signature(short pair, ReadPairBitfield bitfield, dBNode * node);
int db_node_count_signature_matches(ReadPairSignature signature, short pair, ReadPairBitfield bitfield, dBNode * node);
int db_node_set_signature(ReadPairSignature signature, short pair, ReadPairBitfield bitfield, dBNode * node);
boolean db_node_action_add_read_pair(long long count, short pair, ReadPairBitfield bitfield, dBNode * node);
#elif defined ENABLE_MARK_PAIR

uint32_t db_node_get_supernode(dBNode * node);

void db_node_set_supernode_id(uint32_t id, dBNode * node);

#endif

#ifdef SOLID
/**
 * Function that sets certain node as a read-start node, from which we can start the translation. 
 */ 
void db_node_add_starting_base(NucleotideBaseSpace base, Orientation o, dBNode * node);

/**
 * Returns the starting base, provided that there is one and only one base in the given direction
 */
 
NucleotideBaseSpace db_node_get_starting_base(Orientation o, dBNode * node);

  
NucleotideBaseSpace db_node_get_starting_base_search_reverse(Orientation orientation, dBNode * node,short kmer_size);

//Returns the inner structure, for writing the CTX. 
Edges db_node_get_all_starging_bases(dBNode * node);

//Overrides the inner representation, for reading the CTX
void db_node_set_all_starting_bases(Edges e, dBNode * node);
#endif

#endif /* ELEMENT_H_ */
