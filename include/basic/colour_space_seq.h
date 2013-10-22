/**
*
* This set of functions work on conventional Sequences containing just the base space
* TODO: This relies on the binary kmer which has the SOLiD table, we may want to 
* make this independent of the binary kmer if we plan to write mor SOLiD based
* tools.  
*/
#ifndef COLOUR_SEQ_H_
#define COLOUR_SEQ_H_

#ifdef SOLID

/**
* Converts a sequence in colour space to base space. 
*/
Sequence * cs_sequence_to_base_sequence(Sequence * seq_cs, Sequence * seq_base);

/**
* Converts a sequence in base space to colour space. 
*/
Sequence * base_sequence_to_cs_sequence(Sequence * seq_base, Sequence * seq_cs);
#endif
#endif