/*
 * test_alighment.c
 *
 */

#include <CUnit.h>
#include <stdlib.h>
#include <stdio.h>
#include <nucleotide.h>
#include <seq.h>
#include <alignment.h>
#include <test_alignment.h>

void test_perfect_alignment(){
  Alignment * al;
  al = alignment_new(100, 100);
  CU_ASSERT(al != NULL);
  Sequence * a = calloc(1, sizeof(Sequence));
  Sequence * b = calloc(1, sizeof(Sequence));
  
  alloc_sequence(a, 100,100, 0);
  alloc_sequence(b, 100,100, 0);
  sequence_set_name("A", a);
  sequence_set_name("B", b);
  sequence_append("ACTGA", "IIIJI",a);
  sequence_append("ACTGA", "IJIII",b);
  
  CU_ASSERT_STRING_EQUAL("ACTGA", a->seq);
  CU_ASSERT_STRING_EQUAL("ACTGA", b->seq);
  
  alignment_set_sequences(a,b,al);
  
  CU_ASSERT(a == al->a);
  CU_ASSERT(b == al->b);
  
  alignment_align(al, NULL);
  
  CU_ASSERT_STRING_EQUAL("ACTGA", al->a_aligned->seq);
  CU_ASSERT_STRING_EQUAL("ACTGA", al->b_aligned->seq);
  
}

void test_alignment_with_gap_1(){
  Alignment * al;
  al = alignment_new(100, 100);
  CU_ASSERT(al != NULL);
  Sequence * a = calloc(1, sizeof(Sequence));
  Sequence * b = calloc(1, sizeof(Sequence));
  
  alloc_sequence(a, 100,100, 0);
  alloc_sequence(b, 100,100, 0);
  sequence_set_name("A", a);
  sequence_set_name("B", b);
  sequence_append("ACTGA", "IIIII",a);
  sequence_append("CTGA", "IJII",b);
  
  CU_ASSERT_STRING_EQUAL("ACTGA", a->seq);
  CU_ASSERT_STRING_EQUAL("CTGA", b->seq);
  
  alignment_set_sequences(a,b,al);
  
  CU_ASSERT(a == al->a);
  CU_ASSERT(b == al->b);
  
  alignment_align(al, NULL);
  
  CU_ASSERT_STRING_EQUAL("CTGA", al->a_aligned->seq);
  CU_ASSERT_STRING_EQUAL("CTGA", al->b_aligned->seq);
  
}

void test_alignment_with_gap_2(){
  Alignment * al;
  al = alignment_new(100, 100);
  CU_ASSERT(al != NULL);
  Sequence * a = calloc(1, sizeof(Sequence));
  Sequence * b = calloc(1, sizeof(Sequence));
  
  alloc_sequence(a, 100,100, 0);
  alloc_sequence(b, 100,100, 0);
  sequence_set_name("A", a);
  sequence_set_name("B", b);
  sequence_append("ACTG", "IIII",a);
  sequence_append("ACTGA", "IIIII",b);
  
  CU_ASSERT_STRING_EQUAL("ACTG", a->seq);
  CU_ASSERT_STRING_EQUAL("ACTGA", b->seq);
  
  alignment_set_sequences(a,b,al);
  
  CU_ASSERT(a == al->a);
  CU_ASSERT(b == al->b);
  
  alignment_align(al, NULL);
  
  CU_ASSERT_STRING_EQUAL("ACTG", al->a_aligned->seq);
  CU_ASSERT_STRING_EQUAL("ACTG", al->b_aligned->seq);
  
}

void test_alignment_with_gap_3(){
  Alignment * al;
  al = alignment_new(100, 100);
  CU_ASSERT(al != NULL);
  Sequence * a = calloc(1, sizeof(Sequence));
  Sequence * b = calloc(1, sizeof(Sequence));
  
  alloc_sequence(a, 100,100, 0);
  alloc_sequence(b, 100,100, 0);
  sequence_set_name("A", a);
  sequence_set_name("B", b);
  sequence_append("TTACACACTAGGGCCC", "IIIIIIIIIIIII",a);
  sequence_append("TTAGCACACAGGGCCC", "IIIIIIIIIIIII",b);
  
  CU_ASSERT_STRING_EQUAL("TTACACACTAGGGCCC", a->seq);
  CU_ASSERT_STRING_EQUAL("TTAGCACACAGGGCCC", b->seq);
  
  alignment_set_sequences(a,b,al);
  
  CU_ASSERT(a == al->a);
  CU_ASSERT(b == al->b);
  
  alignment_align(al, &default_nucleotide_score);
  
  CU_ASSERT_STRING_EQUAL("TTA-CACACTAGGGCCC", al->a_aligned->seq);
  CU_ASSERT_STRING_EQUAL("TTAGCACAC-AGGGCCC", al->b_aligned->seq);
  
}


void test_initialized_alignment(){
	Alignment * al;
	al = alignment_new(100, 100);
	CU_ASSERT(al != NULL);
	al->a_aligned->seq[0] = 'N';
	al->a_aligned->seq[1] = 'O';
	al->a_aligned->seq[2] = 'I';
	al->a_aligned->seq[3] = 'S';
	al->a_aligned->seq[4] = 'E';
	al->a_aligned->seq[5] = '\0';
	alignment_clean(al);
	CU_ASSERT_STRING_EQUAL("", al->a_aligned->seq);
	alignment_destroy(&al);
	CU_ASSERT(al == NULL);
	
}

void test_alignment_virus(){
  Alignment * al;
  al = alignment_new(5000, 5000);
  CU_ASSERT(al != NULL);
  Sequence * a = calloc(1, sizeof(Sequence));
  Sequence * b = calloc(1, sizeof(Sequence));
  
  alloc_sequence(a, 5000,5000, 0);
  alloc_sequence(b, 5000,5000, 0);
  sequence_set_name("Flu217_1_1273333h12.p1k", a);
  sequence_set_name("B", b);
  sequence_append("TAAGCTTTTGGCCCCCCCCCCCTATCTTGTGCATATTTTCCCC", "",a);
  sequence_append("TCTTGTGCATATTTTCCCCTATGAATGCTCTGAACTATTCAGT", "",b);
  
  
  alignment_set_sequences(a,b,al);
  
  CU_ASSERT(a == al->a);
  CU_ASSERT(b == al->b);
  
  alignment_align(al, NULL);
  
  CU_ASSERT_STRING_EQUAL("TCTTGTGCATATTTTCCCC", al->a_aligned->seq);
  CU_ASSERT_STRING_EQUAL("TCTTGTGCATATTTTCCCC", al->b_aligned->seq);
  printf("A[%d, %d], B[%d, %d]", al->start_a, al->end_a,al->start_b, al->end_b);
  
}

void test_alignment_homopolymer_ins(){
  Alignment * al;
  al = alignment_new(5000, 5000);
  CU_ASSERT(al != NULL);
  Sequence * a = calloc(1, sizeof(Sequence));
  Sequence * b = calloc(1, sizeof(Sequence));
  
  alloc_sequence(a, 5000,5000, 0);
  alloc_sequence(b, 5000,5000, 0);
  sequence_set_name("Reference", a);
  sequence_set_name("B", b);
  sequence_append("TAAGCTTTTGGCTTCTT", "",a);
			
  sequence_append("TAAGCTTTTTGGCTTCTT", "",b);
  
  
  alignment_set_sequences(a,b,al);
  
  CU_ASSERT(a == al->a);
  CU_ASSERT(b == al->b);
  
  alignment_align(al, &half_gap_nucleotide_score);
  alignment_fix_homopolymers(al);
  CU_ASSERT_STRING_EQUAL("TAAGCTTT-TGGCTTCTT", al->a_aligned->seq);
  CU_ASSERT_STRING_EQUAL("TAAGCTTTTTGGCTTCTT", al->b_aligned->seq);
  printf("A[%d, %d], B[%d, %d]", al->start_a, al->end_a,al->start_b, al->end_b);
  CU_ASSERT_STRING_EQUAL("TAAGCTTTTGGCTTCTT", al->consensus->seq);
  
}

void test_alignment_non_homopolymer_ins(){
  Alignment * al;
  al = alignment_new(5000, 5000);
  CU_ASSERT(al != NULL);
  Sequence * a = calloc(1, sizeof(Sequence));
  Sequence * b = calloc(1, sizeof(Sequence));
  
  alloc_sequence(a, 5000,5000, 0);
  alloc_sequence(b, 5000,5000, 0);
  sequence_set_name("Reference", a);
  sequence_set_name("B", b);
  sequence_append("TAAGCTTTTGGCTTCTT", "",a);
  sequence_append("TAAGCATTTTGGCTTCTT", "",b);
  
  
  alignment_set_sequences(a,b,al);

  CU_ASSERT(a == al->a);
  CU_ASSERT(b == al->b);
  
  alignment_align(al, &half_gap_nucleotide_score);
   alignment_fix_homopolymers(al);
    
  CU_ASSERT_STRING_EQUAL("TAAG-CTTTTGGCTTCTT", al->a_aligned->seq);
  CU_ASSERT_STRING_EQUAL("TAAGCATTTTGGCTTCTT", al->b_aligned->seq);
  printf("A[%d, %d], B[%d, %d]", al->start_a, al->end_a,al->start_b, al->end_b);
  CU_ASSERT_STRING_EQUAL("TAAGCTTTTGGCTTCTT", al->consensus->seq);
  
}

void test_alignment_homopolymer_del(){
  Alignment * al;
  al = alignment_new(5000, 5000);
  CU_ASSERT(al != NULL);
  Sequence * a = calloc(1, sizeof(Sequence));
  Sequence * b = calloc(1, sizeof(Sequence));
  
  alloc_sequence(a, 5000,5000, 0);
  alloc_sequence(b, 5000,5000, 0);
  sequence_set_name("Reference", a);
  sequence_set_name("B", b);
  sequence_append("TAAGCTTTTGGCTTCTT", "",a);
  sequence_append("TAAGCTTTGGCTTCTT", "",b);
  
  
  alignment_set_sequences(a,b,al);
  
  CU_ASSERT(a == al->a);
  CU_ASSERT(b == al->b);
  
  alignment_align(al, &half_gap_nucleotide_score);
  alignment_fix_homopolymers(al);
						  
  CU_ASSERT_STRING_EQUAL("TAAGCTTTTGGCTTCTT", al->a_aligned->seq);
  CU_ASSERT_STRING_EQUAL("TAAGCTT-TGGCTTCTT", al->b_aligned->seq);
  printf("A[%d, %d], B[%d, %d]", al->start_a, al->end_a,al->start_b, al->end_b);
  CU_ASSERT_STRING_EQUAL("TAAGCTTTTGGCTTCTT", al->consensus->seq);
  
}

void test_alignment_non_homopolymer_del(){
  Alignment * al;
  al = alignment_new(5000, 5000);
  CU_ASSERT(al != NULL);
  Sequence * a = calloc(1, sizeof(Sequence));
  Sequence * b = calloc(1, sizeof(Sequence));
  
  alloc_sequence(a, 5000,5000, 0);
  alloc_sequence(b, 5000,5000, 0);
  sequence_set_name("Reference", a);
  sequence_set_name("B", b);
  sequence_append("TAAGCTTTTGGCTTCTT", "",a);
  sequence_append("TAAGTTTTGGCTTCTT", "",b);
  
  
  alignment_set_sequences(a,b,al);
  
  CU_ASSERT(a == al->a);
  CU_ASSERT(b == al->b);
  
  alignment_align(al, &half_gap_nucleotide_score);
  alignment_fix_homopolymers(al);
  
  
  CU_ASSERT_STRING_EQUAL("TAAGCTTTTGGCTTCTT", al->a_aligned->seq);
  CU_ASSERT_STRING_EQUAL("TAA-GTTTTGGCTTCTT", al->b_aligned->seq);
  printf("A[%d, %d], B[%d, %d]", al->start_a, al->end_a,al->start_b, al->end_b);
  CU_ASSERT_STRING_EQUAL("TAAGNTTTTGGCTTCTT", al->consensus->seq);
  
}

void test_alignment_homopolymer_mut(){
  Alignment * al;
  al = alignment_new(5000, 5000);
  CU_ASSERT(al != NULL);
  Sequence * a = calloc(1, sizeof(Sequence));
  Sequence * b = calloc(1, sizeof(Sequence));
  
  alloc_sequence(a, 5000,5000, 0);
  alloc_sequence(b, 5000,5000, 0);
  sequence_set_name("Reference", a);
  sequence_set_name("B", b);
  sequence_append("TAAGCTTTTGGCTTCTT", "",a);
  sequence_append("TAAGCTTTCGGCTTCTT", "",b);
  
  
  alignment_set_sequences(a,b,al);
  
  
  CU_ASSERT(a == al->a);
  CU_ASSERT(b == al->b);
  
  alignment_align(al, &half_gap_nucleotide_score);
  alignment_fix_homopolymers(al);
  
  CU_ASSERT_STRING_EQUAL("TAAGCTTTTGGCTTCTT", al->a_aligned->seq);
  CU_ASSERT_STRING_EQUAL("TAAGCTTTCGGCTTCTT", al->b_aligned->seq);
  printf("A[%d, %d], B[%d, %d]", al->start_a, al->end_a,al->start_b, al->end_b);
  CU_ASSERT_STRING_EQUAL("TAAGCTTTCGGCTTCTT", al->consensus->seq);
  
}
