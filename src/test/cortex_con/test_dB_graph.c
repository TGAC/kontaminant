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
#include <CUnit.h>
#include <stdlib.h>
#include <Basic.h>

#include <global.h>
#include <flags.h>
#include <inttypes.h>
#include <nucleotide.h>
#include <binary_kmer.h>
#include <seq.h>
#include <element.h>

#include <db_graph.h>
#include <cleaning.h>
#include <file_reader.h>
#include <perfect_path.h>
#include <y_walk.h>
#include <branches.h>
#include <binary_tree.h>

#ifdef ENABLE_READ_PAIR

#include <read_pair.h>
#endif

#include <test_dB_graph.h>


void test_hash_table_find()
{
	
	//first set up the hash/graph
	int kmer_size = 3;
	int number_of_bits = 4;
	int bucket_size    = 4;
	long long bad_reads = 0; 
	
	int max_retries=10;
	
	dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
	
	//Load the following fasta:
	//    >read1
	//    AAAAAAAA
	//    >read2
	//    GGCT
	//    >read3
	//    TAGG
	
#ifdef ENABLE_READ_PAIR_OLD
	int seq_length = load_fasta_from_filename_into_graph("../data/test/graph/test_dB_graph.fasta",0,0,&bad_reads, 20, db_graph);
#else
	int seq_length = load_fasta_from_filename_into_graph("../data/test/graph/test_dB_graph.fasta",0,&bad_reads, 20, db_graph);
#endif
	//length of total sequence
	CU_ASSERT(seq_length == 16);
	
	//number of kmers
	CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph), 5);
	
	//bad reads
	CU_ASSERT_EQUAL(bad_reads, 0);
	
	//all the kmers and their reverse complements from the reads
	BinaryKmer tmp_kmer1;
	BinaryKmer tmp_kmer2;
	
	dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AAA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("TTT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* test_element3 = hash_table_find(element_get_key(seq_to_binary_kmer("GGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* test_element4 = hash_table_find(element_get_key(seq_to_binary_kmer("GCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* test_element5 = hash_table_find(element_get_key(seq_to_binary_kmer("GCT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* test_element6 = hash_table_find(element_get_key(seq_to_binary_kmer("AGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* test_element7 = hash_table_find(element_get_key(seq_to_binary_kmer("TAG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* test_element8 = hash_table_find(element_get_key(seq_to_binary_kmer("CTA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* test_element9 = hash_table_find(element_get_key(seq_to_binary_kmer("AGG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* test_element10 = hash_table_find(element_get_key(seq_to_binary_kmer("CCT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	
	//some kmers that should not be in the graph
	
	dBNode* test_element11 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* test_element12 = hash_table_find(element_get_key(seq_to_binary_kmer("CCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* test_element13 = hash_table_find(element_get_key(seq_to_binary_kmer("TAT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* test_element14 = hash_table_find(element_get_key(seq_to_binary_kmer("ATA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* test_element15 = hash_table_find(element_get_key(seq_to_binary_kmer("TAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* test_element16 = hash_table_find(element_get_key(seq_to_binary_kmer("ATG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* test_element17 = hash_table_find(element_get_key(seq_to_binary_kmer("TTG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* test_element18 = hash_table_find(element_get_key(seq_to_binary_kmer("AAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* test_element19 = hash_table_find(element_get_key(seq_to_binary_kmer("TGA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* test_element20 = hash_table_find(element_get_key(seq_to_binary_kmer("TCA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	
	//kmers in the graph
	CU_ASSERT(test_element1 != NULL);
	CU_ASSERT(test_element2 != NULL);
	CU_ASSERT(test_element1 == test_element2);
	
	CU_ASSERT(test_element3 != NULL);
	CU_ASSERT(test_element4 != NULL);
	CU_ASSERT(test_element3 == test_element4);
	
	CU_ASSERT(test_element5 != NULL);
	CU_ASSERT(test_element6 != NULL);  
	CU_ASSERT(test_element5 == test_element6);
	
	CU_ASSERT(test_element7 != NULL);
	CU_ASSERT(test_element8 != NULL);
	CU_ASSERT(test_element7 == test_element8);
	
	CU_ASSERT(test_element9 != NULL);
	CU_ASSERT(test_element10 != NULL);
	CU_ASSERT(test_element9 == test_element10);  
	
	//kmers not in the graph
	CU_ASSERT(test_element11 == NULL);
	CU_ASSERT(test_element12 == NULL);
	CU_ASSERT(test_element13 == NULL);
	CU_ASSERT(test_element14 == NULL);
	CU_ASSERT(test_element15 == NULL);
	CU_ASSERT(test_element16 == NULL);
	CU_ASSERT(test_element17 == NULL);
	CU_ASSERT(test_element18 == NULL);
	CU_ASSERT(test_element19 == NULL);
	CU_ASSERT(test_element20 == NULL);
	
	
	hash_table_free(&db_graph);
	CU_ASSERT(db_graph == NULL);
	
}


void test_tip_clipping()
{
	
	//first set up the hash/graph
	int kmer_size = 3;
	int number_of_bits = 4;
	int bucket_size    = 4;
	long long bad_reads = 0;
	path_array_initialise_buffers(kmer_size);
	
	
	BinaryKmer tmp_kmer1;
	BinaryKmer tmp_kmer2;
	
	dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
	
	
	//Load the following fasta:
	
	//>main_trunk
	//GCGTCCCAT
	//>tip
	//CGTTT
#ifdef ENABLE_READ_PAIR_OLD
	int seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_tip.fasta",0,0,&bad_reads, 20, db_graph);
#else
	int seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_tip.fasta",0,&bad_reads, 20, db_graph);
#endif
	
	CU_ASSERT_EQUAL(seq_length,14);
	
	dBNode* node1 = hash_table_find(element_get_key(seq_to_binary_kmer("TTT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	
	CU_ASSERT(node1 != NULL);
	
	db_graph_write_graphviz_file("test_tip_clipping_before.viz", db_graph);
	int tip_length = db_graph_db_node_clip_tip(node1, 10,&db_node_action_set_flag_pruned,db_graph);
	db_graph_write_graphviz_file("test_tip_clipping_after.viz", db_graph);
	CU_ASSERT_EQUAL(tip_length,2);
	
	dBNode* node2 = hash_table_find(element_get_key(seq_to_binary_kmer("GTT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* node3 = hash_table_find(element_get_key(seq_to_binary_kmer("CGT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* node4 = hash_table_find(element_get_key(seq_to_binary_kmer("CCA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	
	
	//check the change in the status
	CU_ASSERT(!db_node_check_flag_not_pruned(node1));
	CU_ASSERT(!db_node_check_flag_not_pruned(node2));
	
	//check status didn't change
	CU_ASSERT(db_node_check_for_flag(node3,ASSIGNED));
	CU_ASSERT(db_node_check_for_flag(node4,ASSIGNED));
	
	//check tip clip works well with db_graph_supernode (ie that the tip was really removed)
	Path * tmp_path = path_new(100, db_graph->kmer_size);
	
	int length_supernode =  db_graph_supernode(node3, &db_node_action_set_flag_visited,
											   tmp_path, db_graph);
	/*int length_supernode = db_graph_supernode(node3,100,&db_node_action_set_flag_visited,
	 nodes_path,orientations_path,labels_path,
	 tmp_seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
	 db_graph);*/
	
	CU_ASSERT_EQUAL(length_supernode,6);
	CU_ASSERT_STRING_EQUAL("GGACGC",tmp_path->seq);
	//fprintf(stderr, "Expected: GGACGC, obtained %s. len: %d\n", tmp_path->seq, length_supernode);
	
	
	//check ends
	dBNode* node5 = hash_table_find(element_get_key(seq_to_binary_kmer("GCG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* node6 = hash_table_find(element_get_key(seq_to_binary_kmer("CAT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	
	//check the status are correct
	
	CU_ASSERT(db_node_check_for_flag(node1,PRUNED));
	CU_ASSERT(db_node_check_for_flag(node2,PRUNED));
	
	CU_ASSERT(db_node_check_for_flag(node3,VISITED));  
	CU_ASSERT(db_node_check_for_flag(node4,VISITED));
	CU_ASSERT(db_node_check_for_flag(node5,VISITED));
	CU_ASSERT(db_node_check_for_flag(node6,VISITED));
	
	
	hash_table_free(&db_graph);
	CU_ASSERT(db_graph == NULL);
	
	
}

void test_node_prunning_low_coverage()
{
	
	//first set up the hash/graph
	int kmer_size = 3;
	int number_of_bits = 4;
	int bucket_size    = 4;
	long long bad_reads = 0;
	BinaryKmer tmp_kmer1;
	BinaryKmer tmp_kmer2;
	
	dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
	
	
	//Load the following fasta:
	
	//>main_trunk
	//GCGTCCCAT
	//>tip
	//CGTTT
	
#ifdef ENABLE_READ_PAIR_OLD
	int seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_tip.fasta",0,0,&bad_reads, 20, db_graph);
#else
	int seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_tip.fasta",0,&bad_reads, 20, db_graph);
#endif
	CU_ASSERT_EQUAL(seq_length,14);
	
	dBNode* node1 = hash_table_find(element_get_key(seq_to_binary_kmer("CCC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* node2 = hash_table_find(element_get_key(seq_to_binary_kmer("TCC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* node3 = hash_table_find(element_get_key(seq_to_binary_kmer("CCA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode* node4 = hash_table_find(element_get_key(seq_to_binary_kmer("ACG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	
	db_graph_write_graphviz_file("test_node_prunning_low_coverage_before.viz", db_graph);
	boolean node_pruned = db_graph_db_node_prune_low_coverage(node1,1,&db_node_action_set_flag_pruned,db_graph);
	db_graph_write_graphviz_file("test_node_prunning_low_coverage_after.viz", db_graph);
	
	
	CU_ASSERT(node_pruned);
	CU_ASSERT(db_node_check_for_flag(node1,PRUNED));
	
	CU_ASSERT(!db_node_edge_exist(node3,Guanine,reverse));
	CU_ASSERT(!db_node_edge_exist(node3,Thymine,forward));
	
	CU_ASSERT(!db_node_edge_exist(node2,Cytosine,reverse));
	CU_ASSERT(!db_node_edge_exist(node2,Cytosine,forward));
	
	CU_ASSERT(!db_node_check_for_flag(node4,PRUNED));
	//check all arrows were removed from node1
	Nucleotide n;
    for (n = Adenine ;n<Undefined; n++) {
        CU_ASSERT(!db_node_edge_exist(node1,n,reverse));
		CU_ASSERT(!db_node_edge_exist(node1,n,forward));
    }
		
	
	
	//nucleotide_iterator(&nucleotide_action);
	hash_table_free(&db_graph);
	
	CU_ASSERT(db_graph == NULL);
}


void test_get_perfect_path() //test db_graph_get_perfect_path
{
	
	//first set up the hash/graph
	int kmer_size = 3;
	int number_of_bits=5;
	int bucket_size   = 10;
	int seq_length;
	long long bad_reads = 0;
	/*dBNode * path_nodes[100];
	 Orientation path_orientations[100];
	 Nucleotide path_labels[100];
	 char tmp_seq[100];
	 
	 double  avg_coverage;
	 int min_coverage, max_coverage;
	 */
	dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
	
	
	//1. Sequence of tests as follows
	//         Each test loads a single specifically designed fasta file into a dB_graph.
	//         The test then picks an element in the graph, and calls get_perfect_path
	//         and checks that it gets the right sequence.
	
	
	// ****
	//1.1 Fasta file that generate a graph with two hairpins, and a single edge (in each rorientation) joining them.
	//  Sequence is :  ACGTAC
	// ****
	
#ifdef ENABLE_READ_PAIR_OLD
	seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_two_self_loops.fasta",0,0,&bad_reads, 20, db_graph);
#else
	seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_two_self_loops.fasta",0,  &bad_reads, 20,  db_graph);
#endif
	CU_ASSERT_EQUAL(seq_length,6);
	CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),2);
	CU_ASSERT_EQUAL(bad_reads,0);
	
	
	// now start at GTA and get all the sequence from there to the end of the supernode, and see
	// if that is right.
	
	BinaryKmer tmp_kmer1;
	BinaryKmer tmp_kmer2;
	
	dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GTA", kmer_size,&tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
	CU_ASSERT(test_element1!=NULL);
	
	boolean is_cycle=false;
	
	Path * tmp_path = path_new(100, db_graph->kmer_size);
	
	int test1_length =  db_graph_supernode(test_element1, &db_node_action_set_flag_visited,
										   tmp_path, db_graph);
	
	/*int test1_length = db_graph_get_perfect_path(test_element1,forward,100,
	 &db_node_action_do_nothing,
	 path_nodes,path_orientations,path_labels,
	 tmp_seq, &avg_coverage, &min_coverage, &max_coverage,
	 &is_cycle,db_graph);*/
    db_graph_write_graphviz_file("generates_graph_with_two_self_loops.viz", db_graph);
    is_cycle = path_is_cycle(tmp_path);
    CU_ASSERT(is_cycle);
    CU_ASSERT_EQUAL(test1_length,4);
	
    CU_ASSERT_STRING_EQUAL(tmp_path->seq,"CGTA");
	
	hash_table_free(&db_graph);
	CU_ASSERT(db_graph == NULL);
	
	
	/*   // **** */
	/*   //1.2 Fasta file that generate a graph with one long supernode, with a conflict at the end */
	/*   //   caused by two outward/exiting edges */
	/*   // **** */
	
	//first set up the hash/graph
	kmer_size = 3;
	number_of_bits= 4;
	bucket_size   = 10;
	
	bad_reads = 0;
	dBGraph * db_graph2 = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
	//db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
#ifdef ENABLE_READ_PAIR_OLD
	seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_one_long_supernode_with_conflict_at_end.fasta",0,0,&bad_reads, 20, db_graph2);
#else
	seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_one_long_supernode_with_conflict_at_end.fasta",0, &bad_reads,20,db_graph2);
#endif
	db_graph_write_graphviz_file("test2.viz", db_graph2);
	CU_ASSERT_EQUAL(seq_length,13);
	CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph2),5);
	CU_ASSERT_EQUAL(bad_reads,0);
	
	test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("ACA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph2);
	CU_ASSERT(test_element1!=NULL);
	char seq[20];
	CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(element_get_kmer(test_element1),3,&seq[0]),"ACA");
	
	//ACA < TGT so forward gives TT
	is_cycle=false;
	path_reset(tmp_path);
	
	//go forward
	/* int test2_length = db_graph_get_perfect_path(test_element1,forward,100,
	 &db_node_action_do_nothing,
	 path_nodes,path_orientations,path_labels,
	 tmp_seq, &avg_coverage, &min_coverage, &max_coverage,
	 &is_cycle,db_graph);*/
	
    int test2_length = db_graph_get_perfect_path(test_element1,undefined,&db_node_action_do_nothing,db_graph2, tmp_path);
	
	
	//fprintf(stderr, "Expected: TT obtained: %s, len: %d", tmp_path->seq, test2_length);
	CU_ASSERT(!path_is_cycle(tmp_path));
	CU_ASSERT_EQUAL(test2_length,2);
	CU_ASSERT_STRING_EQUAL(tmp_path->seq,"TT");
	
	
	//ACA < TGT so backward gives ""
	
	//go backward
	test2_length = db_graph_get_perfect_path(test_element1,reverse,&db_node_action_do_nothing, db_graph2, tmp_path);
	/*test2_length = db_graph_get_perfect_path(test_element1,reverse,100,
	 &db_node_action_do_nothing,
	 path_nodes,path_orientations,path_labels,
	 tmp_seq, &avg_coverage, &min_coverage, &max_coverage,
	 &is_cycle,db_graph);
	 */
	
	CU_ASSERT(!path_is_cycle(tmp_path));
	CU_ASSERT_EQUAL(test2_length,0);
	
	
	CU_ASSERT_STRING_EQUAL(tmp_path->seq,"");
	printf("We got as seq: %s\n",  tmp_path->seq);
	hash_table_free(&db_graph2);
	CU_ASSERT(db_graph2 == NULL);
	
	
	// ****
	//1.3 Fasta file that generate a graph with one long supernode, with a conflict at the end
	//   caused by two INWARD edges in the opposite direction
	// ****
	
	//first set up the hash/graph
	kmer_size = 3;
	number_of_bits = 4;
	bucket_size    = 5;
	bad_reads = 0;
	
	db_graph = hash_table_new(number_of_bits, bucket_size,10,kmer_size);
	
#ifdef ENABLE_READ_PAIR_OLD
	seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_one_long_supernode_with_inward_conflict_at_end.fasta",0,0,&bad_reads, 20, db_graph);
#else
	seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_one_long_supernode_with_inward_conflict_at_end.fasta",0,  &bad_reads, 20,  db_graph);
#endif
	CU_ASSERT_EQUAL(seq_length,13);
	CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),5);
	CU_ASSERT_EQUAL(bad_reads,0);
	
	
	test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("ACA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	CU_ASSERT(test_element1!=NULL);
	
	//ACA < TGT so forward gives TT
	is_cycle=false;
	
	//go forward
	int test3_length = db_graph_get_perfect_path(test_element1,forward,
												 &db_node_action_do_nothing, db_graph, tmp_path);
	
	CU_ASSERT(!path_is_cycle(tmp_path));
	
	printf("_-_-_-_-_-   WE GOT: %s\n", tmp_path->seq);
	CU_ASSERT_EQUAL(test3_length,2);
	CU_ASSERT_STRING_EQUAL(tmp_path->seq,"TT");
	
	
	//ACA < TGT so backwar gives ""
	is_cycle=false;
	
	//go reverse
	test3_length = db_graph_get_perfect_path(test_element1,reverse,
											 &db_node_action_do_nothing, db_graph, tmp_path);
	
	CU_ASSERT(!path_is_cycle(tmp_path));
	
	CU_ASSERT_EQUAL(test3_length,0);
	CU_ASSERT_STRING_EQUAL(tmp_path->seq,"");
	
	hash_table_free(&db_graph);
	CU_ASSERT(db_graph == NULL);
	
	
	// ****
	//1.4 Fasta file that generate a graph with an infinite loop at a single kmer
	//
	// ****
	
	
	//first set up the hash/graph
	kmer_size = 3;
	number_of_bits=8;
	bucket_size   =4;
	
	bad_reads = 0;
	
	db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
#ifdef ENABLE_READ_PAIR_OLD
	seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_infinite_loop.fasta",0,0,&bad_reads, 20, db_graph);
#else
	seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_infinite_loop.fasta",0, &bad_reads,30,db_graph);
#endif
	CU_ASSERT_EQUAL(seq_length,25);
	CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),1);
	CU_ASSERT_EQUAL(bad_reads,0);
	
	test_element1 = hash_table_find(seq_to_binary_kmer("AAA", kmer_size, &tmp_kmer1) ,db_graph);
	CU_ASSERT(test_element1!=NULL);
	
	//forward
	is_cycle=false;
	
	int test4_length = db_graph_get_perfect_path(test_element1,forward,
												 &db_node_action_do_nothing,
												 db_graph, tmp_path);
	CU_ASSERT(path_is_cycle(tmp_path));//this time it should find a cycle
	CU_ASSERT_STRING_EQUAL(tmp_path->seq,"A");
	CU_ASSERT_EQUAL(test4_length,1);
	
	//backward
	is_cycle=false;
	
	test4_length = db_graph_get_perfect_path(test_element1,reverse,
											 &db_node_action_do_nothing,
											 db_graph, tmp_path);
	
	CU_ASSERT(path_is_cycle(tmp_path));//this time it should find a cycle
	CU_ASSERT_EQUAL(test4_length,1);
	CU_ASSERT_STRING_EQUAL(tmp_path->seq,"T");
	
	hash_table_free(&db_graph);
	CU_ASSERT(db_graph == NULL);
	
	
	
	// ****
	// 1.5 check parameters (path nodes,labels,etc) for get_perfect_path
	//
	// ****
	
	kmer_size = 3;
	number_of_bits = 4;
	bucket_size = 4;
	db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
	
	boolean found1, found2, found3;
	dBNode * node1;
	dBNode * node2;
	dBNode * node3;
	dBNode * node4;
	
	
	node1 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CGT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2), &found1, db_graph);
	node2 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("GTT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found2, db_graph);
	
	db_node_add_edge(node1, node2, reverse,reverse, db_graph->kmer_size, 0);
	
	node3 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("TTA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found3, db_graph);
	db_node_add_edge(node2, node3, reverse,reverse, db_graph->kmer_size, 0);
	
	node4 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("TAG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found3, db_graph);
	db_node_add_edge(node3, node4, reverse,reverse, db_graph->kmer_size, 0);
	
	
	//add coverage
	element_update_coverage(node1,0,1);
	element_update_coverage(node2,0,2);
	element_update_coverage(node3,0,1);
	element_update_coverage(node4,0,2);
	
	
	
	/* dBNode * nodes[10];
	 Orientation orientations[10];
	 Nucleotide bases[10];*/
	int test5_length = 0;
	
	
	test5_length = db_graph_get_perfect_path(node1,reverse, 
											 &db_node_action_set_flag_visited,
											 db_graph, tmp_path);
	
	CU_ASSERT_EQUAL(test5_length,3);
	
	//check nodes
	CU_ASSERT_EQUAL(node1,tmp_path->nodes[0]);
	CU_ASSERT_EQUAL(node2,tmp_path->nodes[1]);
	CU_ASSERT_EQUAL(node3,tmp_path->nodes[2]);
	CU_ASSERT_EQUAL(node4,tmp_path->nodes[3]);
	
	//check labels
	CU_ASSERT_EQUAL(tmp_path->labels[0],Thymine);
	CU_ASSERT_EQUAL(tmp_path->labels[1],Adenine);
	CU_ASSERT_EQUAL(tmp_path->labels[2],Guanine);
	
	//check orientations
	CU_ASSERT_EQUAL(tmp_path->orientations[0],reverse);
	CU_ASSERT_EQUAL(tmp_path->orientations[1],reverse);
	CU_ASSERT_EQUAL(tmp_path->orientations[2],reverse);
	
	//check statuses
	CU_ASSERT(db_node_check_for_flag(tmp_path->nodes[0], ASSIGNED));
	CU_ASSERT(db_node_check_for_flag(tmp_path->nodes[1], VISITED));
	CU_ASSERT(db_node_check_for_flag(tmp_path->nodes[2], VISITED));
	CU_ASSERT(db_node_check_for_flag(tmp_path->nodes[3], ASSIGNED));
	
	double avg;
	int min, max;
	
	path_get_statistics(&avg, &min, &max, tmp_path);
	
	
	db_graph_write_graphviz_file("Test1_5.viz", db_graph);
	
	
	//check coverage
	CU_ASSERT_EQUAL(avg,1.5);
	CU_ASSERT_EQUAL(min,1);
	CU_ASSERT_EQUAL(max,2);
	
	
	hash_table_free(&db_graph);
	CU_ASSERT(db_graph == NULL);
}


void test_writing_reading_binary(){
	
	dBNode node1, node2;
	int kmer_size = 3;
	boolean test_read;
	
	char seq[kmer_size];
	
	BinaryKmer tmp_kmer;
	
	element_initialise(&node1,seq_to_binary_kmer("AAA",kmer_size, &tmp_kmer),kmer_size);
	node1.edges[0] = 'a'; // ie 64 = (0110 0100)2
	node1.coverage[0] = 10;
	
	FILE* fp1 = fopen("../data/test/graph/dump_element.ctx", "w");
	db_node_print_binary(fp1,&node1, kmer_size);
	fclose(fp1);
	
	FILE* fp2 = fopen("../data/test/graph/dump_element.ctx", "r");
	
	test_read = db_node_read_binary(fp2,kmer_size,&node2);
	
	CU_ASSERT_EQUAL(test_read, true);
	
	CU_ASSERT(binary_kmer_comparison_operator(node1.kmer, node2.kmer));
	CU_ASSERT_EQUAL(node1.edges[0], node2.edges[0]);
	
	CU_ASSERT_STRING_EQUAL("AAA", binary_kmer_to_seq(&(node1.kmer),kmer_size,seq));
	CU_ASSERT_STRING_EQUAL("AAA", binary_kmer_to_seq(&(node2.kmer),kmer_size,seq));
	
	CU_ASSERT_EQUAL(node1.edges[0], node2.edges[0]);
	CU_ASSERT_EQUAL(node1.edges[0], 'a');
	CU_ASSERT_EQUAL(node2.edges[0], 'a');
	
	CU_ASSERT_EQUAL(node2.coverage[0], 10);
	
	test_read = db_node_read_binary(fp2,kmer_size,&node2);
	CU_ASSERT_EQUAL(test_read, false);
	
	fclose(fp2);
}

void test_detect_and_smoothe_bubble(){
	int kmer_size = 3;
	int number_of_buckets = 4;
	int bucket_size = 3;
	dBGraph * db_graph = hash_table_new(number_of_buckets,bucket_size,10,kmer_size);
	boolean found;
	BinaryKmer tmp_kmer1;
	BinaryKmer tmp_kmer2;
	
	
	dBNode * node1, * node2, * node3, * node4, * node5, * node6, * node7, * node8, * node9;
	
	//start point
	node1 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2), &found, db_graph);
	
	//branch1
	node2 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	
	
	//modify coverage -- used below to clip one branch
	element_update_coverage(node2,0, 5);
	CU_ASSERT_EQUAL(element_get_coverage_all_colours(node2),5);
	
	
	node3 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CAC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	node4 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("ACT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	
	//end point (branch merge here)
	node5 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CTT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	
	//branch 2
	node6 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	
	element_update_coverage(node6, 0,1);
	CU_ASSERT_EQUAL(element_get_coverage_all_colours(node6),1);
	
	node7 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CTC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	node8 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("TCT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	
	//add 3p extension
	node9 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("TTA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	
	//branch1
	db_node_add_edge(node1, node2, forward,forward, db_graph->kmer_size, 0);
	db_node_add_edge(node2, node3, forward,forward, db_graph->kmer_size, 0);
	db_node_add_edge(node3, node4, forward,forward, db_graph->kmer_size, 0);
	db_node_add_edge(node4, node5, forward,reverse, db_graph->kmer_size, 0);
	
	
	//branch2
	db_node_add_edge(node1, node6, forward,reverse, db_graph->kmer_size,0);
	db_node_add_edge(node6, node7, reverse,forward, db_graph->kmer_size,0);
	db_node_add_edge(node7, node8, forward,reverse, db_graph->kmer_size,0);
	db_node_add_edge(node8, node5, reverse,reverse, db_graph->kmer_size,0);
	
	//add 3p extension
	db_node_add_edge(node5, node9, reverse,reverse, db_graph->kmer_size, 0);
	
	CU_ASSERT(db_node_edge_exist(node1,Thymine,forward));
	CU_ASSERT(db_node_edge_exist(node1,Adenine,forward));  
	db_graph_write_graphviz_file("test_detect_and_smoothe_bubble_before.viz", db_graph);
	
	cleaning_remove_bubbles(db_graph->kmer_size + 10,  db_graph);
	db_graph_write_graphviz_file("test_detect_and_smoothe_bubble_after.viz", db_graph);
	
	//the first node shouldn't be marked
	
	CU_ASSERT(db_node_check_for_flag(node1,VISITED));
	CU_ASSERT(db_node_check_for_flag(node2,VISITED));
	CU_ASSERT(db_node_check_for_flag(node3,VISITED));
	CU_ASSERT(db_node_check_for_flag(node4,VISITED));
	
	//This nodes should been pruned
	CU_ASSERT(!db_node_check_for_flag(node1,PRUNED));
	CU_ASSERT(!db_node_check_for_flag(node2,PRUNED));
	CU_ASSERT(!db_node_check_for_flag(node3,PRUNED));
	CU_ASSERT(!db_node_check_for_flag(node4,PRUNED));
	CU_ASSERT(!db_node_check_for_flag(node9,PRUNED));
	
	CU_ASSERT(db_node_check_for_flag(node6,PRUNED));
	CU_ASSERT(db_node_check_for_flag(node7,PRUNED));
	CU_ASSERT(db_node_check_for_flag(node8,PRUNED));
	//the last node shouldn't be marked 
	CU_ASSERT(db_node_check_for_flag(node9,ASSIGNED));
	CU_ASSERT(!db_node_check_for_flag(node9,VISITED));
}


void test_detect_and_remove_low_cov_path(){
	int kmer_size = 3;
	int number_of_buckets = 4;
	int bucket_size = 3;
	dBGraph * db_graph = hash_table_new(number_of_buckets,bucket_size,10,kmer_size);
	boolean found;
	BinaryKmer tmp_kmer1;
	BinaryKmer tmp_kmer2;
	
	
	dBNode * node1, * node2, * node3, * node4, * node5, * node6, * node7, * node8, * node9, * node11, * node12, * node13, * node21, * node22, * node23;
	
	//start point
	node1 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2), &found, db_graph);
	element_update_coverage(node1, 0,5); 
	CU_ASSERT_EQUAL(element_get_coverage_all_colours(node1),5);
	
	//branch1
	node2 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	
	
	//modify coverage -- used below to clip one branch
	element_update_coverage(node2,0, 5);
	CU_ASSERT_EQUAL(element_get_coverage_all_colours(node2),5);
	
	
	node3 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CAC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	element_update_coverage(node3,0, 5);
	CU_ASSERT_EQUAL(element_get_coverage_all_colours(node3),5);
	node4 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("ACT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	element_update_coverage(node4,0, 5);
	CU_ASSERT_EQUAL(element_get_coverage_all_colours(node4),5);
	
	//end point (branch merge here)
	node5 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CTT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	element_update_coverage(node5, 0,5); 
	CU_ASSERT_EQUAL(element_get_coverage_all_colours(node5),5);
	
	//branch 2
	node6 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	element_update_coverage(node6, 0,1); 
	CU_ASSERT_EQUAL(element_get_coverage_all_colours(node6),1);
	
	node7 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CTC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	element_update_coverage(node7, 0,1); 
	CU_ASSERT_EQUAL(element_get_coverage_all_colours(node7),1);
	node8 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("TCT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	element_update_coverage(node8, 0,1); 
	CU_ASSERT_EQUAL(element_get_coverage_all_colours(node8),1);
	
	//add 3p extension
	node9 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("TTA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	element_update_coverage(node9, 0,6); 
	
	CU_ASSERT_EQUAL(element_get_coverage_all_colours(node9),6); 
	
	//unconnected path
	node21 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("ACA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	element_update_coverage(node21, 0,1); 
	CU_ASSERT_EQUAL(element_get_coverage_all_colours(node21),1);
	node22 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CAA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	element_update_coverage(node22, 0,3); 
	CU_ASSERT_EQUAL(element_get_coverage_all_colours(node22),3);
	node23 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("AAC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	element_update_coverage(node23, 0,1); 
	CU_ASSERT_EQUAL(element_get_coverage_all_colours(node23),1);
	
	//unconnected high_cov
	node11 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("AAA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	element_update_coverage(node11, 0,1); 
	CU_ASSERT_EQUAL(element_get_coverage_all_colours(node11),1);
	node12 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("AAT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	element_update_coverage(node12, 0,1); 
	CU_ASSERT_EQUAL(element_get_coverage_all_colours(node12),1);
	node13 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("ATA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
	element_update_coverage(node13, 0,1); 
	CU_ASSERT_EQUAL(element_get_coverage_all_colours(node13),1);
	
	//branch1
	db_node_add_edge(node1, node2, forward,forward, db_graph->kmer_size, 0);
	db_node_add_edge(node2, node3, forward,forward, db_graph->kmer_size, 0);
	db_node_add_edge(node3, node4, forward,forward, db_graph->kmer_size, 0);
	db_node_add_edge(node4, node5, forward,reverse, db_graph->kmer_size, 0);
	
	
	//branch2
	db_node_add_edge(node1, node6, forward,reverse, db_graph->kmer_size,0);
	db_node_add_edge(node6, node7, reverse,forward, db_graph->kmer_size,0);
	db_node_add_edge(node7, node8, forward,reverse, db_graph->kmer_size,0);
	db_node_add_edge(node8, node5, reverse,reverse, db_graph->kmer_size,0);
	
	//add 3p extension
	db_node_add_edge(node5, node9, reverse,reverse, db_graph->kmer_size, 0);
	
	
	//unconnected path
	db_node_add_edge(node11, node12, forward,forward, db_graph->kmer_size, 0);
	db_node_add_edge(node12, node13, forward,forward, db_graph->kmer_size, 0);
	
	//unconnected path high_cov
	db_node_add_edge(node21, node22, forward,forward, db_graph->kmer_size, 0);
	db_node_add_edge(node22, node23, forward,forward, db_graph->kmer_size, 0);
	
	CU_ASSERT(db_node_edge_exist(node1,Thymine,forward));
	CU_ASSERT(db_node_edge_exist(node1,Adenine,forward));  
	
	
	db_graph_write_graphviz_file("test_detect_and_remove_low_cov_path_before.viz", db_graph);
	
	
	cleaning_prune_low_coverage_path(1, 10, db_graph);
	
	
	db_graph_write_graphviz_file("test_detect_and_remove_low_cov_path_after.viz", db_graph);
	
	
	//the first node shouldn't be marked
	CU_ASSERT(!db_node_check_for_flag(node1,PRUNED));
	CU_ASSERT(!db_node_check_for_flag(node2,PRUNED));
	CU_ASSERT(!db_node_check_for_flag(node3,PRUNED));
	CU_ASSERT(!db_node_check_for_flag(node4,PRUNED));
	CU_ASSERT(db_node_check_for_flag(node6,PRUNED));
	CU_ASSERT(db_node_check_for_flag(node7,PRUNED));
	CU_ASSERT(db_node_check_for_flag(node8,PRUNED));
	//the last node shouldn't be marked 
	CU_ASSERT(!db_node_check_for_flag(node9,PRUNED));
	
	//The low coverage floating path should be removed
	CU_ASSERT(db_node_check_for_flag(node11,PRUNED));
	CU_ASSERT(db_node_check_for_flag(node12,PRUNED));
	CU_ASSERT(db_node_check_for_flag(node13,PRUNED));
	
	//The high coverage floating path shouldn't be removed
	CU_ASSERT(!db_node_check_for_flag(node21,PRUNED));
	CU_ASSERT(!db_node_check_for_flag(node22,PRUNED));
	CU_ASSERT(!db_node_check_for_flag(node23,PRUNED));
}  



/*
 void test_db_graph_db_node_has_precisely_n_edges_with_status(){ 
 int kmer_size = 3;
 int number_of_bits = 4;
 int bucket_size = 4;
 dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
 
 boolean found;
 dBNode * node1, * node2, * node3, * node4;
 dBNode  * next_node[4];
 Orientation next_orientation[4];
 Nucleotide next_base[4];
 
 BinaryKmer tmp_kmer1;
 BinaryKmer tmp_kmer2;
 
 node1 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2), &found, db_graph);
 node2 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2), &found, db_graph);
 node3 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2), &found, db_graph);
 
 
 
 db_node_add_edge(node1, node2, forward,reverse, db_graph->kmer_size,0);
 db_node_add_edge(node1, node3, forward,forward, db_graph->kmer_size,0);
 
 boolean one_edge1 = db_graph_db_node_has_precisely_n_edges_with_status(node1,forward,ASSIGNED,1,
 next_node,next_orientation,next_base,db_graph);
 
 
 CU_ASSERT_EQUAL(one_edge1,false);
 
 db_node_action_set_flag(node2,VISITED);
 
 boolean one_edge2 = db_graph_db_node_has_precisely_n_edges_with_status(node1,forward,ASSIGNED,1,
 next_node,next_orientation,next_base,db_graph);
 
 
 CU_ASSERT_EQUAL(one_edge2,true);
 CU_ASSERT_EQUAL(node3,next_node[0]);
 CU_ASSERT_EQUAL(next_orientation[0],forward);
 CU_ASSERT_EQUAL(next_base[0],Guanine);
 
 
 node4 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2), &found, db_graph);
 db_node_add_edge(node1, node4, forward, forward, db_graph->kmer_size, 0);
 
 boolean one_edge3 = db_graph_db_node_has_precisely_n_edges_with_status(node1,forward,ASSIGNED,2,
 next_node,next_orientation,next_base,db_graph);
 
 
 
 CU_ASSERT_EQUAL(one_edge3,true);
 
 //observe that the tests below require to know in wich order the edges are VISITED
 
 CU_ASSERT_EQUAL(node4,next_node[0]);
 CU_ASSERT_EQUAL(node3,next_node[1]);
 
 CU_ASSERT_EQUAL(next_base[0],Adenine);
 CU_ASSERT_EQUAL(next_base[1],Guanine);
 
 CU_ASSERT_EQUAL(next_orientation[0],forward);
 CU_ASSERT_EQUAL(next_orientation[1],forward);
 
 hash_table_free(&db_graph);
 CU_ASSERT(db_graph == NULL);
 }
 
 */

/* TODO we don't have this method anymore, the proper way should be use the Stats program
 void test_get_N50()
 {
 
 //first set up the hash/graph
 int kmer_size = 5;
 int number_of_bits = 5;
 int bucket_size    = 4;
 int max_retries = 10;
 
 dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
 
 int seq_length=0;
 //long long count_kmers = 0;
 long long bad_reads = 0;
 
 
 //1. One fasta file containing two reads:
 //  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA which becomes just a one node supernode
 //  CACGGTATTATTTACCT which is a 13 node supernode.
 // so the N50 is 13
 
 int max_read_length=70;
 seq_length = load_fasta_from_filename_into_graph("../data/test/graph/n50_example1.fasta",0, &bad_reads, max_read_length,  db_graph);
 
 
 CU_ASSERT_EQUAL(seq_length,58);
 //CU_ASSERT_EQUAL(count_kmers,14);
 CU_ASSERT_EQUAL(bad_reads,0);
 
 CU_ASSERT(db_graph_get_N50_of_supernodes(db_graph)==13);
 
 //clean-up
 hash_table_free(&db_graph);
 
 
 */
//***************************************
// example 2:

// >read1
// AAAAAA
// >read2
// CCCCCC
// >read3
// ACGTA
// >read4
// CTATG
// >read5
// TTTAT
// >read6
// GCTTA
// >read7
// AGGCT
// >read8
// CACGGTATT
// 7 singleton super´nodes, and one 5-node supernode. So N50 is 1
/*
 //first set up the hash/graph
 db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
 
 
 //count_kmers=0;
 bad_reads=0;
 seq_length=0;
 seq_length = load_fasta_from_filename_into_graph("../data/test/graph/n50_example2.fasta",0 , &bad_reads, max_read_length,  db_graph);
 
 CU_ASSERT_EQUAL(seq_length,46);
 //  CU_ASSERT_EQUAL(count_kmers,12);
 CU_ASSERT_EQUAL(bad_reads,0);
 
 
 CU_ASSERT(db_graph_get_N50_of_supernodes(db_graph)==1); //remember this leaves the nodes all VISITED
 
 //clean-up
 hash_table_free(&db_graph);
 
 
 
 
 
 
 }
 */
/* TODO write an equivalent function in path, and make the unit test for that
 void test_is_condition_true_for_all_nodes_in_supernode()
 {
 //first set up the hash/graph
 int kmer_size = 3;
 int number_of_bits=4;
 int bucket_size   = 10;
 int seq_length;
 long long bad_reads = 0;
 
 BinaryKmer tmp_kmer1;
 BinaryKmer tmp_kmer2;
 
 
 
 dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
 
 
 //1. Sequence of tests as follows
 //         Each test loads a single specifically designed fasta file into a dB_graph.
 
 
 // ****
 //1.1 Fasta file that generate a graph with two hairpins, and a single edge (in each rorientation) joining them.
 //  Sequence is :  ACGTAC
 // ****
 
 
 seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_two_self_loops.fasta",0, s &bad_reads, 20,  db_graph);
 
 CU_ASSERT_EQUAL(seq_length,6);
 CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),2);
 CU_ASSERT_EQUAL(bad_reads,0);
 
 
 dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GTA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
 CU_ASSERT(test_element1!=NULL);
 dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("ACG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
 CU_ASSERT(test_element2!=NULL);
 
 int limit =50;
 dBNode * nodes_path[limit];
 Orientation orientations_path[limit];
 Nucleotide labels_path[limit];
 char seq[limit+1];
 int length_path;
 double avg_coverage;
 boolean is_cycle;
 int min_coverage, max_coverage;
 
 //  printf("Check if all nodes have status none - this should be true\n");
 CU_ASSERT(db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,						   
 &db_node_check_status_none,
 &db_node_action_do_nothing, 
 nodes_path, orientations_path, labels_path, &length_path,
 seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
 db_graph));
 CU_ASSERT(db_node_check_for_flag(test_element1, none));
 //  printf("set status of one node to VISITED\n");
 db_node_set_flag(test_element1, VISITED);
 CU_ASSERT(db_node_check_for_flag(test_element1, VISITED));
 CU_ASSERT(db_node_check_for_flag(test_element2, ASSIGNED));
 // printf("Checkif all nodes have status none - this should not be true - confirm this\n");
 CU_ASSERT(!db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,
 &db_node_check_status_none,
 &db_node_action_do_nothing, 
 nodes_path, orientations_path, labels_path, &length_path, 
 seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
 db_graph));
 
 
 //  printf("double check statuses of nodes unchanged\n");
 CU_ASSERT(db_node_check_for_flag(test_element1, VISITED));
 CU_ASSERT(db_node_check_for_flag(test_element2, ASSIGNED));
 
 // printf("Check if all nodes have status VISITED - this should not be true - confirm this\n");
 CU_ASSERT(!db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,
 &db_node_check_status_visited,
 &db_node_action_do_nothing, 
 nodes_path, orientations_path, labels_path, &length_path,
 seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
 db_graph));
 
 
 //check again, but this time, set all nodes to VISITED in the process
 //printf("Set all nodes to VISITED while checking them all\n");
 CU_ASSERT(!db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50, 
 &db_node_check_status_visited,
 &db_node_action_set_status_visited_or_visited_and_exists_in_reference, 
 nodes_path, orientations_path, labels_path, &length_path, 
 seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
 db_graph));
 //and now this time it SHOULD BE TRUE
 CU_ASSERT(db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50, 
 &db_node_check_status_visited,
 &db_node_action_do_nothing, 
 nodes_path, orientations_path, labels_path, &length_path, 
 seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
 db_graph));
 
 // and should still be true
 CU_ASSERT(db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,
 &db_node_check_status_visited,
 &db_node_action_do_nothing, 
 nodes_path, orientations_path, labels_path, &length_path, 
 seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
 db_graph));
 
 
 //printf("Nodes currently all VISITED. Set all nodes to none while checking them all to confirm that they are currently all VISITED (before doing the set-to-none)\n");
 CU_ASSERT(db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,
 &db_node_check_status_visited,
 &db_node_action_set_status_none, 
 nodes_path, orientations_path, labels_path, &length_path,
 seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
 db_graph));
 
 CU_ASSERT(db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,
 &db_node_check_status_none,
 &db_node_action_do_nothing, 
 nodes_path, orientations_path, labels_path, &length_path, 
 seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
 db_graph));
 
 db_node_set_flag(test_element2, PRUNED);
 CU_ASSERT(!db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,
 &db_node_check_status_none,
 &db_node_action_do_nothing, 
 nodes_path, orientations_path, labels_path, &length_path, 
 seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
 db_graph));
 
 CU_ASSERT(!db_graph_is_condition_true_for_all_nodes_in_supernode(test_element2, 50,
 &db_node_check_status_none,
 &db_node_action_do_nothing, 
 nodes_path, orientations_path, labels_path, &length_path,
 seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
 db_graph));
 
 
 hash_table_free(&db_graph);
 
 
 }
 
 */

/* TODO find where and for what this is used... 
 void test_read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference()
 {
 
 //first set up the hash/graph
 int kmer_size = 3;
 int number_of_bits=4;
 int bucket_size   = 10;
 long long bad_reads=0;
 int max_read_length=30;
 
 dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
 
 int seq_length = load_fasta_from_filename_into_graph("../data/test/graph/person.fasta", 0,  &bad_reads, max_read_length, db_graph);
 
 read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("../data/test/graph/chrom1.fasta", db_graph);
 read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("../data/test/graph/chrom2.fasta", db_graph);
 
 //Now see if it correctly gets the supernode that does not intersect a chromosome
 dBNode * nodes_path[100];
 Orientation orientations_path[100];
 Nucleotide labels_path[100];
 char seq[100];
 int length_path;
 double avg_coverage;
 boolean is_cycle;
 int min_coverage, max_coverage;
 
 BinaryKmer tmp_kmer1;
 BinaryKmer tmp_kmer2;
 
 
 //element on supernode we know intersects chromosomes
 dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("ACA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
 CU_ASSERT(test_element1!=NULL);
 //elemtn on node that does not intersect chromosomes - is "novel"
 dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
 CU_ASSERT(test_element2!=NULL);
 
 
 CU_ASSERT(!db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50, 
 &db_node_check_status_is_not_exists_in_reference,
 &db_node_action_set_status_visited, 
 nodes_path, orientations_path, labels_path, &length_path, 
 seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
 db_graph));
 
 CU_ASSERT(length_path==3);
 CU_ASSERT_STRING_EQUAL(seq, "TTC");
 
 CU_ASSERT(db_graph_is_condition_true_for_all_nodes_in_supernode(test_element2, 50,
 &db_node_check_status_is_not_exists_in_reference,
 &db_node_action_set_status_visited,
 nodes_path, orientations_path, labels_path, &length_path, 
 seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
 db_graph));
 CU_ASSERT(length_path==2);
 CU_ASSERT_STRING_EQUAL(seq, "CC");
 
 hash_table_free(&db_graph);
 
 
 
 // now a harder example.
 //first set up the hash/graph
 kmer_size = 31;
 number_of_bits=10;
 bucket_size   = 10;
 bad_reads=0;
 max_read_length=2000;
 int max_rehash_tries=10;
 
 db_graph = hash_table_new(number_of_bits,bucket_size,max_rehash_tries,kmer_size);
 
 seq_length = load_fasta_from_filename_into_graph("../data/test/graph/person2.fasta", &bad_reads, max_read_length, db_graph);
 
 read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("../data/test/graph/Homo_sapiens.NCBI36.52.dna.chromosome.1.first_20_lines.fasta", db_graph);
 
 
 //Now see if it correctly gets the supernode that does not intersect a chromosome
 
 //element on supernode we know intersects chromosomes
 test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AACCCTAACCCTAACCCTAACCCTAACCCTA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
 CU_ASSERT(test_element1!=NULL);
 //elemtn on node that does not intersect chromosomes - is "novel"
 test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("GCGGGGCGGGGCGGGGCGGGGCGGGGCCCCC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
 CU_ASSERT(test_element2!=NULL);
 
 
 length_path=0;
 CU_ASSERT(!db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,  
 &db_node_check_status_is_not_exists_in_reference,
 &db_node_action_set_status_visited, 
 nodes_path, orientations_path, labels_path, &length_path, 
 seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
 db_graph));
 
 CU_ASSERT(length_path==6);
 
 CU_ASSERT_STRING_EQUAL(seq,"ACCCTA");
 
 CU_ASSERT(db_graph_is_condition_true_for_all_nodes_in_supernode(test_element2, 50,
 &db_node_check_status_is_not_exists_in_reference,
 &db_node_action_set_status_visited, 
 nodes_path, orientations_path, labels_path, &length_path,
 seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
 db_graph));
 
 CU_ASSERT(length_path==13);
 CU_ASSERT_STRING_EQUAL(seq, "CCCTCACACACAT");
 
 
 hash_table_free(&db_graph);
 
 
 // and now another
 
 //first set up the hash/graph
 kmer_size = 31;
 number_of_bits=20;
 bucket_size   = 10;
 bad_reads=0;
 max_read_length=2000;
 max_rehash_tries=10;
 seq_length=0;
 db_graph = hash_table_new(number_of_bits,bucket_size,max_rehash_tries,kmer_size);
 
 seq_length = load_fasta_from_filename_into_graph("../data/test/graph/person3.fasta", &bad_reads, max_read_length, db_graph);
 
 
 // fasta looks like this:
 *//*
    >read1 overlaps human chrom 1 - it is 9 copies of TAACCC and then TAACC on the end.  
    TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACC
    > read 2 overlaps human chrom 1 - is has a CCCC in the middle which means any 31-mer will contain CCCC and so not overlap the above
    ACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAAC
    > read 3 does not
    GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
    > read 3 does not
    GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
    > read 3 does not
    GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
    > read 4 does not, but has too low coverage
    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
	*/
/*
 
 // This person has the following supernodes:
 
 // >this supernode is the entrety of read2 and is entirely contained in chrom1
 // ACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAAC
 // >node_1 - This supernode lies entirely in chrom1 also. However note the supernode is a loop, so you could start printing at various places....
 // TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT
 // > node 3 - overlaps chrom but has low covg
 // AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
 // >node 4 - no overlap with chromosome
 // GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
 
 
 read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("../data/test/graph/Homo_sapiens.NCBI36.52.dna.chromosome.1.first_20_lines.fasta", db_graph);
 
 char** array_of_supernodes_for_person3= (char**) calloc(10,sizeof(char*));
 array_of_supernodes_for_person3[0]= (char*)calloc(100,sizeof(char));
 array_of_supernodes_for_person3[1]= (char*)calloc(100,sizeof(char));
 array_of_supernodes_for_person3[2]= (char*)calloc(100,sizeof(char));
 
 int number_of_supernodes=0;
 int min_covg_required = 2;
 
 
 //element on supernode we know intersects chromosomes
 test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("ACCCTAACCCTAACCCTAACCCTAACCCTAA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
 CU_ASSERT(test_element1!=NULL);
 //elemtn on node that does not intersect chromosomes - is "novel" - and has coverage 3
 test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("GGGCGGGGCGGGGCGGGGCGGGGCCCCCTCA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
 CU_ASSERT(test_element2!=NULL);
 //element does not intersect chrom but has covg only 1
 dBNode* test_element3 =  hash_table_find(element_get_key(seq_to_binary_kmer("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
 CU_ASSERT(test_element3!=NULL);
 
 
 
 
 
 db_graph_print_supernodes_where_condition_is_true_for_all_nodes_in_supernode(db_graph, &db_node_check_status_is_not_exists_in_reference, min_covg_required, NULL,
 true, array_of_supernodes_for_person3, &number_of_supernodes);
 
 
 CU_ASSERT(number_of_supernodes==1);
 CU_ASSERT( !strcmp(array_of_supernodes_for_person3[0], "CCCGCCCCGCCCC")
 || !strcmp(array_of_supernodes_for_person3[0], "CCCTCACACACAT"));
 
 
 hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);
 
 
 number_of_supernodes=0;
 //here just a quick check of the _at_least_ version of the function
 min_covg_required=1;
 db_graph_print_supernodes_where_condition_is_true_for_at_least_one_node_in_supernode(db_graph,
 &db_node_check_status_exists_in_reference,  min_covg_required, 
 NULL, true, array_of_supernodes_for_person3, 
 &number_of_supernodes);
 
 
 CU_ASSERT(number_of_supernodes==2);
 
 
 //some of read 1 is in chrom1, and the supernode is printed above, just after we loaded the fasta
 //the whole of read2 is a supernode, and is in chrom 1. Note this print function prints only the edges, not the first kmer in the path
 
 
 //one of the supernodes is unambiguous
 CU_ASSERT(    (!strcmp(array_of_supernodes_for_person3[0],"ACCCTAACCCTAAC")) 
 || (!strcmp(array_of_supernodes_for_person3[1],"ACCCTAACCCTAAC"))
 || (!strcmp(array_of_supernodes_for_person3[0],"GTTAGGGTTAGGGT")) 
 || (!strcmp(array_of_supernodes_for_person3[1],"GTTAGGGTTAGGGT")) 
 );
 
 
 //the other one is ambiguous because you could start printing at one of 5 places (is a loop of 5 kmers)
 */
/* - this is what we test for if the print prints the whole supernode including the first 31 bases
 CU_ASSERT( 
 !strcmp(array_of_supernodes_for_person3[0],"TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT") || !strcmp(array_of_supernodes_for_person3[1],"TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA") || !strcmp(array_of_supernodes_for_person3[1],"AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA") || !strcmp(array_of_supernodes_for_person3[1],"AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"TAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTT") || !strcmp(array_of_supernodes_for_person3[1],"TAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTT")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"ACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA") || !strcmp(array_of_supernodes_for_person3[1],"ACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGT") || !strcmp(array_of_supernodes_for_person3[1],"TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGT")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC") || !strcmp(array_of_supernodes_for_person3[1],"CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"GTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG") || !strcmp(array_of_supernodes_for_person3[1],"GTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACC") || !strcmp(array_of_supernodes_for_person3[1],"CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACC")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"GGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGG") || !strcmp(array_of_supernodes_for_person3[1],"GGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGG")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC") || !strcmp(array_of_supernodes_for_person3[1],"CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"GGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAG") || !strcmp(array_of_supernodes_for_person3[1],"GGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAG")
 
 );
 */

/*
 
 CU_ASSERT( !strcmp(array_of_supernodes_for_person3[0],"AACCCT") || !strcmp(array_of_supernodes_for_person3[1],"AACCCT")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"AGGGTT") || !strcmp(array_of_supernodes_for_person3[1],"AGGGTT")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"ACCCTA") || !strcmp(array_of_supernodes_for_person3[1],"ACCCTA")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"TAGGGT") || !strcmp(array_of_supernodes_for_person3[1],"TAGGGT")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"CCCTAA") || !strcmp(array_of_supernodes_for_person3[1],"CCCTAA")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"TTAGGG") || !strcmp(array_of_supernodes_for_person3[1],"TTAGGG")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"CCTAAC") || !strcmp(array_of_supernodes_for_person3[1],"CCTAAC")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"GTTAGG") || !strcmp(array_of_supernodes_for_person3[1],"GTTAGG")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"CTAACC") || !strcmp(array_of_supernodes_for_person3[1],"CTAACC")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"GGTTAG") || !strcmp(array_of_supernodes_for_person3[1],"GGTTAG")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"TAACCC") || !strcmp(array_of_supernodes_for_person3[1],"TAACCC")
 ||
 !strcmp(array_of_supernodes_for_person3[0],"GGGTTA") || !strcmp(array_of_supernodes_for_person3[1],"GGGTTA")
 
 
 );
 
 
 free(array_of_supernodes_for_person3[0]) ;
 free(array_of_supernodes_for_person3[1]) ;
 free(array_of_supernodes_for_person3[2]) ;
 free(array_of_supernodes_for_person3) ;
 
 hash_table_free(&db_graph);
 
 
 
 
 }
 
 */


/*TODO make the unit test  using the colors
 void test_indel_discovery_simple_test_1()
 {
 
 
 //first set up the hash/graph
 int kmer_size = 31;
 int number_of_bits=10;
 int bucket_size   = 5;
 long long bad_reads=0;
 int max_read_length=600;
 
 dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
 
 int seq_length = load_fasta_from_filename_into_graph("../data/test/graph/person_with_sv.fasta", &bad_reads, max_read_length, db_graph);
 CU_ASSERT(seq_length==612);
 
 read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("../data/test/graph/Homo_sapiens.NCBI36.52.dna.chromosome.1.first_20_lines.fasta", db_graph);
 
 BinaryKmer tmp_kmer1;
 BinaryKmer tmp_kmer2;
 
 
 //element on supernode we know intersects chromosome entirely
 dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("TGTGCAGAGGACAACGCAGCTCCGCCCTCGC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
 CU_ASSERT(test_element1!=NULL);
 
 //element on supernode that overlaps at start and end but not middle, with chromosome
 dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("CCGGCGCAGGCGCAGTTGTTGTAGAGGCGCG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
 CU_ASSERT(test_element2!=NULL);
 
 
 dBNode * nodes_path[100];
 Orientation orientations_path[100];
 Nucleotide labels_path[100];
 char seq[100];
 int length_path;
 int min_diff=1; //enought o have one node different in the middle
 double avg_coverage;
 int min_coverage, max_coverage;
 boolean is_cycle;
 
 
 int num_nodes_we_demand_overlap_with_reference_at_start=2;
 int num_nodes_we_demand_overlap_with_reference_at_end=22;
 
 CU_ASSERT(!db_graph_is_condition_true_for_start_and_end_but_not_all_nodes_in_supernode(test_element1, 200,
 &db_node_check_status_exists_in_reference,
 &db_node_action_set_status_visited_or_visited_and_exists_in_reference,
 num_nodes_we_demand_overlap_with_reference_at_start,
 num_nodes_we_demand_overlap_with_reference_at_end, min_diff, 
 nodes_path, orientations_path, labels_path, &length_path, 
 seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
 db_graph));
 
 CU_ASSERT(db_graph_is_condition_true_for_start_and_end_but_not_all_nodes_in_supernode(test_element2, 200, 
 &db_node_check_status_exists_in_reference,
 
 &db_node_action_set_status_visited_or_visited_and_exists_in_reference,
 num_nodes_we_demand_overlap_with_reference_at_start,
 num_nodes_we_demand_overlap_with_reference_at_end, min_diff, 
 nodes_path, orientations_path, labels_path, &length_path, 
 seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
 db_graph));
 
 
 //remove the VISITED markings
 hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);
 
 
 char** array_of_supernodes= (char**) calloc(10,sizeof(char*));
 array_of_supernodes[0]= (char*)calloc(100,sizeof(char));
 array_of_supernodes[1]= (char*)calloc(100,sizeof(char));
 array_of_supernodes[2]= (char*)calloc(100,sizeof(char));
 
 int number_of_supernodes=0;
 int min_covg_required = 2;
 int min_start = 2;
 int min_end = 2;
 
 db_graph_print_supernodes_where_condition_is_true_at_start_and_end_but_not_all_nodes_in_supernode(db_graph, &db_node_check_status_exists_in_reference, min_covg_required,
 min_start, min_end, min_diff, NULL,
 true, array_of_supernodes, &number_of_supernodes);
 
 
 
 
 CU_ASSERT(number_of_supernodes==1);
 
 CU_ASSERT( !strcmp(array_of_supernodes[0], "AGTTGTTGTAGAGGCGCGCCGCGCCGGCGCAGGCGCAGACACATGCTAGCGCGTCGGGGTGGAGGCGT") || !strcmp(array_of_supernodes[0], "TGCGCCTGCGCCGGCGCGGCGCGCCTCTACAACAACTGCGCCTGCGCCGGCGGAGTTGCGTTCTCCTC"));
 //  CU_ASSERT( !strcmp(array_of_supernodes[0], "GAGGAGAACGCAACTCCGCCGGCGCAGGCGCAGTTGTTGTAGAGGCGCGCCGCGCCGGCGCAGGCGCAGACACATGCTAGCGCGTCGGGGTGGAGGCGT") || !strcmp(array_of_supernodes[0], "ACGCCTCCACCCCGACGCGCTAGCATGTGTCTGCGCCTGCGCCGGCGCGGCGCGCCTCTACAACAACTGCGCCTGCGCCGGCGGAGTTGCGTTCTCCTC"));
 
 
 free(array_of_supernodes[0]) ;
 free(array_of_supernodes[1]) ;
 free(array_of_supernodes[2]) ;
 free(array_of_supernodes) ;
 
 
 
 hash_table_free(&db_graph);
 }
 
 */

/*
 //suppose we are given a pair of start-end coords within which we think there is a deletion.
 // Load your fasta, then load precisely that section of the chromosome as reference, and print supernodes that match ref at start and end
 // (WARNING: if you delete the middle of a sequence, you create new kmers in the middle that might not have been there before.
 //           so the nodes in the new supernode are NOT necessatrily all in the reference)
 // Then build a new graph just out of the section of the chromosome, and load the supernodes above as reference. Now look for supernodes
 // in the chromosome graph that match the "reference" supernodes at the start and end only. These are potential deletions.
 void test_deletion_validation()
 {
 
 //first set up the hash/graph
 int kmer_size = 31;
 int number_of_bits=10;
 int bucket_size   = 5;
 long long bad_reads=0;
 int max_read_length=2000;
 
 
 
 //STEP 1: get supernodes from our person which match reference at start and end.
 dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
 int seq_length = load_fasta_from_filename_into_graph("../data/test/graph/person_with_deletion_in_chrom.fasta",  &bad_reads, max_read_length, db_graph);
 read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("../data/test/graph/Homo_sapiens.NCBI36.52.dna.chromosome.1.first_20_lines.fasta", db_graph);
 
 
 FILE* intermediate_output = fopen("../data/test/graph/test_db_graph_intermediate_output_file", "w");
 int min_covg_required = 1;
 int min_start = 1;
 int min_end = 31; //iei 31 bases at start and end
 int min_diff = 2;
 db_graph_print_supernodes_where_condition_is_true_for_at_least_one_node_in_supernode(db_graph, &db_node_check_status_exists_in_reference, min_covg_required,
 intermediate_output, false, NULL, 0);
 
 //  db_graph_print_supernodes_where_condition_is_true_at_start_and_end_but_not_all_nodes_in_supernode(db_graph, &db_node_check_status_exists_in_reference, min_covg_required,
 //                                                                                                  min_start, min_end, min_diff, intermediate_output,
 //												    false, NULL, 0);
 fclose(intermediate_output);
 hash_table_free(&db_graph);
 
 //STEP 2: Load the chromosome as a person, and the previosu supernodes as reference
 db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
 seq_length = load_fasta_from_filename_into_graph("../data/test/graph/Homo_sapiens.NCBI36.52.dna.chromosome.1.first_20_lines.fasta", &bad_reads, max_read_length,db_graph);
 read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("../data/test/graph/test_d_graph_intermediate_output_file", db_graph);
 
 
 
 char** array_of_supernodes= (char**) calloc(10,sizeof(char*));
 array_of_supernodes[0]= (char*)calloc(100,sizeof(char));
 array_of_supernodes[1]= (char*)calloc(100,sizeof(char));
 array_of_supernodes[2]= (char*)calloc(100,sizeof(char));
 
 int number_of_supernodes=0;
 
 //STEP 3: Print supernodes in chromosome that match our person's supernodes (here the reference) at start and end - deletions
 db_graph_print_supernodes_where_condition_is_true_at_start_and_end_but_not_all_nodes_in_supernode(db_graph, &db_node_check_status_exists_in_reference, min_covg_required,
 min_start, min_end, min_diff, NULL,
 true, array_of_supernodes, &number_of_supernodes);
 //db_graph_print_supernodes_where_condition_is_true_for_at_least_one_node_in_supernode(db_graph, &db_node_check_status_exists_in_reference, min_covg_required,
 //										       NULL,true, array_of_supernodes, &number_of_supernodes);
 
 
 
 printf("deletion %s and number fo supernodes is %d\n", array_of_supernodes[0], number_of_supernodes);
 CU_ASSERT_STRING_EQUAL("GAGAGGCGCGCCGCGCCGGCGC", array_of_supernodes[0]);
 
 free(array_of_supernodes[0]);
 free(array_of_supernodes[1]);
 free(array_of_supernodes[2]);
 free(array_of_supernodes);
 hash_table_free(&db_graph);
 
 }
 */




void test_y_walk(){
	int kmer_size = 9;
	int number_of_bits = 10; 
	int bucket_size = 30;
	long long bad_reads = 0; 
	int seq_length;
	dBGraph * db_graph;
	BinaryKmer tmp_kmer1, tmp_kmer2;
	
	path_array_destroy_buffers();
	path_array_initialise_buffers(kmer_size);
	//
	
	db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
	
	seq_length = load_fastq_from_filename_into_graph("../data/test/read_pair/rp_1.fastq",0, &bad_reads, 5, 200, 63, db_graph);
	seq_length = load_fastq_from_filename_into_graph("../data/test/read_pair/rp_2.fastq",0, &bad_reads, 5, 200, 63, db_graph);
	db_graph_write_graphviz_file("test_y_walk.gv", db_graph);
	dBNode * a = hash_table_find(element_get_key(seq_to_binary_kmer("ATTTGACGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * b = hash_table_find(element_get_key(seq_to_binary_kmer("CCATCGGCA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * c = hash_table_find(element_get_key(seq_to_binary_kmer("CCGCTGCCG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * d = hash_table_find(element_get_key(seq_to_binary_kmer("GCATCACCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * e = hash_table_find(element_get_key(seq_to_binary_kmer("CCATCACCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * f = hash_table_find(element_get_key(seq_to_binary_kmer("AATGCCACC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * g = hash_table_find(element_get_key(seq_to_binary_kmer("AATACCACC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * h = hash_table_find(element_get_key(seq_to_binary_kmer("AATCCCACC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * i = hash_table_find(element_get_key(seq_to_binary_kmer("AACAACTTA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * j = hash_table_find(element_get_key(seq_to_binary_kmer("AAGTTACCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);//TODO: validate that I and J are the o
	dBNode * k = hash_table_find(element_get_key(seq_to_binary_kmer("GCATCACCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	
	pathStep start_step[7];
	pathStep end_step[7];
	
	//Unconnected path
	start_step[0].node = i;
	end_step  [0].node = j;
	
	
	//b is the y node
	start_step[1].node = e;
	end_step	[1].node = b;
	
	start_step[2].node = d;
	end_step	[2].node = b;
	
	
	//c is the y node
	start_step[3].node = h;
	end_step	[3].node = c;
	start_step[4].node = g;
	end_step	[4].node = c;
	start_step[5].node = f;
	end_step	[5].node = c;
	start_step[6].node = k;
	end_step	[6].node = c;
	
	
    //Unconnected path
	start_step[0].orientation = reverse;//i
	end_step	[0].orientation = reverse;//j
	
	
	//b is the y node
	start_step[1].orientation = reverse;//e
	end_step  [1].orientation = reverse;//a
	start_step[2].orientation = reverse;//d
	end_step  [2].orientation = reverse;//a
	
	
	//c is the y node
	start_step[3].orientation = forward;//h
	end_step	[3].orientation = forward;//e
	start_step[4].orientation = forward;//g
	end_step	[4].orientation = reverse;//d
	start_step[5].orientation = forward;//f
	end_step	[5].orientation = forward;//c
	start_step[6].orientation = forward;//k
	end_step	[6].orientation = reverse;//d
	dBNode * double_y_node[2];
	
	double_y_node[0] = b;
	double_y_node[1] = c;
	
	Path * tmp_path = path_new(100, db_graph->kmer_size);
	pathStep  ps_last, ps_first;
	int z;
	
	for(z = 0; z < 7; z++){
		printf("___________________________________________TEST %d_____________________________________________________________\n", z);
		path_reset(tmp_path);
		
		if(z < 3){
			db_node_action_set_flag(b, Y_START);
			db_node_action_unset_flag(c, Y_START);
		}else{
			db_node_action_set_flag(c, Y_START);
			db_node_action_unset_flag(b, Y_START);
		}
		
		// y_walk_get_path()
		 y_walk_get_path(start_step[z].node, reverse,  &db_node_action_set_flag_visited, db_graph, false, tmp_path);
		//y_walk_get_path(dBNode * node, Orientation orientation,
		//      void (*node_action) (dBNode * node),
		//     dBGraph * db_graph, boolean both_directions,Path * path)
		//boolean read_pair_get_path(dBNode * node, void (*node_action) (dBNode * node),
		//				 ReadPairDescriptorArray * rpda, dBGraph * db_graph, Path * path)
		
		path_get_last_step(&ps_last, tmp_path);
		path_get_step_at_index(0, &ps_first, tmp_path);
		path_to_fasta(tmp_path, stdout);
		
		switch (z) {
			case 1:
				CU_ASSERT(path_step_equals_without_label(&ps_first, &start_step[1])==true);
				CU_ASSERT(path_step_equals_without_label(&ps_last, &end_step[1])==true);
				break;
			case 0:
				CU_ASSERT(path_step_equals_without_label(&ps_last, &end_step[0])==true);
				CU_ASSERT(path_step_equals_without_label(&ps_first, &start_step[0])==true);
				break;
			case 2:
				CU_ASSERT(path_step_equals_without_label(&ps_first, &start_step[2])==true);
				CU_ASSERT(path_step_equals_without_label(&ps_last, &end_step[2])==true);
				break;
			case 3:
				CU_ASSERT(ps_first.node == start_step[3].node);
				CU_ASSERT(ps_last.node == end_step[3].node);
				break;
			case 4:
				CU_ASSERT((ps_first.node == start_step[4].node)==true);
				CU_ASSERT((ps_last.node == end_step[4].node)==true);
				break;
			case 5:
				CU_ASSERT((ps_first.node == start_step[5].node)==true);
				CU_ASSERT((ps_last.node == end_step[5].node)==true);
				break;
			case 6:
				CU_ASSERT((ps_first.node == start_step[6].node)==true);
				CU_ASSERT((ps_last.node == end_step[6].node)==true);
				break;
			default:
				break;
		}
		
		
	}
	
	hash_table_free(&db_graph);
}

void test_y_walk_from_perfect_path_tests() //test db_graph_get_perfect_path
{
	
	//first set up the hash/graph
	int kmer_size = 3;
	int number_of_bits=5;
	int bucket_size   = 10;
	int seq_length;
	long long bad_reads = 0;
	/*dBNode * path_nodes[100];
	 Orientation path_orientations[100];
	 Nucleotide path_labels[100];
	 char tmp_seq[100];
	 
	 double  avg_coverage;
	 int min_coverage, max_coverage;
	 */
	dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
	
	
	//1. Sequence of tests as follows
	//         Each test loads a single specifically designed fasta file into a dB_graph.
	//         The test then picks an element in the graph, and calls get_perfect_path
	//         and checks that it gets the right sequence.
	
	
	// ****
	//1.1 Fasta file that generate a graph with two hairpins, and a single edge (in each rorientation) joining them.
	//  Sequence is :  ACGTAC
	// ****
	
#ifdef ENABLE_READ_PAIR_OLD
	seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_two_self_loops.fasta",0,0,&bad_reads, 20, db_graph);
#else
	seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_two_self_loops.fasta",0,  &bad_reads, 20,  db_graph);
#endif
	CU_ASSERT_EQUAL(seq_length,6);
	CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),2);
	CU_ASSERT_EQUAL(bad_reads,0);
	
	
	// now start at GTA and get all the sequence from there to the end of the supernode, and see
	// if that is right.
	
	BinaryKmer tmp_kmer1;
	BinaryKmer tmp_kmer2;
	
	dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GTA", kmer_size,&tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
	CU_ASSERT(test_element1!=NULL);
	
	boolean is_cycle=false;
	
	Path * tmp_path = path_new(100, db_graph->kmer_size);
	
//	int test1_length =  db_graph_supernode(test_element1, &db_node_action_set_flag_visited,
	//									   tmp_path, db_graph);
										   
	
	db_graph_write_graphviz_file("test1.viz", db_graph);
	y_walk_get_path(test_element1, forward,  &db_node_action_set_flag_visited, db_graph, false, tmp_path);
	int test1_length = path_get_length(tmp_path);
	
	
    is_cycle = path_is_cycle(tmp_path);
    CU_ASSERT(is_cycle);
    CU_ASSERT_EQUAL(test1_length,4);
	
    CU_ASSERT_STRING_EQUAL(tmp_path->seq,"CGTA");
	
	hash_table_free(&db_graph);
	CU_ASSERT(db_graph == NULL);
	
	
	/*   // **** */
	/*   //1.2 Fasta file that generate a graph with one long supernode, with a conflict at the end */
	/*   //   caused by two outward/exiting edges */
	/*   // **** */
	
	//first set up the hash/graph
	kmer_size = 3;
	number_of_bits= 4;
	bucket_size   = 10;
	
	bad_reads = 0;
	dBGraph * db_graph2 = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
	//db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
#ifdef ENABLE_READ_PAIR_OLD
	seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_one_long_supernode_with_conflict_at_end.fasta",0,0,&bad_reads, 20, db_graph2);
#else
	seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_one_long_supernode_with_conflict_at_end.fasta",0, &bad_reads,20,db_graph2);
#endif
	db_graph_write_graphviz_file("test2.viz", db_graph2);
	CU_ASSERT_EQUAL(seq_length,13);
	CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph2),5);
	CU_ASSERT_EQUAL(bad_reads,0);
	
	test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("ACA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph2);
	CU_ASSERT(test_element1!=NULL);
	char seq[20];
	CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(element_get_kmer(test_element1),3,&seq[0]),"ACA");
	
	//ACA < TGT so forward gives TTG, none is marked as double Y. 
	is_cycle=false;
	path_reset(tmp_path);
	
	//go forward
	y_walk_get_path(test_element1, forward,  &db_node_action_do_nothing, db_graph2, false, tmp_path);	
   int test2_length = path_get_length(tmp_path);
    
    //int test2_length = db_graph_get_perfect_path(test_element1,undefined,&db_node_action_do_nothing,db_graph2, tmp_path);
	
	
	//fprintf(stderr, "Expected: TT obtained: %s, len: %d", tmp_path->seq, test2_length);
	CU_ASSERT(!path_is_cycle(tmp_path));
	CU_ASSERT_EQUAL(test2_length,4);
	CU_ASSERT_STRING_EQUAL(tmp_path->seq,"TTG");
	
	
	//ACA < TGT so backward gives ""
	
	//go backward
	
    y_walk_get_path(test_element1, reverse,  &db_node_action_do_nothing, db_graph2, false, tmp_path);	
     test2_length = path_get_length(tmp_path);
	/*test2_length = db_graph_get_perfect_path(test_element1,reverse,100,
	 &db_node_action_do_nothing,
	 path_nodes,path_orientations,path_labels,
	 tmp_seq, &avg_coverage, &min_coverage, &max_coverage,
	 &is_cycle,db_graph);
	 */
	
	CU_ASSERT(!path_is_cycle(tmp_path));
	CU_ASSERT_EQUAL(test2_length,0);
	
	
	CU_ASSERT_STRING_EQUAL(tmp_path->seq,"");
	printf("We got as seq: %s\n",  tmp_path->seq);
	hash_table_free(&db_graph2);
	CU_ASSERT(db_graph2 == NULL);
	
	
	// ****
	//1.3 Fasta file that generate a graph with one long supernode, with a conflict at the end
	//   caused by two INWARD edges in the opposite direction
	// ****
	
	//first set up the hash/graph
	kmer_size = 3;
	number_of_bits = 4;
	bucket_size    = 5;
	bad_reads = 0;
	
	db_graph = hash_table_new(number_of_bits, bucket_size,10,kmer_size);
	
#ifdef ENABLE_READ_PAIR_OLD
	seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_one_long_supernode_with_inward_conflict_at_end.fasta",0,0,&bad_reads, 20, db_graph);
#else
	seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_one_long_supernode_with_inward_conflict_at_end.fasta",0,  &bad_reads, 20,  db_graph);
#endif
	db_graph_write_graphviz_file("test3.viz", db_graph);
	CU_ASSERT_EQUAL(seq_length,13);
	CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),5);
	CU_ASSERT_EQUAL(bad_reads,0);
	
	
	test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("ACA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	CU_ASSERT(test_element1!=NULL);
	
	//ACA < TGT so forward gives TT
	//go forward
	is_cycle=false;
	y_walk_get_path(test_element1, forward,  &db_node_action_do_nothing, db_graph, false, tmp_path);
	int test3_length = path_get_length(tmp_path);
	
	CU_ASSERT(!path_is_cycle(tmp_path));
	
	CU_ASSERT_EQUAL(test3_length, 4);
	CU_ASSERT_STRING_EQUAL(tmp_path->seq,"TTC");
	
	
	//ACA < TGT so backwar gives ""
	is_cycle=false;
	
	//go reverse
	test3_length = db_graph_get_perfect_path(test_element1,reverse,
											 &db_node_action_do_nothing, db_graph, tmp_path);
	
	CU_ASSERT(!path_is_cycle(tmp_path));
	
	CU_ASSERT_EQUAL(test3_length,0);
	CU_ASSERT_STRING_EQUAL(tmp_path->seq,"");
	
	hash_table_free(&db_graph);
	CU_ASSERT(db_graph == NULL);
	
	
	// ****
	//1.4 Fasta file that generate a graph with an infinite loop at a single kmer
	//
	// ****
	
	
	//first set up the hash/graph
	kmer_size = 3;
	number_of_bits=8;
	bucket_size   =4;
	
	bad_reads = 0;
	
	db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
#ifdef ENABLE_READ_PAIR_OLD
	seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_infinite_loop.fasta",0,0,&bad_reads, 20, db_graph);
#else
	seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_infinite_loop.fasta",0, &bad_reads,30,db_graph);
#endif

	tmp_path->kmer_size = kmer_size;
	
	
	db_graph_write_graphviz_file("test4.viz", db_graph);
	CU_ASSERT_EQUAL(seq_length,25);
	CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),1);
	CU_ASSERT_EQUAL(bad_reads,0);
	
	test_element1 = hash_table_find(seq_to_binary_kmer("AAA", kmer_size, &tmp_kmer1) ,db_graph);
	CU_ASSERT(test_element1!=NULL);
	
	//forward
	is_cycle=false;
	
	
	y_walk_get_path(test_element1, forward,  &db_node_action_do_nothing, db_graph, false, tmp_path);
	int test4_length = path_get_length(tmp_path);
	
	path_to_fasta(tmp_path,stdout);
	CU_ASSERT(path_is_cycle(tmp_path));//this time it should find a cycle
	CU_ASSERT_STRING_EQUAL(tmp_path->seq,"A");
	CU_ASSERT_EQUAL(test4_length,1);
	
	//backward
	is_cycle=false;
	
	y_walk_get_path(test_element1, reverse,  &db_node_action_do_nothing, db_graph, false, tmp_path);
	 test4_length = path_get_length(tmp_path);
	
	CU_ASSERT(path_is_cycle(tmp_path));//this time it should find a cycle
	CU_ASSERT_EQUAL(test4_length,1);
	CU_ASSERT_STRING_EQUAL(tmp_path->seq,"T");
	
	hash_table_free(&db_graph);
	CU_ASSERT(db_graph == NULL);
	
	
	
	// ****
	// 1.5 check parameters (path nodes,labels,etc) for get_perfect_path
	//
	// ****
	
	kmer_size = 3;
	number_of_bits = 4;
	bucket_size = 4;
	db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
	
	boolean found1, found2, found3;
	dBNode * node1;
	dBNode * node2;
	dBNode * node3;
	dBNode * node4;
	
	
	node1 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CGT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2), &found1, db_graph);
	node2 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("GTT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found2, db_graph);
	
	db_node_add_edge(node1, node2, reverse,reverse, db_graph->kmer_size, 0);
	
	node3 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("TTA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found3, db_graph);
	db_node_add_edge(node2, node3, reverse,reverse, db_graph->kmer_size, 0);
	
	node4 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("TAG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found3, db_graph);
	db_node_add_edge(node3, node4, reverse,reverse, db_graph->kmer_size, 0);
	
	
	//add coverage
	element_update_coverage(node1,0,1);
	element_update_coverage(node2,0,2);
	element_update_coverage(node3,0,1);
	element_update_coverage(node4,0,2);
	
	
	
	/* dBNode * nodes[10];
	 Orientation orientations[10];
	 Nucleotide bases[10];*/
	int test5_length = 0;
	
	
	y_walk_get_path(node1, reverse,  &db_node_action_set_flag_visited, db_graph, false, tmp_path);
	test5_length = path_get_length(tmp_path);
	
	CU_ASSERT_EQUAL(test5_length,4);
	
	//check nodes
	CU_ASSERT_EQUAL(node1,tmp_path->nodes[0]);
	CU_ASSERT_EQUAL(node2,tmp_path->nodes[1]);
	CU_ASSERT_EQUAL(node3,tmp_path->nodes[2]);
	CU_ASSERT_EQUAL(node4,tmp_path->nodes[3]);
	
	//check labels
	CU_ASSERT_EQUAL(tmp_path->labels[0],Thymine);
	CU_ASSERT_EQUAL(tmp_path->labels[1],Adenine);
	CU_ASSERT_EQUAL(tmp_path->labels[2],Guanine);
	
	//check orientations
	CU_ASSERT_EQUAL(tmp_path->orientations[0],reverse);
	CU_ASSERT_EQUAL(tmp_path->orientations[1],reverse);
	CU_ASSERT_EQUAL(tmp_path->orientations[2],reverse);
	
	//check statuses
	CU_ASSERT(db_node_check_for_flag(tmp_path->nodes[0], VISITED));
	CU_ASSERT(db_node_check_for_flag(tmp_path->nodes[1], VISITED));
	CU_ASSERT(db_node_check_for_flag(tmp_path->nodes[2], VISITED));
	CU_ASSERT(db_node_check_for_flag(tmp_path->nodes[3], VISITED));
	
	double avg;
	int min, max;
	
	path_get_statistics(&avg, &min, &max, tmp_path);
	
	
	db_graph_write_graphviz_file("Test1_5.viz", db_graph);
	
	
	//check coverage
	CU_ASSERT_EQUAL(avg,1.5);
	CU_ASSERT_EQUAL(min,1);
	CU_ASSERT_EQUAL(max,2);
	
	
	hash_table_free(&db_graph);
	CU_ASSERT(db_graph == NULL);
}




void test_get_all_paths(){
	
	int kmer_size = 9;
	int number_of_bits = 10; 
	int bucket_size = 30;
	long long bad_reads = 0; 
	int seq_length;
	dBGraph * db_graph;
	BinaryKmer tmp_kmer1,tmp_kmer2;
	
	path_array_initialise_buffers(kmer_size);
	
	db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
	
	seq_length = load_fastq_from_filename_into_graph("../data/test/read_pair/rp_1.fastq",0, &bad_reads, 5, 200, 63, db_graph);
	seq_length = load_fastq_from_filename_into_graph("../data/test/read_pair/rp_2.fastq",0, &bad_reads, 5, 200, 63, db_graph);
	
	//y_walk_
	
	dBNode* first_node = hash_table_find(element_get_key(seq_to_binary_kmer("ATTTGACGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	
	CU_ASSERT(first_node != NULL);
	
	PathArray * pa = path_array_get_from_buffer_with_size(4);
	
	perfect_path_get_all_paths_from(first_node, forward, pa, 1000, db_graph);
	
	
	CU_ASSERT_STRING_EQUAL(path_array_get(0, pa)->seq,"");			//A  
	CU_ASSERT_STRING_EQUAL(path_array_get(1, pa)->seq,"CATCGGCA");	//C
	CU_ASSERT_STRING_EQUAL(path_array_get(2, pa)->seq,"");			//G
	CU_ASSERT_STRING_EQUAL(path_array_get(3, pa)->seq,"");			//T

	
	//CU_ASS
	path_array_free_from_buffer(pa);
	pa = path_array_get_from_buffer_with_size(4);
	perfect_path_get_all_paths_from(first_node, reverse, pa, 1000, db_graph);
	int i;
	for(i = 0; i < 4; i++){
	 path_array_get(i, pa)->kmer_size =  kmer_size;
	 path_to_fasta(path_array_get(i, pa),stderr);
	 }
	
	
	CU_ASSERT_STRING_EQUAL(path_array_get(0, pa)->seq,"ACCACC");			//A  
	CU_ASSERT_STRING_EQUAL(path_array_get(1, pa)->seq,"CCCACC");			//C
	CU_ASSERT_STRING_EQUAL(path_array_get(2, pa)->seq,"GCCACC");			//G
	CU_ASSERT_STRING_EQUAL(path_array_get(3, pa)->seq,"");			        //T
	
	
	
	path_array_free_from_buffer(pa);
	
	db_graph_write_graphviz_file("test_read_pair.viz", db_graph);
	
	
	hash_table_free(&db_graph);
	
}

//PathArray * branches_get_all_paths_from(pathStep * first, int max_length, dBGraph * db_graph)

void test_get_all_paths_deep_search(){
	
	int kmer_size = 9;
	int number_of_bits = 10; 
	int bucket_size = 30;
	long long bad_reads = 0; 
	int seq_length;
	dBGraph * db_graph;
	BinaryKmer tmp_kmer1,tmp_kmer2;
	
	path_array_initialise_buffers(kmer_size);
	
	db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
	
	seq_length = load_fastq_from_filename_into_graph("../data/test/read_pair/rp_1.fastq",0, &bad_reads, 5, 200, 63, db_graph);
	seq_length = load_fastq_from_filename_into_graph("../data/test/read_pair/rp_2.fastq",0, &bad_reads, 5, 200, 63, db_graph);
	
	//This kmer should produce 4 paths
    
	dBNode* first_node = hash_table_find(element_get_key(seq_to_binary_kmer("GCATCACCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
    
    
    dBNode* last_node = hash_table_find(element_get_key(seq_to_binary_kmer("AAATCCCAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	
    
    
    pathStep first_step;
    first_step.node = first_node;
    first_step.orientation = reverse;
    first_step.orientation = undefined;
	CU_ASSERT(first_node != NULL);
	
	//PathArray * pa = path_array_get_from_buffer_with_size(4);
	
	//perfect_path_get_all_paths_from(first_node, forward, pa, 1000, db_graph);
	
	PathArray * pa = branches_get_all_paths_from(&first_step, 25, db_graph);
    
    //path_array_to_fasta(stdout, pa);
    
	CU_ASSERT_STRING_EQUAL(path_array_get(0, pa)->seq,"CCGCTGCCGATGGCGTCAAATACCA");			//A  
	CU_ASSERT_STRING_EQUAL(path_array_get(1, pa)->seq,"CCGCTGCCGATGGCGTCAAATCCCA");	//C
	CU_ASSERT_STRING_EQUAL(path_array_get(2, pa)->seq,"CCGCTGCCGATGGCGTCAAATGCCA");			//G
	CU_ASSERT_STRING_EQUAL(path_array_get(3, pa)->seq,"CCGCTGCCGATGGGGTCAA");			//T
    
    Path * seek = path_get_buffer_path();
    pathStep last_step;
    last_step.node = last_node;
    last_step.orientation = reverse;
    last_step.orientation = undefined;
	CU_ASSERT(last_node != NULL);
    
    CU_ASSERT(branches_get_path_between(&first_step, &last_step, 100,4, seek,   db_graph) != NULL);
    CU_ASSERT_STRING_EQUAL(seek->seq,"CCGCTGCCGATGGCGTCAAATCCCAC");	     //C
    path_to_fasta(seek, stdout);
    
    path_reset(seek);
    CU_ASSERT(branches_get_path_between(&first_step, &last_step, 20,4 ,seek,   db_graph) == NULL);
    path_to_fasta(seek, stdout);
    path_free_buffer_path(seek);		
	path_array_free_from_buffer(pa);
	hash_table_free(&db_graph);
	
}

#ifdef ENABLE_READ_PAIR_OLD
void test_read_pair()
{
	int kmer_size = 9;
	int number_of_bits = 10; 
	int bucket_size = 30;
	long long bad_reads = 0; 
	int seq_length;
	dBGraph * db_graph;
	BinaryKmer tmp_kmer1, tmp_kmer2;
	
	path_array_destroy_buffers();
	path_array_initialise_buffers(kmer_size);
	//
	
	db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
	
	seq_length = load_fastq_from_filename_into_graph("../data/test/read_pair/rp_1.fastq",0, 0,&bad_reads, 5, 200, db_graph);
	seq_length = load_fastq_from_filename_into_graph("../data/test/read_pair/rp_2.fastq",0, 1,&bad_reads, 5, 200, db_graph);
	
	
	
	db_graph_write_graphviz_file("test_read_pair.gv", db_graph);
	//So far, we had just been reading the graph. From here, we are naming the end nodes of the different walks and possib
	
	dBNode * a = hash_table_find(element_get_key(seq_to_binary_kmer("ATTTGACGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * b = hash_table_find(element_get_key(seq_to_binary_kmer("CCATCGGCA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * c = hash_table_find(element_get_key(seq_to_binary_kmer("CCGCTGCCG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * d = hash_table_find(element_get_key(seq_to_binary_kmer("GCATCACCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * e = hash_table_find(element_get_key(seq_to_binary_kmer("CCATCACCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * f = hash_table_find(element_get_key(seq_to_binary_kmer("GCATTTGAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * g = hash_table_find(element_get_key(seq_to_binary_kmer("GTATTTGAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * h = hash_table_find(element_get_key(seq_to_binary_kmer("GGATTTGAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * i = hash_table_find(element_get_key(seq_to_binary_kmer("AACAACTTA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * j = hash_table_find(element_get_key(seq_to_binary_kmer("AAGTTACCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);//TODO: validate that I and J are the o
	dBNode * k = hash_table_find(element_get_key(seq_to_binary_kmer("ACCCCATCG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	
	pathStep start_step[7];
	pathStep end_step[7];
	
	//Unconnected path
	start_step[0].node = i;
	end_step  [0].node = j;
	
	
	//b is the y node
	start_step[1].node = e;
	end_step	[1].node = a;
	
	start_step[2].node = d;
	end_step	[2].node = a;
	
	
	//c is the y node
	start_step[3].node = h;
	end_step	[3].node = e;
	start_step[4].node = g;
	end_step	[4].node = d;
	start_step[5].node = f;
	end_step	[5].node = c;
	start_step[6].node = k;
	end_step	[6].node = d;
	
	
    //Unconnected path
	start_step[0].orientation = reverse;//i
	end_step	[0].orientation = reverse;//j
	
	
	//b is the y node
	start_step[1].orientation = reverse;//e
	end_step  [1].orientation = reverse;//a
	start_step[2].orientation = reverse;//d
	end_step  [2].orientation = reverse;//a
	
	
	//c is the y node
	start_step[3].orientation = forward;//h
	end_step	[3].orientation = forward;//e
	start_step[4].orientation = forward;//g
	end_step	[4].orientation = reverse;//d
	start_step[5].orientation = forward;//f
	end_step	[5].orientation = forward;//c
	start_step[6].orientation = forward;//k
	end_step	[6].orientation = reverse;//d
	dBNode * double_y_node[2];
	
	double_y_node[0] = b;
	double_y_node[1] = c;
	
	
	
	Path * tmp_path = path_new(100, db_graph->kmer_size);
	pathStep  ps_last, ps_first;
	int z;
	db_graph_write_graphviz_file("test_read_pair_a.gv", db_graph);
	ReadPairDescriptorArray * rpda = new_read_pair_descriptor_array(2, 10, 1, 0, 3, 3, 0, false);
	/*ReadPairDescriptorArray * new_read_pair_descriptor_array(char capacity, int walk_distance, int max_coverage,
                                                         int min_bits, int start_length, int max_paths,
                                                         int min_kmers, boolean stack_as_n)*/
	
	read_pair_descriptor_array_add_pair(0, 0, 0, 18, 13,3,1,&simple_iterated_score_paths ,rpda	);
	
	db_graph_reset_flags(db_graph);
	db_graph_write_graphviz_file("test_read_pair_b.gv", db_graph);
	
	for(z = 0; z < 7; z++){
		printf("___________________________________________TEST %d_____________________________________________________________\n", z);
		path_reset(tmp_path);
		
		if(z < 3){
			db_node_action_set_flag(b, Y_START);
			db_node_action_unset_flag(c, Y_START);
		}else{
			db_node_action_set_flag(c, Y_START);
			db_node_action_unset_flag(b, Y_START);
		}
		
		read_pair_get_path(start_step[z].node,  &db_node_action_set_flag_visited, rpda, db_graph, tmp_path);
		//boolean read_pair_get_path(dBNode * node, void (*node_action) (dBNode * node),
		//				 ReadPairDescriptorArray * rpda, dBGraph * db_graph, Path * path)
		
		path_get_last_step(&ps_last, tmp_path);
		path_get_step_at_index(0, &ps_first, tmp_path);
		path_to_fasta(tmp_path, stdout);
		
		switch (z) {
			case 1:
				CU_ASSERT(path_step_equals_without_label(&ps_first, &start_step[1])==true);
				CU_ASSERT(path_step_equals_without_label(&ps_last, &end_step[1])==true);
				break;
			case 0:
				CU_ASSERT(path_step_equals_without_label(&ps_last, &end_step[0])==true);
				CU_ASSERT(path_step_equals_without_label(&ps_first, &start_step[0])==true);
				break;
			case 2:
				CU_ASSERT(path_step_equals_without_label(&ps_first, &start_step[2])==true);
				CU_ASSERT(path_step_equals_without_label(&ps_last, &end_step[2])==true);
				break;
			case 3:
				CU_ASSERT(path_step_equals_without_label(&ps_first, &start_step[3])==true);
				CU_ASSERT(path_step_equals_without_label(&ps_last, &end_step[3])==true);
				break;
			case 4:
				CU_ASSERT(path_step_equals_without_label(&ps_first, &start_step[4])==true);
				CU_ASSERT(path_step_equals_without_label(&ps_last, &end_step[4])==true);
				break;
			case 5:
				CU_ASSERT(path_step_equals_without_label(&ps_first, &start_step[5])==true);
				CU_ASSERT(path_step_equals_without_label(&ps_last, &end_step[5])==true);
				break;
			case 6:
				CU_ASSERT(path_step_equals_without_label(&ps_first, &start_step[6])==true);
				CU_ASSERT(path_step_equals_without_label(&ps_last, &end_step[6])==true);
				break;
			default:
				break;
		}
		
		
	}
	db_graph_write_graphviz_file("test_read_pair_c.gv", db_graph);
	//Now, we will remove the read pair information and try to load it back. 
	int y;
	void clean_read_pair(dBNode * node){
		for (y = 0; y <rpda->number_of_pairs; y++) {
			ReadPairDescriptor * rpd = rpda->pair[y];
			node->signature[rpd->colour][rpd->first]  = 0;
			node->signature[rpd->colour][rpd->second] = 0;
			
		}
	}
	hash_table_traverse(&clean_read_pair, db_graph);
	db_graph_write_graphviz_file("test_read_pair_clean.gv", db_graph);
	
	
	FILE * f1 = fopen("../data/test/read_pair/rp_1.fastq", "r");
	FILE * f2 = fopen("../data/test/read_pair/rp_2.fastq", "r");
	
	rpda->pair[0]->max_read_length = 100;
	read_pair_enrich_graph(f1, f2, rpda->pair[0],  db_graph);
	
	fclose(f1);
	fclose(f2);
	
	db_graph_write_graphviz_file("test_read_pair_enriched.gv", db_graph);
	//Cleaning up. 
	destroy_read_pair_descriptor_array(&rpda);
	
	path_destroy(tmp_path);
	hash_table_free(&db_graph);
}


void test_enrich_read_pair()
{
	int kmer_size = 9;
	int number_of_bits = 10; 
	int bucket_size = 30;
	long long bad_reads = 0; 
	int seq_length;
	dBGraph * db_graph;
	BinaryKmer tmp_kmer1, tmp_kmer2;
	
	path_array_destroy_buffers();
	path_array_initialise_buffers(kmer_size);
	//
	
	db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
	
	seq_length = load_fastq_from_filename_into_graph("../data/test/read_pair/rp_1.fastq",0, 0,&bad_reads, 5, 200, db_graph);
	seq_length = load_fastq_from_filename_into_graph("../data/test/read_pair/rp_2.fastq",0, 1,&bad_reads, 5, 200, db_graph);
	
	
	
	db_graph_write_graphviz_file("test_read_pair.viz", db_graph);
	//So far, we had just been reading the graph. From here, we are naming the end nodes of the different walks and possib
	
	dBNode * a = hash_table_find(element_get_key(seq_to_binary_kmer("ATTTGACGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * b = hash_table_find(element_get_key(seq_to_binary_kmer("CCATCGGCA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * c = hash_table_find(element_get_key(seq_to_binary_kmer("CCGCTGCCG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * d = hash_table_find(element_get_key(seq_to_binary_kmer("GCATCACCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * e = hash_table_find(element_get_key(seq_to_binary_kmer("CCATCACCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * f = hash_table_find(element_get_key(seq_to_binary_kmer("GCATTTGAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * g = hash_table_find(element_get_key(seq_to_binary_kmer("GTATTTGAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * h = hash_table_find(element_get_key(seq_to_binary_kmer("GGATTTGAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * i = hash_table_find(element_get_key(seq_to_binary_kmer("AACAACTTA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * j = hash_table_find(element_get_key(seq_to_binary_kmer("AAGTTACCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * k = hash_table_find(element_get_key(seq_to_binary_kmer("ACCCCATCG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	
	pathStep start_step[7];
	pathStep end_step[7];
	
	//Unconnected path
	start_step[0].node = i;
	end_step	[0].node = j;
	
	
	//b is the y node
	start_step[1].node = e;
	end_step	[1].node = a;
	start_step[2].node = d;
	end_step	[2].node = a;
	
	
	//c is the y node
	start_step[3].node = h;
	end_step	[3].node = e;
	start_step[4].node = g;
	end_step	[4].node = d;
	start_step[5].node = f;
	end_step	[5].node = c;
	start_step[6].node = k;
	end_step	[6].node = d;
	
	
    //Unconnected path
	start_step[0].orientation = reverse;//i
	end_step	[0].orientation = reverse;//j
	
	
	//b is the y node
	start_step[1].orientation = reverse;//e
	end_step	[1].orientation = reverse;//a
	start_step[2].orientation = reverse;//d
	end_step	[2].orientation = reverse;//a
	
	
	//c is the y node
	start_step[3].orientation = forward;//h
	end_step	[3].orientation = forward;//e
	start_step[4].orientation = forward;//g
	end_step	[4].orientation = reverse;//d
	start_step[5].orientation = forward;//f
	end_step	[5].orientation = forward;//c
	start_step[6].orientation = forward;//k
	end_step	[6].orientation = reverse;//d
	dBNode * double_y_node[2];
	
	double_y_node[0] = b;
	double_y_node[1] = c;
	
	
	
	Path * tmp_path = path_new(100, db_graph->kmer_size);
	pathStep  ps_last, ps_first;
	int z;
	
	ReadPairDescriptorArray * rpda = new_read_pair_descriptor_array(2, 10, 1, 2, 3, 3, 2, false);
	/*ReadPairDescriptorArray * new_read_pair_descriptor_array(char capacity, int walk_distance, int max_coverage,
                                                         int min_bits, int start_length, int max_paths,
                                                         int min_kmers, boolean stack_as_n)*/
	
	read_pair_descriptor_array_add_pair(0, 0, 0, 5, 5,3,6,&simple_iterated_score_paths ,rpda	);
	db_graph_reset_flags(db_graph);
	//y_walk_dont_mark();
	for(z = 0; z < 7; z++){
		
		path_reset(tmp_path);
		
		if(z < 3){
			db_node_action_set_flag(b, Y_START);
			db_node_action_unset_flag(c, Y_START);
		}else{
			db_node_action_set_flag(c, Y_START);
			db_node_action_unset_flag(b, Y_START);
		}
		
		read_pair_get_path(start_step[z].node, &db_node_action_set_flag_visited, rpda, db_graph, tmp_path);
		path_get_last_step(&ps_last, tmp_path);
		path_get_step_at_index(0, &ps_first, tmp_path);
		path_to_fasta(tmp_path, stdout);
		
		CU_ASSERT(path_step_equals_without_label(&ps_last, &start_step[z])==true);
		CU_ASSERT(path_step_equals_without_label(&ps_last, &end_step[z])==true);
		
	}
	
	destroy_read_pair_descriptor_array(&rpda);
	
	
	
	path_destroy(tmp_path);
	hash_table_free(&db_graph);
}




#endif

#ifdef ENABLE_READ_PAIR
void test_read_pair()
{
	int kmer_size = 9;
	int number_of_bits = 10; 
	int bucket_size = 30;
	long long bad_reads = 0; 
	int seq_length;
	dBGraph * db_graph;
	BinaryKmer tmp_kmer1, tmp_kmer2;
	
	path_array_destroy_buffers();
	path_array_initialise_buffers(kmer_size);
	//
	
	db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
//	load_fastq_from_filename_into_graph(char *filename, short colour,
//                                        long long *bad_reads,
//                                        char quality_cut_off,
//                                        int max_read_length, int fastq_ascii_offset, dBGraph * db_graph)
	seq_length = load_fastq_from_filename_into_graph("../data/test/read_pair/rp_1.fastq", 0,&bad_reads, 5, 200, 64,db_graph);
	seq_length = load_fastq_from_filename_into_graph("../data/test/read_pair/rp_2.fastq", 0,&bad_reads, 5, 200, 64,db_graph);
	
	
	
	db_graph_write_graphviz_file("test_read_pair.gv", db_graph);
	//So far, we had just been reading the graph. From here, we are naming the end nodes of the different walks and possib
	
	dBNode * a = hash_table_find(element_get_key(seq_to_binary_kmer("ATTTGACGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * b = hash_table_find(element_get_key(seq_to_binary_kmer("CCATCGGCA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * c = hash_table_find(element_get_key(seq_to_binary_kmer("CCGCTGCCG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * d = hash_table_find(element_get_key(seq_to_binary_kmer("GCATCACCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * e = hash_table_find(element_get_key(seq_to_binary_kmer("CCATCACCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * f = hash_table_find(element_get_key(seq_to_binary_kmer("GCATTTGAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * g = hash_table_find(element_get_key(seq_to_binary_kmer("GTATTTGAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * h = hash_table_find(element_get_key(seq_to_binary_kmer("GGATTTGAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * i = hash_table_find(element_get_key(seq_to_binary_kmer("AACAACTTA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * j = hash_table_find(element_get_key(seq_to_binary_kmer("AAGTTACCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);//TODO: validate that I and J are the o
	dBNode * k = hash_table_find(element_get_key(seq_to_binary_kmer("ACCCCATCG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	
	pathStep start_step[7];
	pathStep end_step[7];
	
	//Unconnected path
	start_step[0].node = i;
	end_step  [0].node = j;
	
	
	//b is the y node
	start_step[1].node = e;
	end_step	[1].node = a;
	
	start_step[2].node = d;
	end_step	[2].node = a;
	
	
	//c is the y node
	start_step[3].node = h;
	end_step	[3].node = e;
	start_step[4].node = g;
	end_step	[4].node = d;
	start_step[5].node = f;
	end_step	[5].node = c;
	start_step[6].node = k;
	end_step	[6].node = d;
	
	
    //Unconnected path
	start_step[0].orientation = reverse;//i
	end_step	[0].orientation = reverse;//j
	
	
	//b is the y node
	start_step[1].orientation = reverse;//e
	end_step  [1].orientation = reverse;//a
	start_step[2].orientation = reverse;//d
	end_step  [2].orientation = reverse;//a
	
	
	//c is the y node
	start_step[3].orientation = forward;//h
	end_step	[3].orientation = forward;//e
	start_step[4].orientation = forward;//g
	end_step	[4].orientation = reverse;//d
	start_step[5].orientation = forward;//f
	end_step	[5].orientation = forward;//c
	start_step[6].orientation = forward;//k
	end_step	[6].orientation = reverse;//d
	dBNode * double_y_node[2];
	
	double_y_node[0] = b;
	double_y_node[1] = c;
	
	
	
	Path * tmp_path = path_new(100, db_graph->kmer_size);
	pathStep  ps_last, ps_first;
	int z;
	db_graph_write_graphviz_file("test_read_pair_a.gv", db_graph);
	ReadPairDescriptorArray * rpda = new_read_pair_descriptor_array(2, 10, 1, 0, 3, 3, 0, false);
	/*ReadPairDescriptorArray * new_read_pair_descriptor_array(char capacity, int walk_distance, int max_coverage,
     int min_bits, int start_length, int max_paths,
     int min_kmers, boolean stack_as_n)*/
	
	read_pair_descriptor_array_add_pair(0, 0, 0, 18, 13,3,1,&simple_iterated_score_paths ,rpda	);
	FILE * f1 = fopen("../data/test/read_pair/rp_1.fastq", "r");
	FILE * f2 = fopen("../data/test/read_pair/rp_2.fastq", "r");
	
	rpda->pair[0]->max_read_length = 100;
	read_pair_search_enrich_graph(f1, f2, rpda->pair[0],  64, db_graph);
    
	db_graph_reset_flags(db_graph);
	db_graph_write_graphviz_file("test_read_pair_b.gv", db_graph);
	/*
	for(z = 0; z < 7; z++){
		printf("___________________________________________TEST %d_____________________________________________________________\n", z);
		path_reset(tmp_path);
		
		if(z < 3){
			db_node_action_set_flag(b, Y_START);
			db_node_action_unset_flag(c, Y_START);
		}else{
			db_node_action_set_flag(c, Y_START);
			db_node_action_unset_flag(b, Y_START);
		}
		
		read_pair_get_path(start_step[z].node,  &db_node_action_set_flag_visited, rpda, db_graph, tmp_path);
		//boolean read_pair_get_path(dBNode * node, void (*node_action) (dBNode * node),
		//				 ReadPairDescriptorArray * rpda, dBGraph * db_graph, Path * path)
		
		path_get_last_step(&ps_last, tmp_path);
		path_get_step_at_index(0, &ps_first, tmp_path);
		path_to_fasta(tmp_path, stdout);
		
		switch (z) {
			case 1:
				CU_ASSERT(path_step_equals_without_label(&ps_first, &start_step[1])==true);
				CU_ASSERT(path_step_equals_without_label(&ps_last, &end_step[1])==true);
				break;
			case 0:
				CU_ASSERT(path_step_equals_without_label(&ps_last, &end_step[0])==true);
				CU_ASSERT(path_step_equals_without_label(&ps_first, &start_step[0])==true);
				break;
			case 2:
				CU_ASSERT(path_step_equals_without_label(&ps_first, &start_step[2])==true);
				CU_ASSERT(path_step_equals_without_label(&ps_last, &end_step[2])==true);
				break;
			case 3:
				CU_ASSERT(path_step_equals_without_label(&ps_first, &start_step[3])==true);
				CU_ASSERT(path_step_equals_without_label(&ps_last, &end_step[3])==true);
				break;
			case 4:
				CU_ASSERT(path_step_equals_without_label(&ps_first, &start_step[4])==true);
				CU_ASSERT(path_step_equals_without_label(&ps_last, &end_step[4])==true);
				break;
			case 5:
				CU_ASSERT(path_step_equals_without_label(&ps_first, &start_step[5])==true);
				CU_ASSERT(path_step_equals_without_label(&ps_last, &end_step[5])==true);
				break;
			case 6:
				CU_ASSERT(path_step_equals_without_label(&ps_first, &start_step[6])==true);
				CU_ASSERT(path_step_equals_without_label(&ps_last, &end_step[6])==true);
				break;
			default:
				break;
		}
		
		
	}
	db_graph_write_graphviz_file("test_read_pair_c.gv", db_graph);
	//Now, we will remove the read pair information and try to load it back. 
	int y;
	void clean_read_pair(dBNode * node){
		for (y = 0; y <rpda->number_of_pairs; y++) {
			ReadPairDescriptor * rpd = rpda->pair[y];
			node->signature[rpd->colour][rpd->first]  = 0;
			node->signature[rpd->colour][rpd->second] = 0;
			
		}
	}
	hash_table_traverse(&clean_read_pair, db_graph);
	db_graph_write_graphviz_file("test_read_pair_clean.gv", db_graph);
	
	
	FILE * f1 = fopen("../data/test/read_pair/rp_1.fastq", "r");
	FILE * f2 = fopen("../data/test/read_pair/rp_2.fastq", "r");
	
	rpda->pair[0]->max_read_length = 100;
	read_pair_enrich_graph(f1, f2, rpda->pair[0],  db_graph);
	
	fclose(f1);
	fclose(f2);
	
	db_graph_write_graphviz_file("test_read_pair_enriched.gv", db_graph);
	//Cleaning up. 
	destroy_read_pair_descriptor_array(&rpda);
	*/
	path_destroy(tmp_path);
	hash_table_free(&db_graph);
}


void test_enrich_read_pair()
{
	int kmer_size = 9;
	int number_of_bits = 10; 
	int bucket_size = 30;
	long long bad_reads = 0; 
	int seq_length;
	dBGraph * db_graph;
	BinaryKmer tmp_kmer1, tmp_kmer2;
	
	path_array_destroy_buffers();
	path_array_initialise_buffers(kmer_size);
	//
	
	db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
	
	seq_length = load_fastq_from_filename_into_graph("../data/test/read_pair/rp_1.fastq",0, 0,&bad_reads, 5, 200, db_graph);
	seq_length = load_fastq_from_filename_into_graph("../data/test/read_pair/rp_2.fastq",0, 1,&bad_reads, 5, 200, db_graph);
	
	
	
	db_graph_write_graphviz_file("test_read_pair.viz", db_graph);
	//So far, we had just been reading the graph. From here, we are naming the end nodes of the different walks and possib
	
	dBNode * a = hash_table_find(element_get_key(seq_to_binary_kmer("ATTTGACGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * b = hash_table_find(element_get_key(seq_to_binary_kmer("CCATCGGCA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * c = hash_table_find(element_get_key(seq_to_binary_kmer("CCGCTGCCG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * d = hash_table_find(element_get_key(seq_to_binary_kmer("GCATCACCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * e = hash_table_find(element_get_key(seq_to_binary_kmer("CCATCACCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * f = hash_table_find(element_get_key(seq_to_binary_kmer("GCATTTGAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * g = hash_table_find(element_get_key(seq_to_binary_kmer("GTATTTGAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * h = hash_table_find(element_get_key(seq_to_binary_kmer("GGATTTGAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * i = hash_table_find(element_get_key(seq_to_binary_kmer("AACAACTTA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * j = hash_table_find(element_get_key(seq_to_binary_kmer("AAGTTACCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	dBNode * k = hash_table_find(element_get_key(seq_to_binary_kmer("ACCCCATCG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
	
	pathStep start_step[7];
	pathStep end_step[7];
	
	//Unconnected path
	start_step[0].node = i;
	end_step	[0].node = j;
	
	
	//b is the y node
	start_step[1].node = e;
	end_step	[1].node = a;
	start_step[2].node = d;
	end_step	[2].node = a;
	
	
	//c is the y node
	start_step[3].node = h;
	end_step	[3].node = e;
	start_step[4].node = g;
	end_step	[4].node = d;
	start_step[5].node = f;
	end_step	[5].node = c;
	start_step[6].node = k;
	end_step	[6].node = d;
	
	
    //Unconnected path
	start_step[0].orientation = reverse;//i
	end_step	[0].orientation = reverse;//j
	
	
	//b is the y node
	start_step[1].orientation = reverse;//e
	end_step	[1].orientation = reverse;//a
	start_step[2].orientation = reverse;//d
	end_step	[2].orientation = reverse;//a
	
	
	//c is the y node
	start_step[3].orientation = forward;//h
	end_step	[3].orientation = forward;//e
	start_step[4].orientation = forward;//g
	end_step	[4].orientation = reverse;//d
	start_step[5].orientation = forward;//f
	end_step	[5].orientation = forward;//c
	start_step[6].orientation = forward;//k
	end_step	[6].orientation = reverse;//d
	dBNode * double_y_node[2];
	
	double_y_node[0] = b;
	double_y_node[1] = c;
	
	
	
	Path * tmp_path = path_new(100, db_graph->kmer_size);
	pathStep  ps_last, ps_first;
	int z;
	
	ReadPairDescriptorArray * rpda = new_read_pair_descriptor_array(2, 10, 1, 2, 3, 3, 2, false);
	/*ReadPairDescriptorArray * new_read_pair_descriptor_array(char capacity, int walk_distance, int max_coverage,
     int min_bits, int start_length, int max_paths,
     int min_kmers, boolean stack_as_n)*/
	
	read_pair_descriptor_array_add_pair(0, 0, 0, 5, 5,3,6,&simple_iterated_score_paths ,rpda	);
	db_graph_reset_flags(db_graph);
	//y_walk_dont_mark();
	for(z = 0; z < 7; z++){
		
		path_reset(tmp_path);
		
		if(z < 3){
			db_node_action_set_flag(b, Y_START);
			db_node_action_unset_flag(c, Y_START);
		}else{
			db_node_action_set_flag(c, Y_START);
			db_node_action_unset_flag(b, Y_START);
		}
		
		read_pair_get_path(start_step[z].node, &db_node_action_set_flag_visited, rpda, db_graph, tmp_path);
		path_get_last_step(&ps_last, tmp_path);
		path_get_step_at_index(0, &ps_first, tmp_path);
		path_to_fasta(tmp_path, stdout);
		
		CU_ASSERT(path_step_equals_without_label(&ps_last, &start_step[z])==true);
		CU_ASSERT(path_step_equals_without_label(&ps_last, &end_step[z])==true);
		
	}
	
	destroy_read_pair_descriptor_array(&rpda);
	
	
	
	path_destroy(tmp_path);
	hash_table_free(&db_graph);
}
#endif
