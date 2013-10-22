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
 #include <test_dB_graph.h>
#include <test_file_reader.h>
#include <test_cyclic_count.h>
#include <test_graph_element.h>
#include <CUnit.h>
#include <Basic.h>

int  main()
{

  CU_pSuite pSuite = NULL;

  /* initialize the CUnit test registry */
  if (CUE_SUCCESS!=CU_initialize_registry())
    return CU_get_error();
  
  /* add a suite to the registry */
  pSuite = CU_add_suite("Suite_1", NULL, NULL);
  if (NULL == pSuite) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  /* add the tests to the suite */


  if (NULL == CU_add_test(pSuite, "test hash_table_find for dB_graphs",  test_hash_table_find)){
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pSuite, "test of element assignment operator (currently unused)",  test_graph_element_assign)){
    CU_cleanup_registry();
    return CU_get_error();
  }



  if (NULL == CU_add_test(pSuite, "test tip clipping",  test_tip_clipping)){
    CU_cleanup_registry();
    return CU_get_error();
  }


  if (NULL == CU_add_test(pSuite, "test node prune",  test_node_prunning_low_coverage)){
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pSuite, "test getting sliding windows where criterion for breaking a window is that a kmer is not in the graph",  test_getting_sliding_windows_where_you_break_at_kmers_not_in_db_graph )){
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pSuite, "test reading a fastq and dumping only reads that lie within the graph",  test_dumping_of_clean_fasta)){
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pSuite, "test get perfect path",  test_get_perfect_path)){
    CU_cleanup_registry();
    return CU_get_error();
  } 

   /*  
  if (NULL == CU_add_test(pSuite, "test edge coverage",  test_db_graph_db_node_edges_coverage)){
    CU_cleanup_registry();
    return CU_get_error();
  } */

  
  if (NULL == CU_add_test(pSuite, "test writting/reading binary",  test_writing_reading_binary)){
    CU_cleanup_registry();
    return CU_get_error();
  }
 
  
  if (NULL == CU_add_test(pSuite, "test get detect and smoothe bubble",  test_detect_and_smoothe_bubble)){
     CU_cleanup_registry();
     return CU_get_error();
   }
   
     if (NULL == CU_add_test(pSuite, "test get detect and remove low coverage path",  test_detect_and_remove_low_cov_path)){
     CU_cleanup_registry();
     return CU_get_error();
   }


   if (NULL == CU_add_test(pSuite, "test nodes are have coverage correctly marked on file loading",  test_coverage_is_correctly_counted_on_loading_from_file) ){
   CU_cleanup_registry();
   return CU_get_error();
  }

/*
   if (NULL == CU_add_test(pSuite, "test has precisely n edges with status",test_db_graph_db_node_has_precisely_n_edges_with_status)){
     CU_cleanup_registry();
     return CU_get_error();
   }
*/



   if (NULL == CU_add_test(pSuite, "test dumping/reading binary",  test_dump_load_binary)){
     CU_cleanup_registry();
     return CU_get_error();
   }

    
	if (NULL == CU_add_test(pSuite, "test y_walk",  test_y_walk)){
     CU_cleanup_registry();
     return CU_get_error();
   }
   
   if (NULL == CU_add_test(pSuite, "test test_y_walk_from_perfect_path_tests",  test_y_walk_from_perfect_path_tests)){
     CU_cleanup_registry();
     return CU_get_error();
   }


   
   if (NULL == CU_add_test(pSuite, "test test_get_all_paths",  test_get_all_paths)){
     CU_cleanup_registry();
     return CU_get_error();
   }
  /*
    if (NULL == CU_add_test(pSuite, "test test_get_all_paths_deep_search",  test_get_all_paths_deep_search)){
        CU_cleanup_registry();
        return CU_get_error();
    }
    */
   
   

  
/*
   if (NULL == CU_add_test(pSuite, "test calculation of N50",  test_get_N50)){
    CU_cleanup_registry();
    return CU_get_error();
  }
*/

   //if (NULL == CU_add_test(pSuite, "test function for rotating/shifting binary kmers",  test_rotate)){
   // CU_cleanup_registry();
   // return CU_get_error();
   //}
   
/*
   if (NULL == CU_add_test(pSuite, "test is condition true for all nodes in supernode",  test_is_condition_true_for_all_nodes_in_supernode)){
     CU_cleanup_registry();
     return CU_get_error();
   }
*/
  /* if (NULL == CU_add_test(pSuite, "test that can spot supernode that does not intersect any chromosome (with small but real chromosomal data)",  test_read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference))
     {
       CU_cleanup_registry();
       return CU_get_error();
     }
*/
/*
   if (NULL == CU_add_test(pSuite, "Test can pull out supernode that overlaps chromosome1 at start and end but not middle",  test_indel_discovery_simple_test_1)){
     CU_cleanup_registry();
     return CU_get_error();
   }
*/

#ifdef ENABLE_READ_PAIR

if (NULL == CU_add_test(pSuite, "Test get all paths",  test_get_all_paths)){
     CU_cleanup_registry();
     return CU_get_error();
 }

if (NULL == CU_add_test(pSuite, "Test read pair",  test_read_pair)){
     CU_cleanup_registry();
     return CU_get_error();
 }

#endif



  /* Run all tests using the CUnit Basic interface */
  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
 
  CU_cleanup_registry();
  return CU_get_error();


}





