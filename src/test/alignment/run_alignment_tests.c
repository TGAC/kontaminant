

#include <CUnit.h>
#include <Basic.h>
#include <test_alignment.h>
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



 if (NULL == CU_add_test(pSuite, "Testing  alignment initialization", test_initialized_alignment )){
    CU_cleanup_registry();
    return CU_get_error();
  }
  
  if (NULL == CU_add_test(pSuite, "Testing perfect alignment", test_perfect_alignment )){
    CU_cleanup_registry();
    return CU_get_error();
  }
  
   if (NULL == CU_add_test(pSuite, "Testing alignment with gap second seq",  test_alignment_with_gap_1)){
    CU_cleanup_registry();
    return CU_get_error();
  }
  
   if (NULL == CU_add_test(pSuite, "Testing alignment with gap first seq",  test_alignment_with_gap_2)){
    CU_cleanup_registry();
    return CU_get_error();
  }
   
  if (NULL == CU_add_test(pSuite, "Testing alignment with gap on both seqs",  test_alignment_with_gap_3)){
    CU_cleanup_registry();
    return CU_get_error();
  }
  
  if (NULL == CU_add_test(pSuite, "Testing alignment with viral sequence",  test_alignment_virus)){
    CU_cleanup_registry();
    return CU_get_error();
  }
  
  if (NULL == CU_add_test(pSuite, "Testing homopolyper with insertion",  test_alignment_homopolymer_ins)){
    CU_cleanup_registry();
    return CU_get_error();
  }
  
  if (NULL == CU_add_test(pSuite, "Testing homopolymer with deletion",  test_alignment_homopolymer_del)){
    CU_cleanup_registry();
    return CU_get_error();
  }
 
  if (NULL == CU_add_test(pSuite, "Testing non homopolyper with insertion",  test_alignment_non_homopolymer_ins)){
    CU_cleanup_registry();
    return CU_get_error();
  }
  
  if (NULL == CU_add_test(pSuite, "Testing non homopolymer with deletion",  test_alignment_non_homopolymer_del)){
    CU_cleanup_registry();
    return CU_get_error();
  }
  
  if (NULL == CU_add_test(pSuite, "Testing homopolymer with mutation",  test_alignment_homopolymer_mut)){
    CU_cleanup_registry();
    return CU_get_error();
  }
  
  /* Run all tests using the CUnit Basic interface */
  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
 
  CU_cleanup_registry();
  return CU_get_error();


}





