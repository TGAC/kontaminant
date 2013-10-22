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
#include <Basic.h>
#include <inttypes.h>

#include <global.h>
#include <flags.h>
#include <binary_kmer.h>
#include <element.h>
#include <test_graph_element.h>

void test_graph_element_assign()
{
  Element e1, e2;

  BinaryKmer b1;
  binary_kmer_initialise_to_zero(&b1);
  b1[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] = (bitfield_of_64bits) 1;

  BinaryKmer b2;
  binary_kmer_initialise_to_zero(&b2);
  b2[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] = (bitfield_of_64bits) 3;

  
  element_initialise(&e1, &b1, 31);
  element_initialise(&e2, &b2, 31);

  db_node_set_edges(&e1, 0, 1);
  db_node_set_edges(&e2, 0, 2);
  
  db_node_action_set_flag(&e1, VISITED);
  db_node_action_set_flag(&e2, PRUNED);

  element_update_coverage(&e1, 0, 101);
  element_update_coverage(&e2, 0, 202);
  
  element_assign(&e1, &e2);

  CU_ASSERT(e1.coverage[0]==202);
  CU_ASSERT(e1.edges[0]==2 );
  CU_ASSERT(binary_kmer_comparison_operator(e1.kmer,b2) );
  CU_ASSERT(element_check_for_flag(&e1, PRUNED) );
}
