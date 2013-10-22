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
 
#ifndef CLEANING_H_
#define CLEANING_H_

typedef struct {
	 int tip_length;
     int tips_removed;
     int tips_nodes_removed;
    int marked_tips;
    dBGraph * db_graph;
} clip_tip_vars;

int cleaning_prune_low_coverage_path(int min_cov,int limit, dBGraph * db_graph);

int cleaning_remove_low_coverage(int coverage, dBGraph * db_graph);

int cleaning_remove_bubbles(int limit, dBGraph * db_graph);

int cleaning_remove_low_coverage(int coverage, dBGraph * db_graph);

void cleaning_prune_db_node(dBNode * node, dBGraph * db_graph);

int cleaning_remove_spurious_links(int max_difference, int min_coverage, dBGraph * db_graph);

int cleaning_remove_tips(int max_length, int max_it, dBGraph * db_graph);

#endif
