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
 
 void db_graph_identify_branches(int max_length, dBGraph * db_graph); 

void db_graph_walk_branches(char * filename, int total_max_length, int bubble_max_length, int bubble_max_depth, dBGraph * db_graph);

boolean db_graph_compare_paths(PathArray * patharray, pathStep * end, int kmer_size);

int db_graph_found_matched_paths(PathArray * patharray,
					int max_length, Orientation orientation,
					dBNode * start_node, pathStep * end_step,
					char * filename, dBGraph * db_graph);

Path * db_graph_get_surrounding_path(dBNode * node,	Orientation orientation,
					boolean reverse_path, int max_length, dBGraph * db_graph);

void db_graph_walk_display_paths(PathArray * patharray);

void db_graph_prepare_output_files(char * filename);

int db_graph_output_search_paths(char * filename, PathArray * patharray);

char * db_graph_generate_colour_stats_string(Path * p, Path * prefix, Path * suffix, char * type);
