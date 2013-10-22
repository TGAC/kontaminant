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
 
 
 #ifndef PERFECT_PATH_H_
#define PERFECT_PATH_H_



/**
 * Given a dBGraph, print to a file all the   perfect paths 
 */
void perfect_path_print_paths(char * filename, int max_length, int singleton_length,
		boolean with_coverages, dBGraph * db_graph) ;

/**
 * Given a dBGraph, print the supernode containing the given dBNode node
 * Mainly to keep backwards compatibility with the paradignms of db_graph.h
 * and be able to run the former code using the same logic. The node_action
 * is wrapped to a new function that will be executed alongside the step_action
 * of the involved walk. 
 */ 

int perfect_path_get_path(dBNode * node, Orientation orientation,
		void(*node_action)(dBNode * node), dBGraph * db_graph, Path * path);
		
int pefcet_path_get_path_with_callback(dBNode * node, Orientation orientation,
		void(*node_action)(dBNode * node), void(*path_action)(Path * path), dBGraph * db_graph);

WalkingFunctions *perfect_path_get_funtions(WalkingFunctions *
											walking_functions);


/**
 * This version of perfect path, gets the perfect path from a given step. 
 */
int perfect_path_get_path_from_step_with_callback(pathStep  first,
												  void (*node_action) (dBNode * node),
												  void (*path_action) (Path * path),
												  dBGraph * db_graph);

/**
 * Function that applies a function to a path.  that contains a given node.
 */												  
int perfect_path_get_path_with_callback(dBNode * node, Orientation orientation,
					void (*node_action) (dBNode * node),
					void (*path_action) (Path * path),
					dBGraph * db_graph);
					
/**
 *
 * Function that gets all the paths from a given node. It doesn't mark the nodes
 * as visited as it may be used to explore a region without modifing it. It returns
 * the number of paths that were found. 
 * 
 */ 
int perfect_path_get_all_paths_from(dBNode * node, Orientation orientation,
			    PathArray * pa, int limit, dBGraph * db_graph);


/**
 * Identity function that returns the first_step as starting step. 
 */ 
pathStep *get_first_step_identity(pathStep * first_step, dBGraph * db_graph);

int perfect_path_get_path_with_callback_with_args(dBNode * node, Orientation orientation, void (*node_action) (dBNode * , void *), void * node_args,  void (*path_action) (Path * , void *), void * path_args,dBGraph * db_graph);


/**
 * 
 */
void perfect_path_base_callback(Path * p);

/**
 * Function to find a perfect path from a node and apply an action to a node, a path passing a arguments
 * to the action funcrions
 */

int perfect_path_get_path_from_step_with_callback_with_args(pathStep  first,
                                                            void (*node_action) (dBNode * node, void *),
                                                            void * node_args,
                                                            void (*path_action) (Path * path, void *), 
                                                            void * path_args,
                                                            dBGraph * db_graph);
#endif
