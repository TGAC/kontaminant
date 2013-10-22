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
 
 
#ifndef Y_WALK_H_
#define Y_WALK_H_

void mark_double_y(dBGraph * db_graph);
void y_walk_print_paths(char *filename, int max_length,int singleton_length, 
			boolean with_coverages, boolean with_viz, dBGraph * db_graph);

Path *y_walk_get_path(dBNode * node, Orientation orientation,
		      void (*node_action) (dBNode * node),
		      dBGraph * db_graph, boolean both_directions,Path * path);

void y_walk_dont_mark();//Not to be used in production, but to control the unit testing. 
#endif 
