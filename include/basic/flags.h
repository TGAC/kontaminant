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
 
#ifndef FLAGS_H_
#define FLAGS_H_
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

//#include <global.h>

//Global flags
#define  ALL_OFF  		   	 0

#ifndef BUILD_RAINBOWS
//TODO make a more facny set of ifs to decide to use or not the short version of Flags
#ifndef SHORT_FLAGS
#define SHORT_FLAGS
#endif
#endif
/**
 * This will be used to mark from which paths the node has been traversed
 * The list of defines give the value for each mask.
 *
 */
#ifndef SHORT_FLAGS
typedef uint32_t Flags;

//Block of flags for the element
#define  VISITED_A  		  (1 <<  0) //x0000 0001
#define  VISITED_C  		  (1 <<  1) //x0000 0002
#define  VISITED_G  		  (1 <<  2) //x0000 0004
#define  VISITED_T  		  (1 <<  3) //x0000 0008
#define  EDGES_FLAGS		  (VISITED_A 	| VISITED_C 	| VISITED_G 	| VISITED_T)     //x0000 000F
#define  VISITED_REV_A  	  (1 <<  4) //x0000 0010
#define  VISITED_REV_C  	  (1 <<  5) //x0000 0020
#define  VISITED_REV_G  	  (1 <<  6) //x0000 0040
#define  VISITED_REV_T        (1 <<  7) //x0000 0080
#define  EDGES_REV_FLAGS	  (VISITED_REV_A | VISITED_REV_C | VISITED_REV_G | VISITED_REV_T) //x0000 00F0
#define  VISITED_FORWARD  	  (1 <<  8) //x0000 0100
#define  VISITED_REVERSE  	  (1 <<  9) //x0000 0200
#define  CURRENT_PATH_FORWARD (1 << 10) //x0000 0400
#define  CURRENT_PATH_REVERSE (1 << 11) //x0000 0800
#define  END_NODE_FORWARD  	  (1 << 12) //x0000 1000
#define  END_NODE_REVERSE  	  (1 << 13) //x0000 2000
#define  BRANCH_NODE_FORWARD  (1 << 14) //x0000 4000
#define  BRANCH_NODE_REVERSE  (1 << 15) //x0000 8000
#define  PLAIN_NODE 		  (1 << 16) //x0001 0000
#define  PRINT_FORWARD 		  (1 << 17) //x0002 0000
#define  PRINT_REVERSE 		  (1 << 18) //x0004 0000
#define  STARTING_FORWARD	  (1 << 19) //x0008 0000
#define  TIP_START            (1 << 20) //x0010 0000
#define  IGNORE_START_NODE	  (1 << 21) //x0020 0000
#define  X_NODE			 	  (1 << 22) //x0040 0000
#define  VISITED			  (1 << 23) //x0080 0000
#define  PRUNED				  (1 << 24) //x0100 0000
#define  ASSIGNED			  (1 << 25) //x0200 0000
#define  Y_START			  (1 << 26) //x0400 0000

#else
typedef uint16_t Flags;
#define  ASSIGNED			  (1 << 0)  //x0001
#define  VISITED			  (1 << 1)  //x0002 
#define  PRUNED				  (1 << 2)  //x0003 
#define  BRANCH_NODE_FORWARD  (1 << 3)  //x0004 
#define  BRANCH_NODE_REVERSE  (1 << 4)  //x0010
#define  X_NODE			 	  (1 << 5)  //x0020
#define  END_NODE_FORWARD  	  (1 << 6)  //x0040
#define  END_NODE_REVERSE  	  (1 << 7)  //x0080
#define  STARTING_FORWARD	  (1 << 8)  //x0100 Used to mark where a new SOLiD read starts. 
#define  TIP_START            (1 << 9)  //x0200 To mark where the tip clip starts
#define  IGNORE_START_NODE	  (1 << 10) //x0400
#define  CURRENT_PATH_FORWARD (1 << 11) //x0800
#define  CURRENT_PATH_REVERSE (1 << 12) //x1000
#define	 Y_START			  (1 << 13) //x2000 
#define  VISITED_FORWARD  	  (1 << 14) //x4000
#define  VISITED_REVERSE  	  (1 << 15) //x8000
#endif





boolean flags_check_for_flag(Flags f, Flags * db);

boolean flags_check_for_any_flag(Flags f, Flags * db);

void flags_action_clear_flags(Flags * db);

void flags_action_set_flag(Flags f, Flags * db);

void flags_action_unset_flag(Flags f, Flags * db);

//Maybe the following should live somewhere else
void benchmark_start_counter();
double benchmark_get_counter();
double benchmark_stop();
#endif /* FLAGS_H_ */
