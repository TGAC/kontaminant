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
#include <time.h>
#include <stdint.h>
#include <global.h>
#include <flags.h>

#ifdef THREADS
//#include <pthread.h>
/*#include <semaphore.h>

#static boolean init = true;
#static sem_t flag_sem;
#static void init_sem(){
#	init = false;
#	sem_init(&flag_sem, 0, 1);
#}

#static void lock_flag(){
	
#	static struct timespec ts;
#	ts.tv_nsec = 50000;//TODO make this dinamic, random wait?
#	int locked;
#	int count = 0;
#	do{
#		locked = sem_timedwait(&flag_sem, &ts);
#		if(locked != 0){
#			ts.tv_nsec = rand()%500000;
#			count++;
#			if(count > 1000){
#				fprintf(stderr, "[flags.lock_flag] to mcuh waiting for the lock!");
#				exit(-1);
#			}
#		}
#	}while(locked != 0);
}

static void unlock_flag(){
		sem_post(&flag_sem);
	
}
 */
#endif

boolean flags_check_for_flag(Flags f, Flags * db)
{
	//if(init)
	//	init_sem();
	//lock_flag();
	//printf("Flag %d & %d\n", *db, f);
	boolean r =(*db & f) == f;
	 
//	#unlock_flag();
	return r;
}

boolean flags_check_for_any_flag(Flags f, Flags * db)
{
	//#if(init)
	//#	init_sem();
	boolean r = (*db & f) > 0; 
	//#unlock_flag();
	return r;
}

void flags_action_clear_flags(Flags * db)
{
//#ifdef THREADS
//	if(init)
//		init_sem();
//	lock_flag();
//#endif
	*db = ALL_OFF;
//#ifdef THREADS
//	unlock_flag();
//#endif
}

void flags_action_set_flag(Flags f, Flags * db)
{
#ifdef THREADS
//	if(init)
//		init_sem();
//	lock_flag();
#endif
	*db = *db | f;
#ifdef THREADS
//	unlock_flag();
#endif
}

void flags_action_unset_flag(Flags f, Flags * db)
{
#ifdef THREADS
//	if(init)
//		init_sem();
//	lock_flag();
#endif
	*db = *db & ~f;
#ifdef THREADS
//	unlock_flag();
#endif
}

/**
 * From here, the functions are for benchmarking. Once they are mature enough, we can move them to a new file. 
 * 
 * 
 * 
 */ 

// Keep track of most recent reading of cycle counter 
static clock_t start;



void benchmark_start_counter()
{

 start = clock();
}



double benchmark_get_counter()
{

	return (double) clock() - start;

}



double benchmark_stop(){
	return  benchmark_get_counter() /*/ (MHZ * 1e6)*/;
}
