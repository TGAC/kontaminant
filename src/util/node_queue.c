/*
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
 
/*----------------------------------------------------------------------*
 * File:    nodequeue.c                                                 *
 * Purpose: A queue for nodes!                                          *
 * Author:  Richard Leggett                                             *
 * History: 12-Aug-10: RML: Created                                     *
 *----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "binary_kmer.h"
#include "flags.h"
#include "element.h"
#include "node_queue.h"

/*----------------------------------------------------------------------*
 * General purpose queue functions                                      *
 *----------------------------------------------------------------------*/
Queue* queue_new(int n)
{
	Queue* q = calloc(1, sizeof(Queue));
    if (!q) {
        printf("Couldn't get memory for queue\n");
        exit(1);
    } else {
        //printf("Allocated!\n");
    }
	q->max_size = n;
	q->items = calloc(n, sizeof(QueueItem*));
    if (!q->items) {
        printf("Couldn't get memory for queue items\n");
        exit(1);
    }
	q->number_of_items = 0;
	return q;
}

void* queue_push(Queue* q, void* item)
{
	if (q->number_of_items == q->max_size) {
		return NULL;
	}
	
    q->items[q->number_of_items++] = item;
    
	return item;
}

void* queue_pop(Queue* q)
{
    void* item = 0;
	int i;
	
	if (q->number_of_items > 0) {
		item = q->items[0];
        
		for (i=1; i<q->number_of_items; i++) {
			q->items[i-1] = q->items[i];
		}
		
		q->number_of_items--;
	}
    
	return item;
}

void queue_free(Queue *q)
{
	if (q) {
		int i;
		for (i=0; i<q->number_of_items; i++) {
			QueueItem* qi = queue_pop(q);
            
			if (qi) {
			    free(qi);
			} 
		}
		if (q->items) {
			free(q->items);
		}
		
		free(q);
	}
}

/*----------------------------------------------------------------------*
 * Queue wrappers for nodes                                             *
 *----------------------------------------------------------------------*/
QueueItem* queue_push_node(Queue* q, dBNode* n, int d)
{
	QueueItem* item = malloc(sizeof(QueueItem));
	
	if (!item) {
		return 0;
	}
	
	item->node = n;
	item->depth = d;
    
	//BinaryKmer tmp;
	//char seq[kmer_size];	
	//binary_kmer_assignment_operator(tmp, n->kmer);
	//binary_kmer_to_seq(&tmp, kmer_size, seq);	
	//printf("Pushing %s depth %d\n", seq, d);
	
	return queue_push(q, item);
}

dBNode* queue_pop_node(Queue* q, int* d)
{
	QueueItem* item = queue_pop(q);
	dBNode* node = 0;
	
	if (item) {
		*d = item->depth;
		node = item->node;		
        
		//BinaryKmer tmp;
		//char seq[kmer_size];	
		//binary_kmer_assignment_operator(tmp, node->kmer);
		//binary_kmer_to_seq(&tmp, kmer_size, seq);
		//printf("Popped %s depth %d\n", seq, *d);
        
		free(item);
	}
	
	return node;
}
