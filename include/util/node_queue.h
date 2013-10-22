/*----------------------------------------------------------------------*
 * File:    nodequeue.h                                                 *
 * Purpose: A queue for nodes!                                          *
 * Author:  Richard Leggett                                             *
 * History: 12-Aug-10: RML: Created                                     *
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 * Structures                                                           *
 *----------------------------------------------------------------------*/
typedef struct {
	dBNode* node;
	int depth;
} QueueItem;

typedef struct {
	int max_size;
	int number_of_items;
	int head;
	int tail;
	void** items;
} Queue;

/*----------------------------------------------------------------------*
 * General purpose queue functions                                      *
 *----------------------------------------------------------------------*/
Queue* queue_new(int n);
void* queue_push(Queue* q, void* item);
void* queue_pop(Queue* q);
void queue_free(Queue *q);

/*----------------------------------------------------------------------*
 * Queue wrappers for nodes                                             *
 *----------------------------------------------------------------------*/
QueueItem* queue_push_node(Queue* q, dBNode* n, int d);
dBNode* queue_pop_node(Queue* q, int* d);
