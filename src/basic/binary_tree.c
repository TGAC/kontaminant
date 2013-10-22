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

#include <binary_tree.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

static TreeElement * tree_element_buffer_get_slot(TreeElementBuffer * teb);

#if 0
static TreeElement *bt_find_min(TreeElement * te)
{
	
	TreeElement *next = te->left;
	while (next != NULL) {
		next = next->left;
	}
	return next;
}

static TreeElement *bt_find_max(TreeElement * te)
{
	
	TreeElement *next = te->right;
	while (next != NULL) {
		next = next->right;
	}
	return next;
}

static TreeElement *bt_find_sustitute(TreeElement * te)
{
	
	return NULL;
}
#endif

static TreeElement *binary_tree_get_root(BinaryTree * bt)
{
	return bt->root;
}

BinaryTree *binary_tree_new(int initial_capacity,
							int (*compare_element) (void *a, void *b))
{
	BinaryTree *tmp = calloc(1, sizeof(BinaryTree));
	tmp->capacity = 0;
	tmp->size = 0;
	tmp->compare_element = compare_element;
	//tmp->elements = calloc(initial_capacity, sizeof(TreeElement));
	tmp->elements = tree_element_buffer_new(initial_capacity);
	return tmp;
	
}

void binary_tree_destroy(BinaryTree ** bt)
{
	free((*bt)->elements);
	free(*bt);
	bt = NULL;
	
}

void binary_tree_clean(BinaryTree * bt)
{
	int i,j;
	TreeElement * te;
	for(j=0; j<bt->elements->number_of_blocks;j++){
		for (i = 0; i < bt->elements->block_size; i++) {	//Pesimistic approach, with the view of setting everything to 0
			te = bt->elements->tree_element_arays[j];
			te[i].left = NULL;
			te[i].right = NULL;
			te[i].head = NULL;
			te[i].element = NULL;
			
		}
	}
	
}

static TreeElement *binary_tree_get_next_free_slot(BinaryTree * bt)
{
	TreeElement* tmp_te = tree_element_buffer_get_slot(bt->elements);
	
	bt->size++;
	return tmp_te;
	
}

static void binary_tree_undo_last_element(BinaryTree * bt){
	TreeElementBuffer * teb = bt->elements;
	int min_slot = teb->first_block_with_free_elements;
	teb->free_elements[min_slot]++;
	bt->size--;
	assert(teb->free_elements[min_slot] <= teb->block_size);
}

int binary_tree_add_leaf(TreeElement * adding, BinaryTree * bt)
{
	int compared = 0;
	TreeElement *current = binary_tree_get_root(bt);
	if (current == NULL) {
		bt->root = adding;
		return 1;
	}
	TreeElement *next = current;
	TreeSide ts = BT_CENTRE;
	int tries = bt->size;
	do {
		current = next;
		compared =
		bt->compare_element(adding->element, current->element);
		assert(compared != 0);
		if (compared > 0) {
			next = current->right;
			//fprintf(stderr, "[binary_tree_add_leaf]%pRIGHT....\n", next);
			ts = BT_RIGHT;
		} else if (compared < 0) {
			//      fprintf(stderr, "[binary_tree_add_leaf]%pLEFT....\n", next);
			next = current->left;
			ts = BT_LEFT;
		} else {
			fprintf(stderr,
					"[binary_tree_add_leaf]The element is already in the list....\n");
			return 0;
		}
		if (tries-- == 0) {
			fprintf(stderr, "[binary_tree_add_leaf]To much!...\n");
			exit(-1);
		}
	} while (next != NULL);
	
	adding->head = current;	//Keep track of the parent to simplify other algorithms. 
	switch (ts) {
		case BT_LEFT:
			current->left = adding;
			break;
		case BT_RIGHT:
			current->right = adding;
			break;
		default:
			fprintf(stderr,
					"[binary_tree_add_leaf]Something went wrong with the switch!\n");
			exit(-1);
			break;
	}
	return 1;
	
}

void binary_tree_add_element(void *element, BinaryTree * bt)
{
	TreeElement *te = binary_tree_get_next_free_slot(bt);
	te->element = element;
	int added =  binary_tree_add_leaf(te, bt);
	if(added == 0){
		binary_tree_undo_last_element(bt);
		te->element= NULL;
		te->left = NULL;
		te->right = NULL;
		te->head = NULL;
		
	}
	
}

/**
 * binary_tree_get_size
 * Returns the size of the binary three
 * It just calls the size in the data structure.
 */ 
int binary_tree_get_size(BinaryTree * bt){
	assert(bt!=NULL);
	return bt->size;
}

void binary_tree_increase_capacity(BinaryTree * bt){
	fprintf(stderr, "[binary_tree_increase_capacity] not yet implemented");
	exit(-1);
}

void binary_tree_remove_element(void *element, BinaryTree * bt)
{
	//TODO! the logic to remove the node, without destroying the content...
	fprintf(stderr, "[binary_tree_remove_element] not yet implemented");
	exit(-1);
}

TreeElement *binary_tree_find(void *element, BinaryTree * bt)
{
	assert(bt!=NULL);
	assert(element !=NULL);
	assert(bt->compare_element != NULL);
	int compared = 0;
	TreeElement *current = binary_tree_get_root(bt);	//TODO: be able to rellocate the head to anywhere. 
	
	TreeSide ts = BT_CENTRE;
	TreeElement *found = NULL;
	while (current != NULL && found == NULL) {
		assert(current->element !=NULL);
		compared = bt->compare_element(element, current->element);
		if (compared > 0) {
			current = current->right;
			ts = BT_RIGHT;
		} else if (compared < 0) {
			current = current->left;
			ts = BT_LEFT;
		} else {
			found = current;
		}
		
	}
	return found;
}

void *binary_tree_find_or_add(void *element,  void * (*copy_function)(void * el) ,BinaryTree * bt)
{
	assert(bt!=NULL);
	assert(element !=NULL);
	assert(bt->compare_element != NULL);
	int compared = 0;
	TreeElement *current = binary_tree_get_root(bt);	//TODO: be able to rellocate the head to anywhere. 
	
	TreeSide ts = BT_CENTRE;
	TreeElement *found = NULL;
	while (current != NULL && found == NULL) {
		assert(current->element !=NULL);
		compared = bt->compare_element(element, current->element);
		if (compared > 0) {
			current = current->right;
			ts = BT_RIGHT;
		} else if (compared < 0) {
			current = current->left;
			ts = BT_LEFT;
		} else {
			found = current;
		}
	}
    
    
    
    if (found == NULL) {
        void * element2 = copy_function(element);
        binary_tree_add_element(element2, bt);
        return element2;
    }else{
        return found->element;
    }
}

static void walk(TreeElement * te, void (*f) (void *e))
{
	if (te->left != NULL) {
		walk(te->left, f);
	}
	f(te->element);
	if (te->right != NULL) {
		walk(te->right, f);
	}
	
}

static void reverse_walk(TreeElement * te, void (*f) (void *e))
{
	if (te->right != NULL) {
		walk(te->right, f);
	}
	f(te->element);
	if (te->left != NULL) {
		walk(te->left, f);
	}
	
}

void binary_tree_sorted_walk(void (*f) (void *e), BinaryTree * bt)
{
	assert(bt != NULL);
	assert(f != NULL);
	TreeElement * te = binary_tree_get_root(bt);
	if (te == NULL) {
		return;
	}
	walk(te, f);
}

void binary_tree_reverse_sorted_walk(void (*f) (void *e), BinaryTree * bt)
{
	assert(bt != NULL);
	assert(f != NULL);
	TreeElement * te = binary_tree_get_root(bt);
	if (te == NULL) {
		return;
	}
	reverse_walk(te, f);
}

TreeElementBuffer * tree_element_buffer_new(int block_size){
	TreeElementBuffer * teb = calloc(1, sizeof(TreeElementBuffer));
	teb->tree_element_arays = calloc(_TREE_ELEMNT_BUFFER_DEFAULT_SIZE, sizeof(TreeElement *));
	teb->free_elements = calloc(_TREE_ELEMNT_BUFFER_DEFAULT_SIZE, sizeof(int));
	teb->block_size = block_size;
	teb->number_of_blocks = 0;
	teb->tree_element_arays_size = _TREE_ELEMNT_BUFFER_DEFAULT_SIZE;
	teb->first_block_with_free_elements = -1;
	return teb;
}

//Checks if the array is empty. It doesn´t do any lock checking, so beware if used
//in a multithreaded environment. To make it work in a multithreaded environment, the calling
//function should be atomic. 
#if 0
static TreeElement * tree_element_array_get_empty(TreeElement * te, int size){
	int i;
	
	TreeElement * temp = NULL;
	for (i = 0; i < size && temp == NULL; i++) {
		if (te[i].element == NULL) {
			temp = &te[i];
		}
	}
	return te;	
}
#endif

static TreeElement * tree_element_array_alloc(int size){
	TreeElement * te = calloc(size, sizeof(TreeElement));
	if (te==NULL) {
		fprintf(stderr, "[tree_element_array_alloc] not enough memory to allocate more tree elements\n");
		exit(-1);
	}
	return te;
}

static int find_next_free_slot(int from, TreeElementBuffer * teb){
	int i, ret = -1;
	if(from >= teb->number_of_blocks){
		return teb->number_of_blocks;
	}
	for(i = from; i < teb->number_of_blocks && ret == -1; i++){
		if(teb->free_elements[i] > 0){
			ret = i;
		}
	}
	
	return ret != -1? ret:i;//If we couldnt find, tag that is the next block, to create a new one...
	
}

static void binary_tree_array_grow(TreeElementBuffer *teb){
	//int new_size = teb->tree_element_arays_size * 2;
	int new_size = (int)(ceil((double)teb->tree_element_arays_size * 1.1));
	int old_size=teb->tree_element_arays_size ;
	int * free_elements = realloc(teb->free_elements, new_size * sizeof(int) );
	int i;
	if(free_elements == NULL){
		fprintf(stderr, "[binary_tree_array_grow] Unable to  allocate more memory for the binary tree (free_elements).\n");
		exit(-1);
	}
	
	
	
	//teb->free_elements = free_elements;
	
	
	TreeElement ** tree_element_arays = realloc(teb->tree_element_arays, new_size * sizeof(TreeElement *) );
	
	if(tree_element_arays == NULL){
		fprintf(stderr, "[binary_tree_array_grow] Unable to  allocate more memory for the binary tree (tree_element_aray).\n");
		exit(-1);
	}
	
	for (i=old_size; i<new_size; i++) {
		free_elements[i] = 0;
		tree_element_arays[i] = NULL;
	}
	
	teb->free_elements = free_elements;
	teb->tree_element_arays = tree_element_arays;
	teb->tree_element_arays_size = new_size;
}

static void binary_tree_alloc(TreeElementBuffer *teb){
	
	if(teb->number_of_blocks == teb->tree_element_arays_size){
		//fprintf(stderr, "Array capacity (before): %d\n", teb->tree_element_arays_size);
		binary_tree_array_grow(teb);
		//	fprintf(stderr, "Array capacity(after): %d\n", teb->tree_element_arays_size);
	}
	
	TreeElement * new_te = tree_element_array_alloc(teb->block_size);
	
	teb->tree_element_arays[teb->number_of_blocks] = new_te;
	teb->free_elements[teb->number_of_blocks] = teb->block_size;
	teb->number_of_blocks++;
	// fprintf(stderr, "Number of blocks: %d\n", teb->number_of_blocks);
}

static TreeElement * tree_element_array_find_free( int size, int free_elements, TreeElement * tea){
	int i = size - free_elements;
	int j = i; //To set the limit if we can't find a free element from the number of free elements. Avoid iterating the array twice. 
	TreeElement * te = NULL;
	TreeElement * tmp = NULL;
	for(; i < size && te == NULL; i++){ //We assume that the array has been allocated in order, hence we can start searching at the end of it. If we can't find anything free, we search again. 
		tmp = &tea[i];
		if(tmp->element == NULL){
			te = tmp;
		}
	}
	
	for(i = 0; i<j && te == NULL ;i++){
		tmp = &tea[i];
		if(tmp->element == NULL){
			te = tmp;
		}
	}
	
	assert(te!=NULL);
	return te;
}



static TreeElement * tree_element_buffer_get_slot(TreeElementBuffer * teb){
	TreeElement * te = NULL;
	TreeElement * tea = NULL;
	int min_slot = teb->first_block_with_free_elements;
	if(min_slot == -1 || teb->free_elements[min_slot] == 0){
		min_slot = find_next_free_slot(min_slot, teb);
		if(min_slot ==  teb->number_of_blocks){
			//fprintf(stderr, "ALLOCATING NEW TREE BLOCK\n");
			binary_tree_alloc(teb);
			
		}
		teb->first_block_with_free_elements = min_slot;
	}
	
	tea = teb->tree_element_arays[min_slot];
	te = tree_element_array_find_free(teb->block_size, teb->free_elements[min_slot], tea);
	teb->free_elements[min_slot]--;
	//te = teb->tree_element_arays[teb->number_of_blocks-1];
	return te;
	
}

