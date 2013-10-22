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
 
#ifndef _TREE_ELEMNT_H
#define _TREE_ELEMNT_H

#define _TREE_ELEMNT_BUFFER_DEFAULT_SIZE 10

typedef enum {

	BT_LEFT = -1 ,BT_RIGHT = 1, BT_CENTRE = 0
} TreeSide;

typedef struct tree{	
	struct tree * left; 
	struct tree * right;
	struct tree * head;
	void * element;
}TreeElement;

//To make the buffer multithreaded, it may be worthy to add a lock for the structure. 
typedef struct {
	TreeElement ** tree_element_arays;
	int * free_elements;
	int block_size;
	int first_block_with_free_elements; //The index of the smallest block with at least one free element. Starts in -1 when none has been allocated. 
	int number_of_blocks;
	int tree_element_arays_size;//The size of the first dimension of tree_elements_arrays
}TreeElementBuffer;

typedef struct {
	int capacity;
	int size;
	int last;
	TreeElementBuffer * elements;//The head of the tree is always the node in the possition 0;
//	TreeElement * elements;
	TreeElement * root;
	int(*compare_element)(void * a, void * b);// A function with a binary comparator  
}BinaryTree;

BinaryTree * binary_tree_new(int initial_capacity,int(*compare_element)(void * a, void * b) );

void binary_tree_destroy(BinaryTree ** bt);

void binary_tree_clean(BinaryTree * bt);

//Returns if it was able to add the leaf or not. 
int binary_tree_add_leaf(TreeElement * adding, BinaryTree * bt);

void binary_tree_add_element(void * element, BinaryTree * bt);

void binary_tree_remove_element(void * element, BinaryTree * bt);

TreeElement * binary_tree_find(void * element,  BinaryTree * tree);

void *binary_tree_find_or_add(void *element,void * (*copy_function)(void * el) , BinaryTree * bt);

void binary_tree_sorted_walk(void (*f)(void *e), BinaryTree * bt);

int binary_tree_get_size(BinaryTree * bt);

void binary_tree_reverse_sorted_walk(void (*f) (void *e), BinaryTree * bt);

TreeElementBuffer * tree_element_buffer_new(int block_size);



#endif


