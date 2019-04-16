#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef LINKED_LIST
#define LINKED_LIST
	//==============================================================================
	/*
	* Header file to contain general functions, variables, structs
	* These can be used in different linked list operations
	* Avoids redefinitions elsewhere
	*/
	//==============================================================================
	// INCLUDES
#include "common_functions.h"
// STRUCTURES
// Enum States
	enum CellState
	{
		FROZEN, NARROWBAND, UNKNOWN
	};
	//==============================================================================
	// Structure Objects
	typedef struct
	{
		enum cellState state;
		dataType arrival;
		size_t position;
		size_t xpos;
		size_t ypos;
		size_t zpos;
	} Obj_Structure;
	//==============================================================================
	// Linked List
	struct Node
	{
		enum CellState state;
		dataType arrival;
		size_t position;
		size_t xpos;
		size_t ypos;
		size_t zpos;
		struct Node* next;
	};
	//==============================================================================
	// PROTOTYPES
	// 1.a Searches linearly for an element in the linked list and updates the element value
	void searchUpdate(struct Node ** head, struct Node * objAdd, const size_t position);
	// 1.b Binary Search Function for a Linked list
	void binarySearchLinkedList(struct Node **head_ref, dataType newValue, size_t position);
	struct Node *middleNode(struct Node *startNode, struct Node *endNode);
	//==============================================================================
	// 2. Heap Sort Function for a Linked List
	struct MaxHeap
	{
		int size;
		size_t ** metadata;
		dataType* array;
	};
	void swap(dataType* a, dataType* b);
	void swapMeta(size_t ** a, size_t ** b);
	void maxHeapify(struct MaxHeap* maxHeap, int idx);
	struct MaxHeap* createAndBuildHeap(struct Node *head);
	void heap_sort(struct Node **head);
	//==============================================================================
	/* Linked List Specific FUnctions*/
	// 3. Gets the Nth element from the list
	struct Node * getElement(struct Node* head, int index);
	//==============================================================================
	// 4.a Pushes/Inserts an object at the beginning of the list
	void push(struct Node** head_ref, Obj_Structure objAdd);
	// 4.b Pushes/Inserts an element at the bottom or after another element in the list and maintains the ascending ordering as it adds
	void pushNodeBottom(struct Node** head_ref, struct Node * objsBand);
	// 4.c Pushes a Node at the beginning of the list
	void pushNode(struct Node** head_ref, struct Node * objsBand);
	//==============================================================================
	// 5 Pops/Removes the first element from a linked list
	int pop(struct Node ** head);
	//==============================================================================
	// 6 Prints the elements of the list
	void printList(struct Node *node);
	//==============================================================================
	// 7 Returns the size of the linked list - number of elements in the list
	int listSize(struct Node * head);
	//==============================================================================
	// 8 Deletes the list - Free memory used by a list
	void deleteList(struct Node * head);
	//==============================================================================
#endif // !LINKED_LIST

#ifdef __cplusplus
}
#endif
