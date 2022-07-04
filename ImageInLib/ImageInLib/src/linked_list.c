#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "linked_list.h"

// 1.a Searches linearly for an element in the linked list and updates the element value
void searchUpdate(struct Node ** head, struct Node * objAdd, const size_t position)
{
	struct Node *current = (*head);

	while (current != NULL)
	{
		if (current->position == position)
		{
			current = objAdd;
			return;
		}
		current = current->next;
	}
	// Delete list after
	//deleteList(current); // Important in memory optimisation
}
// 1.b Binary Search Function for a Linked list
void binarySearchLinkedList(struct Node **head_ref, dataType newValue, size_t position)
{
	struct Node *startNode = (*head_ref), *endNode = NULL;
	do
	{
		struct Node *middle = middleNode(startNode, endNode);
		if (middle == NULL)
		{
			// Searched Element Not Present
			break;
		}
		if (middle->position == position)
		{
			// Update
			middle->arrival = newValue;
			break;
		}
		else if (middle->position < position)
		{
			// need to search in upper half
			startNode = middle->next;
		}
		else
		{
			// need to search in lower half
			endNode = middle;
		}
	} while (endNode == NULL || endNode->next != startNode);
}
struct Node *middleNode(struct Node *startNode, struct Node *endNode)
{
	if (startNode == NULL)
	{
		// If the linked list is empty
		return NULL;
	}
	struct Node *slowPtr = startNode, *fastPtr = startNode;
	while (fastPtr != endNode)
	{
		fastPtr = fastPtr->next;

		if (fastPtr != endNode)
		{
			slowPtr = slowPtr->next;
			fastPtr = fastPtr->next;

			/*
			Note that for each loop iteration,
			slowPtr moves just one location
			while fastPtr moves two nodes at a time.
			*/
		}
	}
	return slowPtr;
	/*
	At the end , the slowPtr will be
	pointing to the middle node
	*/
}
// 2. Heap Sort Function for a Linked List
void swap(dataType* a, dataType* b)
{
	dataType t = *a;
	*a = *b;
	*b = t;
}
void swapMeta(size_t ** a, size_t ** b)
{
	size_t * t = *a;
	*a = *b;
	*b = t;
}
void maxHeapify(struct MaxHeap* maxHeap, int idx)
{
	int largest = idx;  // Initialize largest as root
	int left = (idx << 1) + 1;  // left = 2*idx + 1
	int right = (idx + 1) << 1; // right = 2*idx + 2

	// See if left child of root exists and is greater than root
	if (left < maxHeap->size && maxHeap->array[left] > maxHeap->array[largest])
		largest = left;

	// See if right child of root exists and is greater than the largest so far
	if (right < maxHeap->size &&
		maxHeap->array[right] > maxHeap->array[largest])
		largest = right;

	// Change root, if needed
	if (largest != idx)
	{
		swap(&maxHeap->array[largest], &maxHeap->array[idx]);
		swapMeta(&maxHeap->metadata[largest], &maxHeap->metadata[idx]);
		// Also swap the corresponding metadata
		maxHeapify(maxHeap, largest);
	}
}
// Creates heap from sent list
struct MaxHeap* createAndBuildHeap(struct Node *head)
{
	struct Node *current; // Set as the current Node
	current = head;
	int size = listSize(current);

	int i;
	struct MaxHeap* maxHeap = (struct MaxHeap*) malloc(sizeof(struct MaxHeap));
	// Temp array for arrivals
	dataType* tmpPtr = malloc(sizeof(dataType)*size);
	size_t **meta = malloc(sizeof(dataType)*size);
	for (i = 0; i < size; i++)
	{
		meta[i] = malloc(sizeof(dataType) * 5);
	}
	int count = 0;
	while (current != NULL)
	{
		tmpPtr[count] = current->arrival;
		size_t position[5] = { current->xpos,current->ypos,current->zpos,1,current->position };
		for (i = 0; i < 5; i++)
		{
			meta[count][i] = position[i];
		}
		current = current->next;
		count++;
	}
	maxHeap->size = count; // initialize size of heap
	maxHeap->array = tmpPtr; // Assign address of first element of array
	maxHeap->metadata = meta; // Assign address of first element of metadata
	// Start from bottommost and rightmost internal mode and heapify all internal modes in bottom up way
	for (i = (maxHeap->size - 2) / 2; i >= 0; --i)
		maxHeapify(maxHeap, i);
	return maxHeap;
}
void heap_sort(struct Node **head)
{
	// Build a heap from the input data.
	struct MaxHeap* maxHeap = createAndBuildHeap((*head));
	// Repeat following steps while heap size is greater than 1.
	// The last element in max heap will be the minimum element
	// New band
	struct Node * newBand = NULL; // Holds all the Objects
	struct Node* new_node = (struct Node*) malloc(sizeof(struct Node));
	int counts = maxHeap->size - 1, count = counts;
	size_t i;
	while (maxHeap->size > 1)
	{
		// The largest item in Heap is stored at the root. Replace
		// it with the last item of the heap followed by reducing the
		// size of heap by 1.
		swap(&maxHeap->array[0], &maxHeap->array[maxHeap->size - 1]);
		// Swap meta data
		swapMeta(&maxHeap->metadata[0], &maxHeap->metadata[maxHeap->size - 1]);
		--maxHeap->size;  // Reduce heap size

		// Finally, heapify the root of tree.
		maxHeapify(maxHeap, 0);
	}
	// Update the new band
	deleteList((*head)); // Memory Control - if not uses large memry resources!
	while (counts >= 0)
	{
		// Create a new list - pushing at the top
		new_node->arrival = maxHeap->array[counts];
		new_node->xpos = maxHeap->metadata[counts][0];
		new_node->ypos = maxHeap->metadata[counts][1];
		new_node->zpos = maxHeap->metadata[counts][2];
		new_node->state = NARROWBAND;
		new_node->position = maxHeap->metadata[counts][4];
		// Push to new band
		pushNode(&newBand, new_node);

		counts--;
	}
	// Update the band
	*head = newBand;
	new_node = NULL;
	// Free Memory After
	free(maxHeap->array);
	for (i = 0; i < count + 1; i++)
	{
		free(maxHeap->metadata[i]);
	}
	free(maxHeap->metadata);
	maxHeap = NULL;
}
// 3. Gets the Nth element from the list
struct Node * getElement(struct Node* head, int index)
{
	struct Node* current = head, *temp = (struct Node*) malloc(sizeof(struct Node));
	// the index of the node we're currently looking at
	int count = 0;
	while (current != NULL)
	{
		if (count == index)
		{
			temp->arrival = current->arrival;
			temp->position = current->position;
			temp->state = current->state;
			temp->xpos = current->xpos;
			temp->ypos = current->ypos;
			temp->zpos = current->zpos;

			temp->next = NULL;
			return(temp);
		}

		count++;
		current = current->next;
	}

	/* if we get to this line, the caller was asking for a non-existent element*/
	// Therefore return uninitialized Object

	return NULL;
	/*deleteList(current); // <-- Memory watchers!
	free(temp);*/  // <-- Memory watchers!
}
// 4.a Pushes/Inserts an element at the beginning of the list
void push(struct Node** head_ref, Obj_Structure objAdd)
{
	/* allocate node */
	struct Node* new_node = (struct Node*) malloc(sizeof(struct Node));
	/* put in the data */
	new_node->arrival = objAdd.arrival;
	new_node->position = objAdd.position;
	new_node->state = objAdd.state;
	new_node->xpos = objAdd.xpos;
	new_node->ypos = objAdd.ypos;
	new_node->zpos = objAdd.zpos;
	///* link the old list off the new node */
	new_node->next = (*head_ref);
	/* move the head to point to the new node */
	(*head_ref) = new_node;
	// Delete list after being used
	//new_node = NULL; // <-- Memory watchers!
}
// 4.b Pushes/Inserts an element at the bottom or after another element in the list and maintains the ascending ordering as it adds
void pushNodeBottom(struct Node** head_ref, struct Node * objsBand)
{
	/* allocate node */
	struct Node* new_node = (struct Node*) malloc(sizeof(struct Node));
	// Create a temp
	struct Node* temp_node = (struct Node*) malloc(sizeof(struct Node));
	/* put in the data */
	new_node->arrival = objsBand->arrival;
	new_node->position = objsBand->position;
	new_node->state = objsBand->state;
	new_node->xpos = objsBand->xpos;
	new_node->ypos = objsBand->ypos;
	new_node->zpos = objsBand->zpos;
	new_node->next = NULL;
	// Check for first iteration
	if (*head_ref == NULL)
	{
		*head_ref = new_node; // Adds at the beginning
		return;
	}
	else
	{
		// Loop through and find the last
		struct Node *current = (*head_ref);
		struct Node *prev = NULL;
		while (current->next != NULL)
		{
			prev = current;
			// Check if current arrival is <- new arrival
			if (current->arrival >= new_node->arrival)
			{
				prev->next = new_node;
				return;
			}
			else if (current->next->arrival >= new_node->arrival) // Checks if the next is the last in the list
			{
				new_node->next = current->next;

				prev->next = new_node;
				return;
			}
			else // The new value is greater than any value in the list. Will be added at the end of the list
			{
				current = current->next;
			}
		}
		// At the last element
		// Check if new value < || >
		if (current->arrival <= new_node->arrival)
		{
			// Large or the same value as the current last object, add at the end
			current->next = new_node;
			return;
		}
		else
		{
			new_node->next = current;
			prev->next = new_node;
			return;
		}
	}
}
// 4.c Pushes a Node at the beginning of the list
void pushNode(struct Node** head_ref, struct Node * objsBand)
{
	/* allocate node */
	struct Node* new_node = (struct Node*) malloc(sizeof(struct Node));
	/* put in the data */
	new_node->arrival = objsBand->arrival;
	new_node->position = objsBand->position;
	new_node->state = objsBand->state;
	new_node->xpos = objsBand->xpos;
	new_node->ypos = objsBand->ypos;
	new_node->zpos = objsBand->zpos;

	new_node->next = (*head_ref);
	/* move the head to point to the new node */
	(*head_ref) = new_node;
	// Delete list after being used
	//new_node = NULL; // <-- Memory watchers!
}
// 5 Pops/Removes the first element from a linked list
int pop(struct Node ** head)
{
	struct Node * next_node = NULL;

	if (*head == NULL) {
		return -1;
	}
	next_node = (*head)->next;
	size_t retval = (*head)->position; // Removes by position
	free(*head);
	*head = next_node;
	return 0;
	//free(next_node); // <-- Memory watchers!
}
// 6 Prints the elements of the list
void printList(struct Node *node)
{
	while (node != NULL)
	{
		printf("T = %6.4lf state = %2d xPos = %2zd yPos = %2zd zpos = %2zd position = %6zd \n",
			node->arrival, node->state, node->xpos, node->ypos, node->zpos, node->position);
		node = node->next;
	}
}
// 7 Returns the size of the linked list - number of elements in the list
int listSize(struct Node * head)
{
	struct Node *current = head;
	int size = 0;
	while (current != NULL)
	{
		++size;
		current = current->next;
	}
	return size;
}
// 8 Deletes the list - Free memory used by a list
void deleteList(struct Node * head)
{
	struct Node *current, *next; // Set as the current Node
	current = head;
	next = head;

	while (current) {
		next = current->next;
		free(current);
		current = next;
	}
}