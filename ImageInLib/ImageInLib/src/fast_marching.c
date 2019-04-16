#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
//==============================================================================
#include "fast_marching.h"
//==============================================================================
// Local Prototypes
// Solves the Eikonal Equation 3D
dataType eikonalSolve3D(/*dataType uniform_speed, dataType grid_spacex, dataType grid_spacey*/ dataType T1, dataType T2, dataType T3);
// Solves the Quadratic Equation
dataType quadraticSolve(dataType a, dataType b, dataType c);
// Swaps two pointers values
void swapper(dataType *num1, dataType *num2);
// Sorts 3 pointers
void smallSort(dataType *a, dataType *b, dataType *c);
// Sorts 2 pointer values
void smallSort2D(dataType *a, dataType *b);
//==============================================================================
// Fast Marching Method
void fastMarching3D(struct Node * band, Obj_Structure ** object, Point3D points[], Arrival_Time *known, size_t imageHeight, size_t imageLength, size_t imageWidth, size_t countPoints)
{
	// 1. Initialization
	size_t i;
	for (i = 0; i < countPoints; i++)
	{
		const int x = (int)(points[i].x + 0.5);
		const int y = (int)(points[i].y + 0.5);
		const int z = (int)(points[i].z + 0.5);

		// 2D representation
		size_t xpos = x_new(x, y, imageLength);
		// Initialize Object structure
		object[z][xpos].arrival = known[i].T; // Either from provided arrival times or Zero
		object[z][xpos].state = FROZEN; // Known point values initialized as Frozen
		object[z][xpos].position = x_flat(x, y, z, imageLength, imageWidth); // 3D to 1D representation
		object[z][xpos].xpos = x;
		object[z][xpos].ypos = y;
		object[z][xpos].zpos = z;
		// Create a list to hold the object Neighbors
		struct Node * neighbours = NULL;
		// Update neighbors position, coordinates
		if (z - 1 > -1 && z - 1 < imageHeight && x > -1 && x < imageLength && y > -1 && y < imageWidth)
		{
			// Z top
			object[z - 1][xpos].position = x_flat(x, y, z - 1, imageLength, imageWidth);
			object[z - 1][xpos].xpos = x;
			object[z - 1][xpos].ypos = y;
			object[z - 1][xpos].zpos = z - 1;
			push(&neighbours, object[z - 1][xpos]);
		}
		if (z + 1 > -1 && z + 1 < imageHeight && x > -1 && x < imageLength && y > -1 && y < imageWidth)
		{
			// Z Bottom
			object[z + 1][xpos].position = x_flat(x, y, z + 1, imageLength, imageWidth);
			object[z + 1][xpos].xpos = x;
			object[z + 1][xpos].ypos = y;
			object[z + 1][xpos].zpos = z + 1;
			push(&neighbours, object[z + 1][xpos]);
		}
		if (z > -1 && z < imageHeight && x - 1 > -1 && x - 1 < imageLength && y > -1 && y < imageWidth)
		{
			// X left
			size_t xposL = x_new(x - 1, y, imageLength);
			object[z][xposL].position = x_flat(x - 1, y, z, imageLength, imageWidth);
			object[z][xposL].xpos = x - 1;
			object[z][xposL].ypos = y;
			object[z][xposL].zpos = z;
			push(&neighbours, object[z][xposL]);
		}
		if (z > -1 && z < imageHeight && x + 1 > -1 && x + 1 < imageLength && y > -1 && y < imageWidth)
		{
			// X Right
			size_t xposR = x_new(x + 1, y, imageLength);
			object[z][xposR].position = x_flat(x + 1, y, z, imageLength, imageWidth);
			object[z][xposR].xpos = x + 1;
			object[z][xposR].ypos = y;
			object[z][xposR].zpos = z;
			push(&neighbours, object[z][xposR]);
		}
		if (z > -1 && z < imageHeight && x > -1 && x < imageLength && y - 1 > -1 && y - 1 < imageWidth)
		{
			// Y Begin
			size_t xposB = x_new(x, y - 1, imageLength);
			object[z][xposB].position = x_flat(x, y - 1, z, imageLength, imageWidth);
			object[z][xposB].xpos = x;
			object[z][xposB].ypos = y - 1;
			object[z][xposB].zpos = z;
			push(&neighbours, object[z][xposB]);
		}
		if (z > -1 && z < imageHeight && x > -1 && x < imageLength && y + 1 > -1 && y + 1 < imageWidth)
		{
			// Y Begin
			size_t xposE = x_new(x, y + 1, imageLength);
			object[z][xposE].position = x_flat(x, y + 1, z, imageLength, imageWidth);
			object[z][xposE].xpos = x;
			object[z][xposE].ypos = y + 1;
			object[z][xposE].zpos = z;
			push(&neighbours, object[z][xposE]);
		}
		while (neighbours != NULL)
		{
			if (neighbours->state == FROZEN)
			{
				neighbours = neighbours->next;// Goes to the next neighbor in the list
				continue; // Skips calculating arrival time for this neighbor
			}
			else
			{
				// Solve T's
				dataType T01, T02, T03;
				// Access Object Coordinate Positions
				size_t posx = neighbours->xpos;
				size_t posy = neighbours->ypos;
				size_t posz = neighbours->zpos;
				size_t xy = x_new(posx, posy, imageLength);
				// Y Begin
				Obj_Structure objBegin;
				if (posz > -1 && posz < imageHeight && posy - 1 > -1 && posy - 1 < imageWidth && posx > -1 && posx < imageLength)
				{
					size_t yB = x_new(posx, posy - 1, imageLength);
					objBegin = object[posz][yB];
				}
				else
				{
					objBegin.arrival = INFINITY;
				}
				// Y End
				Obj_Structure objEnd;
				if (posz > -1 && posz < imageHeight && posy + 1 > -1 && posy + 1 < imageWidth && posx > -1 && posx < imageLength)
				{
					size_t yE = x_new(posx, posy + 1, imageLength);
					objEnd = object[posz][yE];
				}
				else
				{
					objEnd.arrival = INFINITY;
				}
				// X Left
				Obj_Structure objLeft;
				if (posz > -1 && posz < imageHeight && posy > -1 && posy < imageWidth && posx - 1 > -1 && posx - 1 < imageLength)
				{
					size_t xL = x_new(posx - 1, posy, imageLength);
					objLeft = object[posz][xL];
				}
				else
				{
					objLeft.arrival = INFINITY;
				}
				// X Right
				Obj_Structure objRight;
				if (posz > -1 && posz < imageHeight && posy > -1 && posy < imageWidth && posx + 1 > -1 && posx + 1 < imageLength)
				{
					size_t xR = x_new(posx + 1, posy, imageLength);
					objRight = object[posz][xR];
				}
				else
				{
					objRight.arrival = INFINITY;
				}
				// Z Top
				Obj_Structure objTop;
				if (posz - 1 > -1 && posz - 1 < imageHeight && posy > -1 && posy < imageWidth && posx > -1 && posx < imageLength)
				{
					objTop = object[posz - 1][xy];
				}
				else
				{
					objTop.arrival = INFINITY;
				}
				// Z Bottom
				Obj_Structure objBottom;
				if (posz + 1 > -1 && posz + 1 < imageHeight && posy > -1 && posy < imageWidth && posx > -1 && posx < imageLength)
				{
					objBottom = object[posz + 1][xy];
				}
				else
				{
					objBottom.arrival = INFINITY;
				}
				T01 = min(objBegin.arrival, objEnd.arrival);
				T02 = min(objLeft.arrival, objRight.arrival);
				T03 = min(objTop.arrival, objBottom.arrival);

				neighbours->arrival = eikonalSolve3D(T01, T02, T03);
				if (neighbours->state == NARROWBAND)
				{
					// Update Band
					searchUpdate(&band, neighbours, neighbours->position);
				}
				else
				{
					// Change state
					neighbours->state = NARROWBAND;
					// Insert to band
					pushNode(&band, neighbours);
				}
				// Update Object From the neighbor
				if (object[posz][xy].state != FROZEN)
				{
					object[posz][xy].position = neighbours->position;
					object[posz][xy].arrival = neighbours->arrival;
					object[posz][xy].state = neighbours->state;
					object[posz][xy].xpos = neighbours->xpos;
					object[posz][xy].ypos = neighbours->ypos;
					object[posz][xy].zpos = neighbours->zpos;
				}
			}
			neighbours = neighbours->next;
		}
	}
	// 2. Begin Iterations
	while (band != NULL)
	{
		// Gets the first element
		struct Node * objFirst = getElement(band, 0);
		// Change state to Frozen
		size_t objx = objFirst->xpos, objy = objFirst->ypos, objz = objFirst->zpos;
		size_t objxy = x_new(objx, objy, imageLength);
		object[objz][objxy].state = FROZEN;
		// Remove from band
		pop(&band);
		// Create a list to hold the object Neighbors
		struct Node * neighbourx = NULL;
		// Access the neighbors
		if (objz > -1 && objz < imageHeight && objy - 1 > -1 && objy - 1 < imageWidth && objx > -1 && objx < imageLength)
		{
			// Y Begin
			size_t yBg = x_new(objx, objy - 1, imageLength);
			push(&neighbourx, object[objz][yBg]);
		}
		if (objz > -1 && objz < imageHeight && objy + 1 > -1 && objy + 1 < imageWidth && objx > -1 && objx < imageLength)
		{
			// Y End
			size_t yEd = x_new(objx, objy + 1, imageLength);
			push(&neighbourx, object[objz][yEd]);
		}
		if (objz > -1 && objz < imageHeight && objy > -1 && objy < imageWidth && objx - 1 > -1 && objx - 1 < imageLength)
		{
			// X Left
			size_t xLt = x_new(objx - 1, objy, imageLength);
			push(&neighbourx, object[objz][xLt]);
		}
		if (objz > -1 && objz < imageHeight && objy > -1 && objy < imageWidth && objx + 1 > -1 && objx + 1 < imageLength)
		{
			// X Right
			size_t xRt = x_new(objx + 1, objy, imageLength);
			push(&neighbourx, object[objz][xRt]);
		}
		if (objz - 1 > -1 && objz - 1 < imageHeight && objy > -1 && objy < imageWidth && objx > -1 && objx < imageLength)
		{
			// Top
			push(&neighbourx, object[objz - 1][objxy]);
		}
		if (objz + 1 > -1 && objz + 1 < imageHeight && objy > -1 && objy < imageWidth && objx > -1 && objx < imageLength)
		{
			// Bottom
			push(&neighbourx, object[objz + 1][objxy]);
		}
		while (neighbourx != NULL)
		{
			if (neighbourx->state == FROZEN)
			{
				neighbourx = neighbourx->next;// Goes to the next neighbor in the list
				continue;
			}
			else
			{
				// Solve T
				dataType T11, T12, T13;
				// Access Object Coordinate Positions
				size_t pox = neighbourx->xpos;
				size_t poy = neighbourx->ypos;
				size_t poz = neighbourx->zpos;
				size_t pxy = x_new(pox, poy, imageLength);

				// Begin
				Obj_Structure neibourBegin;
				if (poz > -1 && poz < imageHeight && poy - 1 > -1 && poy - 1 < imageWidth && pox > -1 && pox < imageLength)
				{
					size_t poB = x_new(pox, poy - 1, imageLength);
					neibourBegin = object[poz][poB];
					if (neibourBegin.arrival != neibourBegin.arrival)
					{
						neibourBegin.arrival = INFINITY;
					}
				}
				else
				{
					neibourBegin.arrival = INFINITY;
				}
				// End
				Obj_Structure neibourEnd;
				if (poz > -1 && poz < imageHeight && poy + 1 > -1 && poy + 1 < imageWidth && pox > -1 && pox < imageLength)
				{
					size_t poE = x_new(pox, poy + 1, imageLength);
					neibourEnd = object[poz][poE];
				}
				else
				{
					neibourEnd.arrival = INFINITY;
				}
				//Left
				Obj_Structure neibourLeft;
				if (poz > -1 && poz < imageHeight && poy > -1 && poy < imageWidth && pox - 1 > -1 && pox - 1 < imageLength)
				{
					size_t poL = x_new(pox - 1, poy, imageLength);
					neibourLeft = object[poz][poL];
				}
				else
				{
					neibourLeft.arrival = INFINITY;
				}
				// Right
				Obj_Structure neibourRight;
				if (poz > -1 && poz < imageHeight && poy > -1 && poy < imageWidth && pox + 1 > -1 && pox + 1 < imageLength)
				{
					size_t poR = x_new(pox + 1, poy, imageLength);
					neibourRight = object[poz][poR];
				}
				else
				{
					neibourRight.arrival = INFINITY;
				}
				// Top
				Obj_Structure neibourTop;
				if (poz - 1 > -1 && poz - 1 < imageHeight && poy > -1 && poy < imageWidth && pox > -1 && pox < imageLength)
				{
					neibourTop = object[poz - 1][pxy];
				}
				else
				{
					neibourTop.arrival = INFINITY;
				}
				// Bottom
				Obj_Structure neibourBottom;
				if (poz + 1 > -1 && poz + 1 < imageHeight && poy > -1 && poy < imageWidth && pox > -1 && pox < imageLength)
				{
					neibourBottom = object[poz + 1][pxy];
				}
				else
				{
					neibourBottom.arrival = INFINITY;
				}

				T11 = min(neibourBegin.arrival, neibourEnd.arrival);
				T12 = min(neibourLeft.arrival, neibourRight.arrival);
				T13 = min(neibourTop.arrival, neibourBottom.arrival);

				dataType tmp = eikonalSolve3D(T11, T12, T13);

				if (neighbourx->arrival == INFINITY)
				{
					neighbourx->arrival = tmp;
				}
				else if (neighbourx->arrival < tmp)
				{
					neighbourx->arrival = tmp;
					// Update the band
					searchUpdate(&band, neighbourx, neighbourx->position);
					//searchSort(&band, neighbourx->arrival, neighbourx->position);
					//binarySearchLinkedList(&band, neighbourx->arrival, neighbourx->position);
				}
				if (neighbourx->state == UNKNOWN)
				{
					// Change state
					neighbourx->state = NARROWBAND;
					// Insert to band
					pushNodeBottom(&band, neighbourx);
				}
				// Update Object From the neighbor
				if (object[poz][pxy].state != FROZEN)
				{
					object[poz][pxy].position = neighbourx->position;
					object[poz][pxy].arrival = neighbourx->arrival;
					object[poz][pxy].state = neighbourx->state;
					object[poz][pxy].xpos = neighbourx->xpos;
					object[poz][pxy].ypos = neighbourx->ypos;
					object[poz][pxy].zpos = neighbourx->zpos;
				}
			}
			neighbourx = neighbourx->next;
		}
	}
	// 3. Finalization - Exit the function after calculating the new values
}
//==============================================================================
dataType eikonalSolve3D(dataType T1, dataType T2, dataType T3)
{
	dataType uniform_speed = 1.0;
	dataType grid_spacex = 1.0, grid_spacey = 1.0, grid_spacez = 1.0;
	dataType T, Tp, grid_T;
	dataType a, b, c;
	if (T1 == INFINITY && T2 == INFINITY)
	{
		// Use T3
		a = 1.0 / (pow(grid_spacez, 2));
		b = (-2 * T3) / (pow(grid_spacez, 2));
		c = (pow(T3, 2) / (pow(grid_spacez, 2))) - 1;
		if (pow(b, 2) < (4 * a*c))
		{
			return T3; // Solution needs to be at least as large as the inputs
			printf("Pause"); // Negative - -nan(inf) - check
		}
		else
		{
			Tp = quadraticSolve(a, b, c); // Calculated Arrival Time
			if (Tp < T3 || Tp == T3)
			{
				return T3; // Solution needs to be at least as large as the inputs
			}
			else
			{
				return Tp;
			}
		}
	}
	else if (T1 == INFINITY && T3 == INFINITY)
	{
		// Use T2
		a = 1.0 / (pow(grid_spacey, 2));
		b = (-2 * T2) / (pow(grid_spacey, 2));
		c = (pow(T2, 2) / (pow(grid_spacey, 2))) - 1;
		if (pow(b, 2) < (4 * a*c))
		{
			return T2; // Solution needs to be at least as large as the inputs
			printf("Pause"); // Negative - -nan(inf) - check
		}
		else
		{
			Tp = quadraticSolve(a, b, c); // Calculated Arrival Time
			if (Tp < T2 || Tp == T2)
			{
				return T2; // Solution needs to be at least as large as the inputs
			}
			else
			{
				return Tp;
			}
		}
	}
	else if (T2 == INFINITY && T3 == INFINITY)
	{
		// Use T1
		a = 1.0 / (pow(grid_spacex, 2));
		b = (-2 * T1) / (pow(grid_spacex, 2));
		c = (pow(T1, 2) / (pow(grid_spacex, 2))) - 1;
		if (pow(b, 2) < (4 * a*c))
		{
			return T1; // Solution needs to be at least as large as the inputs
			printf("Pause"); // Negative - -nan(inf) - check
		}
		else
		{
			Tp = quadraticSolve(a, b, c); // Calculated Arrival Time
			if (Tp < T1 || Tp == T1)
			{
				return T1; // Solution needs to be at least as large as the inputs
			}
			else
			{
				return Tp;
			}
		}
	}
	else if (T1 == INFINITY)
	{
		// Use T2 and T3
		a = (1 / (pow(grid_spacez, 2)) + 1 / (pow(grid_spacey, 2)))*uniform_speed;
		b = -2 * (T3 / (pow(grid_spacez, 2)) + T2 / (pow(grid_spacey, 2)))*uniform_speed;
		c = (pow(T3, 2) / (pow(grid_spacez, 2)) + pow(T2, 2) / (pow(grid_spacey, 2)))*uniform_speed - 1;
		// Check the Existence
		if (pow(b, 2) < (4 * a*c))
		{
			if (T2 < T3)
			{
				T = T2;
				grid_T = grid_spacey;
			}
			else
			{
				T = T3;
				grid_T = grid_spacez;
			}
			// Use Minimum to calculate new variables
			a = 1.0 / (pow(grid_T, 2));
			b = (-2 * T) / (pow(grid_T, 2));
			c = (pow(T, 2) / (pow(grid_T, 2))) - 1;
			if (pow(b, 2) < (4 * a*c))
			{
				return T; // minimal value
				printf("Pause"); // Negative - -nan(inf) - check
			}
			else
			{
				Tp = quadraticSolve(a, b, c); // Calculated Arrival Time
				if (Tp < T || Tp == T)
				{
					return T; // Solution needs to be at least as large as the inputs
				}
				else
				{
					return Tp;
				}
			}
		}
		else
		{
			Tp = quadraticSolve(a, b, c); // Calculated Arrival Time
			// Check if Tp < T2, T3
			smallSort2D(&T2, &T3);
			if (Tp < T2 || Tp < T3)
			{
				// Discard the largest and try again
				// Use T2
				a = 1.0 / (pow(grid_spacey, 2));
				b = (-2 * T2) / (pow(grid_spacey, 2));
				c = (pow(T2, 2) / (pow(grid_spacey, 2))) - 1;
				if (pow(b, 2) < (4 * a*c))
				{
					return T2;  // minimal value
					printf("Pause"); // Negative - -nan(inf)
				}
				else
				{
					Tp = quadraticSolve(a, b, c); // Calculated Arrival Time
					if (Tp < T2 || Tp == T2)
					{
						return T2; // Solution needs to be at least as large as the inputs
					}
					else
					{
						return Tp;
					}
				}
			}
			else // As large as the inputs
			{
				return Tp;
			}
		}
	}
	else if (T2 == INFINITY)
	{
		// Use T1 and T3
		a = (1 / (pow(grid_spacex, 2)) + 1 / (pow(grid_spacez, 2)))*uniform_speed;
		b = -2 * (T1 / (pow(grid_spacex, 2)) + T3 / (pow(grid_spacez, 2)))*uniform_speed;
		c = (pow(T1, 2) / (pow(grid_spacex, 2)) + pow(T3, 2) / (pow(grid_spacez, 2)))*uniform_speed - 1;
		// Check the Existence
		if (pow(b, 2) < 4 * a*c)
		{
			if (T1 < T3)
			{
				T = T1;
				grid_T = grid_spacex;
			}
			else
			{
				T = T3;
				grid_T = grid_spacez;
			}
			// Use Minimum to calculate new variables
			a = 1.0 / (pow(grid_T, 2));
			b = (-2 * T) / (pow(grid_T, 2));
			c = (pow(T, 2) / (pow(grid_T, 2))) - 1;
			if (pow(b, 2) < (4 * a*c))
			{
				return T; // minimal value
				printf("Pause"); // Negative - -nan(inf) - check
			}
			else
			{
				Tp = quadraticSolve(a, b, c); // Calculated Arrival Time
				if (Tp < T || Tp == T)
				{
					return T; // Solution needs to be at least as large as the inputs
				}
				else
				{
					return Tp;
				}
			}
		}
		else
		{
			Tp = quadraticSolve(a, b, c); // Calculated Arrival Time
										  // Check if Tp < T1, T3
			smallSort2D(&T1, &T3);
			if (Tp < T1 || Tp < T3)
			{
				// Discard the largest and try again
				// Use T1
				a = 1.0 / (pow(grid_spacex, 2));
				b = (-2 * T1) / (pow(grid_spacex, 2));
				c = (pow(T1, 2) / (pow(grid_spacex, 2))) - 1;
				if (pow(b, 2) < (4 * a*c))
				{
					return T1; // minimal value
					printf("Pause"); // Negative - -nan(inf)
				}
				else
				{
					Tp = quadraticSolve(a, b, c); // Calculated Arrival Time
					if (Tp < T1 || Tp == T1)
					{
						return T1; // Solution needs to be at least as large as the inputs
					}
					else
					{
						return Tp;
					}
				}
			}
			else // As large as the inputs
			{
				return Tp;
			}
		}
	}
	else if (T3 == INFINITY)
	{
		// Use T2 and T1
		a = (1 / (pow(grid_spacex, 2)) + 1 / (pow(grid_spacey, 2)))*uniform_speed;
		b = -2 * (T1 / (pow(grid_spacex, 2)) + T2 / (pow(grid_spacey, 2)))*uniform_speed;
		c = (pow(T1, 2) / (pow(grid_spacex, 2)) + pow(T2, 2) / (pow(grid_spacey, 2)))*uniform_speed - 1;
		// Check the Existence
		if (pow(b, 2) < 4 * a*c)
		{
			if (T1 < T2)
			{
				T = T1;
				grid_T = grid_spacex;
			}
			else
			{
				T = T2;
				grid_T = grid_spacey;
			}
			// Use Minimum to calculate new variables
			a = 1.0 / (pow(grid_T, 2));
			b = (-2 * T) / (pow(grid_T, 2));
			c = (pow(T, 2) / (pow(grid_T, 2))) - 1;
			if (pow(b, 2) < (4 * a*c))
			{
				return T;
				printf("Pause"); // Negative - -nan(inf) - check
			}
			else
			{
				Tp = quadraticSolve(a, b, c); // Calculated Arrival Time
				if (Tp < T || Tp == T)
				{
					return T; // Solution needs to be at least as large as the inputs
				}
				else
				{
					return Tp;
				}
			}
		}
		else
		{
			Tp = quadraticSolve(a, b, c); // Calculated Arrival Time
										  // Check if Tp < T2, T1
			smallSort2D(&T1, &T2);
			if (Tp < T1 || Tp < T2)
			{
				// Discard the largest and try again
				// Use T1
				a = 1.0 / (pow(grid_spacex, 2));
				b = (-2 * T1) / (pow(grid_spacex, 2));
				c = (pow(T1, 2) / (pow(grid_spacex, 2))) - 1;
				if (pow(b, 2) < (4 * a*c))
				{
					return T1;
					printf("Pause"); // Negative - -nan(inf)
				}
				else
				{
					Tp = quadraticSolve(a, b, c); // Calculated Arrival Time
					if (Tp < T1 || Tp == T1)
					{
						return T1; // Solution needs to be at least as large as the inputs
					}
					else
					{
						return Tp;
					}
				}
			}
			else // As large as the inputs
			{
				return Tp;
			}
		}
	}
	else
	{
		// Use All
		a = (1 / (pow(grid_spacex, 2)) + 1 / (pow(grid_spacey, 2)) + 1 / (pow(grid_spacez, 2)))*uniform_speed;
		b = -2 * (T1 / (pow(grid_spacex, 2)) + T2 / (pow(grid_spacey, 2)) + T3 / (pow(grid_spacez, 2)))*uniform_speed;
		c = (pow(T1, 2) / (pow(grid_spacex, 2)) + pow(T2, 2) / (pow(grid_spacey, 2)) + pow(T3, 2) / (pow(grid_spacez, 2)))*uniform_speed - 1;
		if (pow(b, 2) < 4 * a*c)
		{
			Tp = quadraticSolve(a, b, c); // Calculated Arrival Time
			if (T1 == T2 && T1 == T3 && T1 < Tp)
			{
				return T1;
			}
			else
			{
				smallSort(&T1, &T2, &T3);
				// Use only the smallest two
				// Use T2 and T1
				a = (1 / (pow(grid_spacex, 2)) + 1 / (pow(grid_spacey, 2)))*uniform_speed;
				b = -2 * (T1 / (pow(grid_spacex, 2)) + T2 / (pow(grid_spacey, 2)))*uniform_speed;
				c = (pow(T1, 2) / (pow(grid_spacex, 2)) + pow(T2, 2) / (pow(grid_spacey, 2)))*uniform_speed - 1;
				// Check the Existence
				if (pow(b, 2) < 4 * a*c)
				{
					if (T1 < T2)
					{
						T = T1;
						grid_T = grid_spacex;
					}
					else
					{
						T = T2;
						grid_T = grid_spacey;
					}
					// Use Minimum to calculate new variables
					a = 1.0 / (pow(grid_T, 2));
					b = (-2 * T) / (pow(grid_T, 2));
					c = (pow(T, 2) / (pow(grid_T, 2))) - 1;
					if (pow(b, 2) < (4 * a*c))
					{
						return T;
						printf("Pause"); // Negative - -nan(inf) - check
					}
					else
					{
						Tp = quadraticSolve(a, b, c); // Calculated Arrival Time
						if (Tp < T || Tp == T)
						{
							return T; // Solution needs to be at least as large as the inputs
						}
						else
						{
							return Tp;
						}
					}
				}
				else
				{
					Tp = quadraticSolve(a, b, c); // Calculated Arrival Time
												  // Check if Tp < T2, T1
					smallSort2D(&T1, &T2);
					if (Tp < T1 || Tp < T2)
					{
						// Discard the largest and try again
						// Use T1
						a = 1.0 / (pow(grid_spacex, 2));
						b = (-2 * T1) / (pow(grid_spacex, 2));
						c = (pow(T1, 2) / (pow(grid_spacex, 2))) - 1;
						if (pow(b, 2) < (4 * a*c))
						{
							return T1;
							printf("Pause"); // Negative - -nan(inf)
						}
						else
						{
							Tp = quadraticSolve(a, b, c); // Calculated Arrival Time
							if (Tp < T1 || Tp == T1)
							{
								return T1; // Solution needs to be at least as large as the inputs
							}
							else
							{
								return Tp;
							}
						}
					}
					else // As large as the inputs
					{
						return Tp;
					}
				}
			}
		}
		else
		{
			Tp = quadraticSolve(a, b, c); // Calculated Arrival Time
			if (Tp >= T1 && Tp >= T2 && Tp >= T3)
			{
				return Tp; // As large as the inputs
			}
			else
			{
				smallSort(&T1, &T2, &T3);
				// Use only the smallest two
				// Use T2 and T1
				a = (1 / (pow(grid_spacex, 2)) + 1 / (pow(grid_spacey, 2)))*uniform_speed;
				b = -2 * (T1 / (pow(grid_spacex, 2)) + T2 / (pow(grid_spacey, 2)))*uniform_speed;
				c = (pow(T1, 2) / (pow(grid_spacex, 2)) + pow(T2, 2) / (pow(grid_spacey, 2)))*uniform_speed - 1;
				// Check the Existence
				if (pow(b, 2) < 4 * a*c)
				{
					if (T1 < T2)
					{
						T = T1;
						grid_T = grid_spacex;
					}
					else
					{
						T = T2;
						grid_T = grid_spacey;
					}
					// Use Minimum to calculate new variables
					a = 1.0 / (pow(grid_T, 2));
					b = (-2 * T) / (pow(grid_T, 2));
					c = (pow(T, 2) / (pow(grid_T, 2))) - 1;
					if (pow(b, 2) < (4 * a*c))
					{
						return T;
						printf("Pause");
					}
					else
					{
						Tp = quadraticSolve(a, b, c); // Calculated Arrival Time
						if (Tp < T || Tp == T)
						{
							return T; // Solution needs to be at least as large as the inputs
						}
						else
						{
							return Tp;
						}
					}
				}
				else
				{
					Tp = quadraticSolve(a, b, c); // Calculated Arrival Time
												  // Check if Tp < T2, T1
					smallSort2D(&T1, &T2);
					if (Tp < T1 || Tp < T2)
					{
						// Discard the largest and try again
						// Use T1
						a = 1.0 / (pow(grid_spacex, 2));
						b = (-2 * T1) / (pow(grid_spacex, 2));
						c = (pow(T1, 2) / (pow(grid_spacex, 2))) - 1;
						if (pow(b, 2) < (4 * a*c))
						{
							return T1;
							printf("Pause"); // Negative - -nan(inf)
						}
						else
						{
							Tp = quadraticSolve(a, b, c); // Calculated Arrival Time
							if (Tp < T1 || Tp == T1)
							{
								return T1; // Solution needs to be at least as large as the inputs
							}
							else
							{
								return Tp;
							}
						}
					}
					else // As large as the inputs
					{
						return Tp;
					}
				}
			}
		}
	}
}
//==============================================================================
dataType quadraticSolve(dataType a, dataType b, dataType c)
{
	dataType resultA = (-1 * b - sqrt(pow(b, 2) - 4 * a*c)) / (2 * a);
	dataType resultB = (-1 * b + sqrt(pow(b, 2) - 4 * a*c)) / (2 * a);
	return max(resultA, resultB);
}
//==============================================================================
void swapper(dataType * num1, dataType * num2)
{
	dataType temp = *num1;

	*num1 = *num2;
	*num2 = temp;
}
//==============================================================================
void smallSort(dataType * a, dataType * b, dataType * c)
{
	if (*a > *c)swapper(a, c);
	if (*a > *b)swapper(a, b);
	if (*b > *c)swapper(b, c);
}
//==============================================================================
void smallSort2D(dataType *a, dataType *b)
{
	if (*a > *b)swapper(a, b);
}
//==============================================================================