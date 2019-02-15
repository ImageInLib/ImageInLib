/*
* Author: Markjoe Olunna UBA
* Purpose: ImageInLife project - 4D Image Segmentation Methods
* Language:  C
*/
#include <stdio.h>
#include <stdlib.h>
#include "printing_function.h"

bool printing_function3D(unsigned char **array3DPtr, const size_t xDim, const size_t yDim, const size_t zDim)
{
	size_t k;

	//checks if the memory was allocated
	if (array3DPtr == NULL)
		return false;

	// Addition of salt and pepper noise to 3D image
	for (k = 0; k < zDim; k++)
	{
		printing_function2D(array3DPtr[k], xDim, yDim);
	}
	return true;
}

bool printing_function2D(unsigned char *array2DPtr, const size_t xDim, const size_t yDim)
{
	size_t i, j;
	//checks if the memory was allocated
	if (array2DPtr == NULL)
		return false;

	//accessing the array and printing of values to the console
	for (i = 0; i <= (xDim + 1); i++) {
		for (j = 0; j <= (yDim + 1); j++) {


			if ((i == 0) && (j == 0)) {
				printf("%c", 218);//top left corner
			}//end if

			else if ((i == (xDim + 1)) && (j == yDim + 1)) {
				printf("%c", 217);//bottom right corner
			}//end if

			else if ((i == 0) && (j == (yDim + 1))) {
				printf("%c", 191);//top right corner
			}//end if

			else if ((i == (xDim + 1)) && (j == 0)) {
				printf("%c", 192);//bottom left corner
			}//end if

			else if ((i >= 0 && i <= (xDim + 1)) && (j == 0)) {
				printf("%c", 179);// boundary for first column
			}//end if

			else if ((i >= 0 && i <= (xDim + 1)) && (j == (yDim + 1))) {
				printf("%c", 179); // boundary for last column
			}//end if

			else if ((i == 0) && (j >= 0 && j <= (yDim + 1))) {
				printf("%c", 196); // boundary for first row
			}//end if

			else if ((i == (xDim + 1)) && (j >= 0 && j <= (yDim + 1))) {
				printf("%c", 196);// boundary for last row
			}//end if
			else if ((*(array2DPtr + ((i - 1) * yDim) + (j - 1)) == 255)
				|| (*(array2DPtr + ((i - 1) * yDim) + (j - 1)) == 0)) {
				printf("%c", 219); //printing a shaded rectangle
			}//end if
			else {
				printf("%d", *(array2DPtr + ((i - 1) * yDim) + (j - 1)));//printing a blank space
			}
		}
		printf("\n");
	}
	return true;
}
