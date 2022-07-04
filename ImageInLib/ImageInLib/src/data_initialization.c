/*
* Author: Markjoe Olunna UBA
* Purpose: ImageInLife project - 4D Image Segmentation Methods
* Language:  C
*/
#include <stdlib.h>
#include "data_initialization.h"
#include "common_functions.h"

//function for initialization of 3D array with some constant value.
//xDim is the x dimension, yDim is the y dimension and zDim is the z dimension
//value is the initial constant value
bool initialize3dArrayUC(unsigned char ** array3DPtr, const size_t xDim, const size_t yDim,
	const size_t zDim, const unsigned char value)
{
	size_t k;//loop counter for z dimension

			 //checks if the memory was allocated
	if (array3DPtr == NULL)
		return false;

	// filling the 2d array with number=value
	for (k = 0; k < zDim; k++)
	{
		initialize2dArrayUC(array3DPtr[k], xDim, yDim, value);
	}

	return true;
}

//function for initialization of 2D array with some constant value.
bool initialize2dArrayUC(unsigned char * array2DPtr, const size_t xDim, const size_t yDim, const unsigned char value)
{
	const size_t dim2D = xDim * yDim;
	size_t i;

	if (array2DPtr == NULL)
		return false;

	// filling the 1d array with number=value
	for (i = 0; i < dim2D; i++) {
		array2DPtr[i] = value;
	}


	return true;
}

bool initialize3dArrayD(dataType ** array3DPtr, const size_t xDim, const size_t yDim,
	const size_t zDim, dataType value)
{
	size_t k;//loop counter for z dimension

			 //checks if the memory was allocated
	if (array3DPtr == NULL)
		return false;

	// filling the 2d array with number=value
	for (k = 0; k < zDim; k++)
	{
		initialize2dArrayD(array3DPtr[k], xDim, yDim, value);
	}

	return true;
}

//function for initialization of 2D array with some constant value.
bool initialize2dArrayD(dataType * array2DPtr, const size_t xDim, const size_t yDim, dataType value)
{
	const size_t dim2D = xDim * yDim;
	size_t i;

	if (array2DPtr == NULL)
		return false;

	// filling the 1d array with number=value
	for (i = 0; i < dim2D; i++) {
		array2DPtr[i] = value;
	}


	return true;
}