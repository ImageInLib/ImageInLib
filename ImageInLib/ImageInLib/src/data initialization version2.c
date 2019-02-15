/*
* Author: Markjoe Olunna UBA
* Purpose: ImageInLife project - 4D Image Segmentation Methods
* Language:  C
*/
#include <stdlib.h>
#include <stdbool.h>
#include "common_functions.h"

bool InitializeArray_2D(dataType * array2DPtr, const size_t xDim,
	const size_t yDim, const dataType value);
bool InitializeArray_3D(dataType ** array3DPtr, const size_t xDim, const size_t yDim,
	const size_t zDim, const dataType value);

bool InitializeArray_3D(dataType ** array3DPtr, const size_t xDim, const size_t yDim,
	const size_t zDim, const dataType value)
{
	size_t k;

	if (array3DPtr == NULL)
		return false;

	// filling the 2d array with number=value
	for (k = 0; k < zDim; k++)
	{
		InitializeArray_2D(array3DPtr[k], xDim, yDim, value);
	}

	return true;
}

bool InitializeArray_2D(dataType * array2DPtr, const size_t xDim, const size_t yDim, const dataType value)
{
	const size_t dim2D = xDim * yDim;
	size_t i;

	if (array2DPtr == NULL)
		return false;

	// filling the 2d array with number=value
	for (i = 0; i < dim2D; i++) {
		array2DPtr[i] = value;
	}


	return true;
}