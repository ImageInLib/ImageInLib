/*
* Author: Markjoe Olunna UBA
* Purpose: ImageInLife project - 4D Image Segmentation Methods
* Language:  C
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "image_norm.h"
#include "common_functions.h"

double norm3dDataArrayD(dataType ** dataArray3DPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType h)
{
	size_t i, k;//loop counter for z dimension
	const size_t dim2D = xDim * yDim;
	double norm, sumPower = 0;

	//checks if the memory was allocated
	if (dataArray3DPtr == NULL)
		return false;

	/*Computation of norm*/
	for (k = 0; k < zDim; k++)
	{
		for (i = 0; i < dim2D; i++)
		{
			sumPower += (pow(dataArray3DPtr[k][i], 2) * (h * h * h));
		}
	}
	norm = sqrt(sumPower);
	return norm;
}