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

dataType timespacel2norm3dDataArrayD(dataType ** dataArray3DPtr1, dataType ** dataArray3DPtr2, const size_t xDim, const size_t yDim, const size_t zDim, dataType h)
{
	size_t i, k;//loop counter for z dimension
	const size_t dim2D = xDim * yDim;
	dataType productOfSpaceSumAndtime = 0, sumPower = 0;
	dataType tau = h * h;
	dataType hhh = h * h * h;

	//checks if the memory was allocated
	if (dataArray3DPtr1 == NULL || dataArray3DPtr2 == NULL)
		return false;

	/*Computation of norm*/
	for (k = 0; k < zDim; k++)
	{
		for (i = 0; i < dim2D; i++)
		{
			sumPower += (dataType)(pow(dataArray3DPtr1[k][i] - dataArray3DPtr2[k][i], 2) * hhh);
		}
	}
	productOfSpaceSumAndtime = sumPower * tau;// sumPower * (step * tau) or sumPower * tau

	return productOfSpaceSumAndtime;
}