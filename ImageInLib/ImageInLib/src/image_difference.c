/*
* Author: Markjoe Olunna UBA
* Purpose: ImageInLife project - 4D Image Segmentation Methods
* Language:  C
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "image_difference.h"
#include "common_functions.h"

bool imageDifference3dDataArrayD(dataType ** image3DPtr1, dataType ** image3DPtr2, dataType ** outputImage3DPtr,
	const size_t xDim, const size_t yDim, const size_t zDim)
{
	size_t i, k;//loop counter for z dimension
	const size_t dim2D = xDim * yDim;

	//checks if the memory was allocated
	if (image3DPtr1 == NULL || image3DPtr2 == NULL || outputImage3DPtr == NULL)
		return false;

	/* storage of the absolute value of difference of the processed image and original image*/
	for (k = 0; k < zDim; k++)
	{
		for (i = 0; i < dim2D; i++) 
		{
			outputImage3DPtr[k][i] = (dataType)fabs(image3DPtr1[k][i] - image3DPtr2[k][i]);
		}
	}

	return true;
}