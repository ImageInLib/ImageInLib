/*
* Author: Markjoe Olunna UBA
* Purpose: ImageInLife project - 4D Image Segmentation Methods
* Language:  C
*/
#include <stdio.h>
#include <stdlib.h>
#include "thresholding.h"
#include "common_functions.h"

bool thresholding2dFunctionD(dataType * image2DPtr, const size_t xDim, const size_t yDim,
	const dataType thresholdvalue, const dataType inputbgvalue, const dataType inputfgvalue);
bool thresholding2dFunctionUC(unsigned char * image2DPtr, const size_t xDim, const size_t yDim,
	const unsigned char thresholdvalue, const unsigned char inputbgvalue, const unsigned char inputfgvalue);


bool thresholding3dFunctionUC(unsigned char ** image3DPtr, const size_t xDim, const size_t yDim,
	const size_t zDim, const unsigned char thresholdvalue, const unsigned char inputbgvalue,
	const unsigned char inputfgvalue)
{
	size_t k;//loop counter for z dimension

			 //checks if the memory was allocated
	if (image3DPtr == NULL)
		return false;

	// filling the 2d array with number=value
	for (k = 0; k < zDim; k++)
	{
		thresholding2dFunctionUC(image3DPtr[k], xDim, yDim, thresholdvalue, inputbgvalue, inputfgvalue);
	}

	return true;
}

bool thresholding3dFunctionD(dataType ** image3DPtr, const size_t xDim, const size_t yDim,
	const size_t zDim, const dataType thresholdvalue, const dataType inputbgvalue,
	const dataType inputfgvalue)
{
	size_t k;//loop counter for z dimension

			 //checks if the memory was allocated
	if (image3DPtr == NULL)
		return false;

	// filling the 2d array with number=value
	for (k = 0; k < zDim; k++)
	{
		thresholding2dFunctionD(image3DPtr[k], xDim, yDim, thresholdvalue, inputbgvalue, inputfgvalue);
	}

	return true;
}

//function for thresholding of 2D array with some constant value.
bool thresholding2dFunctionUC(unsigned char * image2DPtr, const size_t xDim, const size_t yDim,
	const unsigned char thresholdvalue, const unsigned char inputbgvalue, const unsigned char inputfgvalue)
{
	const size_t dim2D = xDim * yDim;
	size_t i;

	//checks for successful allocation of memory
	if (image2DPtr == NULL)
		return false;

	// Thresholding of 2D slice with thresholdvalue
	for (i = 0; i < dim2D; i++) {
		if (image2DPtr[i] <= thresholdvalue)
			image2DPtr[i] = inputbgvalue;
		else
			image2DPtr[i] = inputfgvalue;
	}


	return true;
}
//function for thresholding of 2D array with some constant value.
bool thresholding2dFunctionD(dataType * image2DPtr, const size_t xDim, const size_t yDim,
	const dataType thresholdvalue, const dataType inputbgvalue, const dataType inputfgvalue)
{
	const size_t dim2D = xDim * yDim;
	size_t i;

	//checks for successful allocation of memory
	if (image2DPtr == NULL)
		return false;

	// Thresholding of 2D slice with thresholdvalue
	for (i = 0; i < dim2D; i++) {
		if (image2DPtr[i] <= thresholdvalue)
			image2DPtr[i] = inputfgvalue;
		else
			image2DPtr[i] = inputbgvalue;
	}


	return true;
}

bool thresholding3dFunctionN(dataType** image3DPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType thres_min, dataType thres_max, dataType backGround, dataType overGround) {
	size_t i, j, k;

	if (image3DPtr == NULL)
		return false;

	for (k = 0; k < zDim; k++) {
		for (i = 0; i < xDim; i++) {
			for (j = 0; j < yDim; j++) {
				if (image3DPtr[k][x_new(i, j, xDim)] >= thres_min && image3DPtr[k][x_new(i, j, xDim)] <= thres_max) {
					image3DPtr[k][x_new(i, j, xDim)] = overGround;
				}
				else {
					image3DPtr[k][x_new(i, j, xDim)] = backGround;
				}
			}
		}
	}
}