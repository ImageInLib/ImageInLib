
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

bool thresholding3dFunctionN(dataType** image3DPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType thres_min, dataType thres_max, dataType backGround, dataType forGround) {
	size_t i, j, k;

	if (image3DPtr == NULL)
		return false;

	for (k = 0; k < zDim; k++) {
		for (i = 0; i < xDim; i++) {
			for (j = 0; j < yDim; j++) {
				if (image3DPtr[k][x_new(i, j, xDim)] >= thres_min && image3DPtr[k][x_new(i, j, xDim)] <= thres_max) {
					image3DPtr[k][x_new(i, j, xDim)] = forGround;
				}
				else {
					image3DPtr[k][x_new(i, j, xDim)] = backGround;
				}
			}
		}
	}
}

bool thresholdingOTSU(dataType** image3DPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType background, dataType foreground) {

	size_t i, j, k;

	//Find min and max data
	dataType min_data = 1000000, max_data = -10000000;
	for (k = 0; k < zDim; k++) {
		for (i = 0; i < xDim; i++) {
			for (j = 0; j < yDim; j++) {
				if (image3DPtr[k][x_new(i, j, xDim)] < min_data)
					min_data = image3DPtr[k][x_new(i, j, xDim)];
				if (image3DPtr[k][x_new(i, j, xDim)] > max_data)
					max_data = image3DPtr[k][x_new(i, j, xDim)];
			}
		}
	}

	//Compute histogram
	size_t totalClass = (size_t)(max_data - min_data + 1);
	size_t* histogram = (size_t*)malloc(totalClass * sizeof(size_t));
	for (i = 0; i < totalClass; i++) 
		histogram[i] = 0;
	for (k = 0; k < zDim; k++) {
		for (i = 0; i < xDim; i++) {
			for (j = 0; j < yDim; j++) {
				histogram[(short)image3DPtr[k][x_new(i, j, xDim)]]++;
			}
		}
	}

	//Number of cells
	size_t numberOfCells = 0;
	for (i = 0; i < totalClass; i++) {
		numberOfCells = numberOfCells + histogram[i];
	}

	//Compute probability
	dataType* Proba = (dataType*)malloc(totalClass * sizeof(dataType));
	for (i = 0; i < totalClass; i++) Proba[i] = 0;
	dataType sumProba = 0;
	for (i = 0; i < totalClass; i++) {
		Proba[i] = (dataType)histogram[i] / numberOfCells;
		//sumProba = sumProba + Proba[i];  // just for verification
	}

	dataType* interClassVariance = (dataType*)malloc(totalClass * sizeof(dataType));
	for (i = 0; i < totalClass; i++) interClassVariance[i] = 0;
	size_t T, optimalThresholdValue = 0;
	dataType sigma = 0, sum_weight;
	dataType weightClass_b, meanClass_b, varClass_b;
	dataType weightClass_f, meanClass_f, varClass_f;

	for (T = 0; T < totalClass; T++) {
		weightClass_b = 0, meanClass_b = 0, varClass_b = 0;
		weightClass_f = 0, meanClass_f = 0, varClass_f = 0;
		//compute class weight
		for (i = 0; i <= T; i++) {
			weightClass_b = weightClass_b + Proba[i];
		}
		for (i = T + 1; i < totalClass; i++) {
			weightClass_f = weightClass_f + Proba[i];
		}
		sum_weight = weightClass_b + weightClass_f;
		//compute class mean
		for (i = 0; i <= T; i++) {
			meanClass_b = meanClass_b + i * Proba[i];
		}
		meanClass_b = meanClass_b / weightClass_b;
		for (i = T + 1; i < totalClass; i++) {
			meanClass_f = meanClass_f + i * Proba[i];
		}
		meanClass_f = meanClass_f / weightClass_f;
		//compute class variance
		for (i = 0; i <= T; i++) {
			varClass_b = varClass_b + (i - meanClass_b) * (i - meanClass_b) * Proba[i];
		}
		varClass_b = varClass_b / weightClass_b;
		for (i = T + 1; i < totalClass; i++) {
			varClass_f = varClass_f + (i - meanClass_f) * (i - meanClass_f) * Proba[i];
		}
		varClass_f = varClass_f / weightClass_f;

		//compute inter-class variance
		sigma = weightClass_b * varClass_b + weightClass_f * varClass_f;
		//sigma = sqrt(sigma);
		interClassVariance[T] = sigma;
	}

	dataType minVariance = 1e10;
	for (T = 0; T < totalClass; T++) {
		if (interClassVariance[T] < minVariance) {
			minVariance = interClassVariance[T];
		}
	}
	for (T = 0; T < totalClass; T++) {
		if (interClassVariance[T] == minVariance) {
			optimalThresholdValue = T;
		}
	}

	//Threshold
	for (k = 0; k < zDim; k++) {
		for (i = 0; i < xDim; i++) {
			for (j = 0; j < yDim; j++) {
				if (image3DPtr[k][x_new(i, j, xDim)] < optimalThresholdValue) {
					image3DPtr[k][x_new(i, j, xDim)] = background;
				}
				else {
					image3DPtr[k][x_new(i, j, xDim)] = foreground;
				}
			}
		}
	}

	free(histogram); free(interClassVariance);

	return true;
}