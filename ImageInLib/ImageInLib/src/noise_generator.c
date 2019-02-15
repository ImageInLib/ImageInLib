/*
* Authors: Polycarp Omondi Okock and Markjoe Olunna UBA
* Purpose: ImageInLife project - 4D Image Segmentation Methods
* Language:  C
*/
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "common_math.h"
#include "noise_generator.h"
#include "common_functions.h"

// Local Functions Prototypes
/*
* Calculates the upper value b in a continuous uniform function
* v is the variance
*/
dataType upperBValue(dataType v);
/*
* Uniform number generator between a - b
* lower value is a
* upper value is b
*/
dataType randomUniformNumber(dataType lowerValue, dataType upperValue);

// Function for introduction of additive noise to data(3D)
bool additive3dNoise_UC(unsigned char ** array3DPtr, const size_t xDim, const size_t yDim,
	const size_t zDim, const int C)
{
	size_t i, k;
	const size_t dim2D = xDim * yDim;
	int j;

	//checks if the memory was allocated
	if (array3DPtr == NULL)
		return false;

	srand((unsigned)time(NULL)); //seed for randon number generator

	// Addition of salt and pepper noise to 3D image
	for (k = 0; k < zDim; k++)
	{
		for (i = 0; i < dim2D; i++)
		{
			j = array3DPtr[k][i] - (C / 2) + rand() % 51;
			if (j < 0)
			{
				array3DPtr[k][i] = 0;
			}//end if
			else if (j > 255)
			{
				array3DPtr[k][i] = 255;
			}//end else
			else {
				array3DPtr[k][i] = j;
			}
		}
	}

	return true;
}

bool additive2dNoise_UC(unsigned char * array2DPtr, const size_t xDim, const size_t yDim, const int C, bool flag)
{
	const size_t dim2D = xDim * yDim;
	size_t i;
	int j;

	if (flag == true)
	srand((unsigned)time(NULL)); //seed for randon number generator

	//checks if the memory was allocated
	if (array2DPtr == NULL)
		return false;


	// Addition of salt and pepper noise to 2D image
	for (i = 0; i < dim2D; i++)
	{
		j = array2DPtr[i] - (C / 2) + rand() % 51;
		if (j < 0)
		{
			array2DPtr[i] = 0;
		}//end if
		else if (j > 255)
		{
			array2DPtr[i] = 255;
		}//end else
		else {
			array2DPtr[i] = j;
		}

	}
	return true;
}

bool additive3dNoise_D(dataType ** array3DPtr, const size_t xDim, const size_t yDim,
	const size_t zDim, dataType C)
{
	size_t i, k;
	const size_t dim2D = xDim * yDim;
	size_t j;

	//checks if the memory was allocated
	if (array3DPtr == NULL)
		return false;

	srand((unsigned)time(NULL)); //seed for randon number generator

								 // Addition of salt and pepper noise to 3D image
	for (k = 0; k < zDim; k++)
	{
		for (i = 0; i < dim2D; i++)
		{
			j = (size_t)(array3DPtr[k][i] - (C / 2) + rand() % 51 + 0.5);
			if (j < 0)
			{
				array3DPtr[k][i] = 0.;
			}//end if
			else if (j > 255)
			{
				array3DPtr[k][i] = 255.;
			}//end else
			else {
				array3DPtr[k][i] = (dataType)j;
			}

		}
	}
	return true;
}

bool additive2dNoise_D(dataType * array2DPtr, const size_t xDim, const size_t yDim, const int C, bool flag)
{
	const size_t dim2D = xDim * yDim;
	size_t i, j;

	if (flag == true)
		srand((unsigned)time(NULL)); //seed for randon number generator

									 //checks if the memory was allocated
	if (array2DPtr == NULL)
		return false;


	// Addition of salt and pepper noise to 2D image
	for (i = 0; i < dim2D; i++)
	{
		j = (size_t)(array2DPtr[i] - (C / 2) + rand() % 51 + 0.5);
		if (j < 0)
		{
			array2DPtr[i] = 0;
		}//end if
		else if (j > 255)
		{
			array2DPtr[i] = 255;
		}//end else
		else {
			array2DPtr[i] = (dataType)j;
		}
	}
	return true;
}


// Function for addition of salt and pepper noise to data
bool saltAndPepper3dNoise_UC(unsigned char ** array3DPtr, const size_t xDim, const size_t yDim,
	const size_t zDim, double K)
{
	size_t i, k, s;
	const size_t dim2DK = (int)((xDim * yDim * K) + 0.5);
	const size_t dim2D = xDim * yDim;

	//checks if the memory was allocated
	if (array3DPtr == NULL)
		return false;

	srand((unsigned)time(NULL)); //seed for randon number generator

								 //checks for correctness of the density (ie: percent noise, on [0,1] of salt & pepper noise
	if ((K < 0) && (K > 1))
		return false;

	// Addition of salt and pepper noise to 3D image
	for (k = 0; k < zDim; k++)
	{
		for (s = 0; s < dim2DK; s++) {
			i = rand() % (dim2D + 1);
			array3DPtr[k][i] = (rand() % 2) * 255;
		}
	}
	return true;
}

bool saltAndPepper2dNoise_UC(unsigned char * array2DPtr, const size_t xDim, const size_t yDim, double K, bool flag)
{
	//checks to make sure the density is in an allowable range(ie: percent noise, on [0,1])
	if ((K < 0) && (K > 1))
		return false;

	if (flag == true)
	srand((unsigned)time(NULL)); //seed for randon number generator

	const size_t dim2DK = (int)((xDim * yDim * K) + 0.5);
	const size_t dim2D = xDim * yDim;
	size_t i, k;
	srand((unsigned)time(NULL)); //seed for randon number generator

								 //checks if the memory was allocated
	if (array2DPtr == NULL)
		return false;

	// Addition of salt and pepper noise to 2D image
	for (k = 0; k < dim2DK; k++) {
		i = rand() % (dim2D + 1);
		array2DPtr[i] = (rand() % 2) * 255;
	}
	return true;
}

bool saltAndPepper3dNoise_D(dataType ** array3DPtr, const size_t xDim, const size_t yDim,
	const size_t zDim, dataType K)
{
	size_t i, k, s;
	const size_t dim2DK = (int)((xDim * yDim * K) + 0.5);
	const size_t dim2D = xDim * yDim;

	//checks if the memory was allocated
	if (array3DPtr == NULL)
		return false;

	srand((unsigned)time(NULL)); //seed for randon number generator

	//checks for correctness of the density (ie: percent of noise, on [0,1] of salt & pepper noise
	if ((K < 0) && (K > 1))
		return false;

	// Addition of salt and pepper noise to 3D image
	for (k = 0; k < zDim; k++)
	{
		for (s = 0; s < dim2DK; s++) {
			i = rand() % (dim2D + 1);
			array3DPtr[k][i] = (rand() % 2) * 255;
		}
	}

	return true;
}

bool saltAndPepper2dNoise_D(dataType * array2DPtr, const size_t xDim, const size_t yDim, dataType K, bool flag)
{
	//checks to make sure the density is in an allowable range(ie: percent of noise, on [0,1])
	if ((K < 0) && (K > 1))
		return false;

	if (flag == true)
		srand((unsigned)time(NULL)); //seed for randon number generator

	const size_t dim2DK = (int)((xDim * yDim * K) + 0.5);
	const size_t dim2D = xDim * yDim;
	size_t i, k;
	srand((unsigned)time(NULL)); //seed for randon number generator

	//checks if the memory was allocated
	if (array2DPtr == NULL)
		return false;

	// Addition of salt and pepper noise to 2D image
	for (k = 0; k < dim2DK; k++) {
		i = rand() % (dim2D + 1);
		array2DPtr[i] = (rand() % 2) * 255;
	}
	return true;
}


/*
* Multiplicative noise adds noise to imageDataPtr
* 2D
*/
void addMultiplicativeNoise(dataType ** imageDataPtr, size_t imageHeight, size_t imageWidth, dataType variance)
{
	size_t i, j;
	dataType lower = 0.0, upper = upperBValue(variance); // a, b values

	for (i = 0; i < imageHeight; i++)
	{
		for (j = 0; j < imageWidth; j++)
		{
			imageDataPtr[i][j] = imageDataPtr[i][j] + (randomUniformNumber(lower, upper) / (upper / 2.))*imageDataPtr[i][j];
		}
	}
}
/*
* Structural noise generate adds noise to imageDataPtr
* 2D
*/
void addStructuralNoise(dataType ** imageDataPtr, int imageHeight, int imageWidth)
{
	size_t i, j;

	// Create a mesh grid for the height and width
	dataType **periodicNoise = malloc(imageHeight * sizeof(dataType));

	for (i = 0; i < imageHeight; i++)
	{
		periodicNoise[i] = malloc(imageWidth * sizeof(dataType));
		for (j = 0; j < imageWidth; j++)
		{
			// This formula is modifiable, any sinusoid equation combination can be used!
			periodicNoise[i][j] = 100*cos(2*M_PI*((12. * j) / 512. + (16. * i) / 512.)) + 100.;
			imageDataPtr[i][j] = (imageDataPtr[i][j] + (periodicNoise[i][j] / 2.)) / 2.;
		}
	}
}

dataType upperBValue(dataType v)
{
	return sqrt(12 * v);
}

dataType randomUniformNumber(dataType lowerValue, dataType upperValue)
{
	if (lowerValue == upperValue)
	{
		return lowerValue;
	}
	else
	{
		return lowerValue + (rand() / (RAND_MAX / (upperValue - lowerValue)));
	}
}