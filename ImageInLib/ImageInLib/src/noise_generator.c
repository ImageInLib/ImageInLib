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
double upperBValue(double v);
/*
* Uniform number generator between a - b
* lower value is a
* upper value is b
*/
double randomUniformNumber(double lowerValue, double upperValue);

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

bool additive3dNoise_D(dataType ** array3DPtr, const size_t xDim, const size_t yDim, const size_t zDim, int C, dataType fgMin, dataType bgMax)
{
	size_t i, k;
	const size_t dim2D = xDim * yDim;
	dataType j;

	//checks if the memory was allocated
	if (array3DPtr == NULL)
		return false;

	srand((unsigned)time(NULL)); //seed for randon number generator
	const dataType range_half = (dataType)(ceil((0.5 * (C - 1))));

	for (k = 0; k < zDim; k++)
	{
		for (i = 0; i < dim2D; i++)
		{
			j = (dataType)(rand() % C); // 0 - C-1
			j = j - range_half;
			array3DPtr[k][i] = array3DPtr[k][i] + j;
			if (array3DPtr[k][i] < fgMin)
			{
				array3DPtr[k][i] = fgMin;
			}//end if
			else if (array3DPtr[k][i] > bgMax)
			{
				array3DPtr[k][i] = bgMax;
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

bool saltAndPepper3dNoise_D(dataType** array3DPtr, const size_t xDim, const size_t yDim,
	const size_t zDim, dataType density, const dataType pepper)
{
	size_t i, j, k, l = 0, m, n, xd, s;

	//checks for correctness of the density (ie: percent noise, on [0,1] of salt & pepper noise
	if ((density < 0) && (density > 1))
		return false;

	//checks if the memory was allocated
	if (array3DPtr == NULL)
		return false;

	//Number of voxels affected by noise
	const size_t dim3DK = (size_t)((xDim * yDim * zDim * density) + 0.5);
	const size_t dim2D = xDim * yDim;

	// Generate Random points
	srand(time(NULL)); //seed for randon number generator
	Random3dPoints* generated_points = malloc(sizeof(Random3dPoints) * dim3DK);
	bool loop = true;
	Random3dPoints* tmpRdPts;
	do
	{
		//Coordinates of the affected voxel
		i = (rand() % (xDim));
		j = (rand() % (yDim));
		k = (rand() % (zDim));

		//intensity of the affected voxel
		n = (rand() % 2) * pepper;

		xd = x_new(i, j, xDim);
		tmpRdPts = &generated_points[l];
		tmpRdPts->xd = xd;
		tmpRdPts->k = k;
		tmpRdPts->p = n;
		//========================
		l = l + 1;
		if (l == dim3DK)
		{
			loop = false;
		}
	} while ((loop) && (l <= dim3DK));

	// Addition of salt and pepper noise to 3D image
	for (s = 0; s < dim3DK; s++) {
		size_t dx = generated_points[s].xd;
		size_t dk = generated_points[s].k;
		dataType res = (dataType)generated_points[s].p;
		array3DPtr[dk][dx] = res;
	}
	return true;
}

bool saltAndPepper2dNoise_D(dataType* array2DPtr, const size_t xDim, const size_t yDim, dataType density, const dataType pepper)
{
	size_t i, j, l = 0, n, xd, s;

	//checks to make sure the density is in an allowable range(ie: percent of noise, on [0,1])
	if ((density < 0) && (density > 1))
		return false;

	//checks if the memory was allocated
	if (array2DPtr == NULL)
		return false;

	const size_t dim2DK = (size_t)((xDim * yDim * density) + 0.5);
	const size_t dim2D = xDim * yDim;

	// Generate Random points
	srand(time(NULL)); //seed for randon number generator
	Random2dPoints* generated_points = malloc(sizeof(Random2dPoints) * dim2DK);
	bool loop = true;
	Random2dPoints* tmpRdPts;

	do {
		//Coordinates of the affected voxel
		i = (rand() % (xDim));
		j = (rand() % (yDim));

		//intensity of the affected voxel
		n = (rand() % 2) * pepper;

		xd = x_new(i, j, xDim);
		tmpRdPts = &generated_points[l];
		tmpRdPts->xd = xd;
		tmpRdPts->p = n;
		//========================
		l = l + 1;
		if (l == dim2DK)
		{
			loop = false;
		}

	} while ((loop) && (l < dim2DK));

	// Addition of salt and pepper noise to 2D image
	for (s = 0; s < dim2DK; s++) {
		size_t dx = generated_points[s].xd;
		dataType res = (dataType)generated_points[s].p;
		array2DPtr[dx] = res;
	}
	return true;
}


/*
* Multiplicative noise adds noise to imageDataPtr
* 2D
*/
void addMultiplicativeNoise(dataType ** imageDataPtr, size_t imageHeight, size_t imageLength, size_t imageWidth, float variance, dataType fgMin, dataType bgMax)
{
	size_t i, j;
	// a, b values
	double upper = 0.5 * upperBValue(variance);
	double lower = -1.0 * upper;

	size_t dim2D = imageLength * imageWidth;

	for (i = 0; i < imageHeight; i++)
	{
		for (j = 0; j < dim2D; j++)
		{
			imageDataPtr[i][j] = imageDataPtr[i][j] + (dataType)(randomUniformNumber(lower, upper)*imageDataPtr[i][j]);
			imageDataPtr[i][j] = (imageDataPtr[i][j] < fgMin) ? fgMin : imageDataPtr[i][j];
			imageDataPtr[i][j] = (imageDataPtr[i][j] > bgMax) ? bgMax : imageDataPtr[i][j];
		}
	}
}
/*
* Structural noise generate adds noise to imageDataPtr
* 2D
*/
void addStructuralNoise(dataType ** imageDataPtr, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType fgMin, dataType bgMax)
{
	size_t i, j, k, xd;
	size_t dim2D = imageLength * imageWidth;
	// Create a mesh grid for the height and width
	dataType ** periodicNoise = (dataType **)malloc(imageHeight * sizeof(dataType*));

	for (k = 0; k < imageHeight; k++)
	{
		periodicNoise[k] = (dataType *)malloc(dim2D * sizeof(dataType));
		for (i = 0; i < imageLength; i++)
		{
			for (j = 0; j < imageWidth; j++)
			{
				xd = x_new(i, j, imageLength);
				// This formula is modifiable, any sinusoid equation combination can be used!
				periodicNoise[k][xd] = (dataType)(100 * cos(2 * M_PI * ((12. * j) / 512. + (16. * i) / 512.)) + 100.);
				imageDataPtr[k][xd] = (imageDataPtr[k][xd] + (dataType)((periodicNoise[k][xd] / 2.)) / 2);

				imageDataPtr[k][xd] = (imageDataPtr[k][xd] < fgMin) ? fgMin : imageDataPtr[k][xd];
				imageDataPtr[k][xd] = (imageDataPtr[k][xd] > bgMax) ? bgMax : imageDataPtr[k][xd];
			}
		}
	}
}

double upperBValue(double v)
{
	return (double)sqrt(v);
}

double randomUniformNumber(double lowerValue, double upperValue)
{
	if (lowerValue == upperValue)
	{
		return lowerValue;
	}
	else
	{
		double range = 2.0 * (rand() / (RAND_MAX / (upperValue - lowerValue)));
		return lowerValue + (range * upperValue);
	}
}