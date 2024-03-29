/*
* Author: Markjoe Olunna UBA
* Purpose: ImageInLife project - 4D Image Segmentation Methods
* Language:  C
*/
#include <math.h>
#include <stdlib.h>
#include "data_initialization.h"
#include "distance_function.h"
#include "common_functions.h"

bool fastSweepingFunction_3D(dataType ** distance3DPtr, dataType ** curve3DPtr, const size_t xDim, const size_t yDim, const size_t zDim,
	const dataType h, const dataType largeValue, const dataType fgroundValue)
{
	size_t j, i, k, i_n;//i, j and k are loop counters. i_n is also a loop counter given by i_n = i + j * xDim
	size_t sweepNumber = 0, sweepDirection;
	const size_t dim2D = xDim * yDim;
	dataType ** temp3dPtr = (dataType **)malloc(sizeof(dataType*) * zDim);// temporary array used in the computation of distance
	dataType ** tmpcurve3DPtr = (dataType **)malloc(sizeof(dataType*) * zDim);// temporary array to copy curve3DPtr
	const size_t memory_alloc_size = sizeof(dataType) * dim2D;
	for (i = 0; i < zDim; i++) {
		temp3dPtr[i] = (dataType *)malloc(memory_alloc_size);
		tmpcurve3DPtr[i] = (dataType *)malloc(memory_alloc_size);
	}

	//checks if the memory was allocated
	for (i = 0; i < zDim; i++)
	{
		if (temp3dPtr[i] == NULL)
			return false;
	}

	//checks if the memory was allocated
	if (temp3dPtr == NULL)
		return false;

	if (distance3DPtr == NULL)
		return false;

	if (curve3DPtr == NULL)
		return false;

	//Checks for 2D case
	if (zDim == 1)
		return false;
	//==============================================================================
	// curve3DPtr is in the range 0-1 binary
	const double threshold = 0.5, range = 0.5 / 255.0;
	//Initialization stage
	for (k = 0; k < zDim; k++)
	{
		for (i = 0; i < xDim; i++)
		{
			for (j = 0; j < yDim; j++)
			{
				i_n = i + j * xDim;

				if (curve3DPtr[k][i_n] < threshold + range)
				{
					tmpcurve3DPtr[k][i_n] = distance3DPtr[k][i_n] = fgroundValue;
				}
				else
				{
					distance3DPtr[k][i_n] = largeValue;
					tmpcurve3DPtr[k][i_n] = curve3DPtr[k][i_n];
				}
			}
		}
	}

	//Main loop for height sweeps
	while (sweepNumber < 8)
	{
		sweepDirection = sweepNumber;
		switch (sweepDirection)
		{
		case 0:
			for (k = 0; k < zDim; k++)
			{
				for (i = 0; i < xDim; i++)
				{
					for (j = (yDim - 1); j >= 0; j--)
					{
						i_n = i + j * xDim;
						if (tmpcurve3DPtr[k][i_n] != fgroundValue)
						{
							compute3dDistance(i_n, k, xDim, yDim, zDim, distance3DPtr, temp3dPtr, dim2D, h);
						}

						if (j == 0) //because if j is unsigned type, then j(0) - 1 == type maximum value
							break;
					}
				}
			}
			sweepNumber++;
			break;

		case 1:
			for (k = 0; k < zDim; k++)
			{
				for (i = (xDim - 1); i >= 0; i--)
				{
					for (j = (yDim - 1); j >= 0; j--)
					{
						i_n = i + j * xDim;
						if (tmpcurve3DPtr[k][i_n] != fgroundValue)
						{
							compute3dDistance(i_n, k, xDim, yDim, zDim, distance3DPtr, temp3dPtr, dim2D, h);
						}

						if (j == 0) //because if j is unsigned type, then j(0) - 1 == type maximum value
							break;
					}

					if (i == 0) //because if i is unsigned type, then i(0) - 1 == type maximum value
						break;
				}
			}
			sweepNumber++;
			break;

		case 2:
			for (k = 0; k < zDim; k++)
			{
				for (i = (xDim - 1); i >= 0; i--)
				{
					for (j = 0; j < yDim; j++)
					{
						i_n = i + j * xDim;
						if (tmpcurve3DPtr[k][i_n] != fgroundValue)
						{
							compute3dDistance(i_n, k, xDim, yDim, zDim, distance3DPtr, temp3dPtr, dim2D, h);
						}
					}

					if (i == 0) //because if i is unsigned type, then i(0) - 1 == type maximum value
						break;
				}
			}
			sweepNumber++;
			break;

		case 3:
			for (k = (zDim - 1); k >= 0; k--)
			{
				for (i = 0; i < xDim; i++)
				{
					for (j = 0; j < yDim; j++)
					{
						i_n = i + j * xDim;
						if (tmpcurve3DPtr[k][i_n] != fgroundValue)
						{
							compute3dDistance(i_n, k, xDim, yDim, zDim, distance3DPtr, temp3dPtr, dim2D, h);
						}
					}
				}
				if (k == 0) //because if k is unsigned type, then k(0) - 1 == type maximum value
					break;
			}
			sweepNumber++;
			break;
		case 4:
			for (k = (zDim - 1); k >= 0; k--)
			{
				for (i = 0; i < xDim; i++)
				{
					for (j = (yDim - 1); j >= 0; j--)
					{
						i_n = i + j * xDim;
						if (tmpcurve3DPtr[k][i_n] != fgroundValue)
						{
							compute3dDistance(i_n, k, xDim, yDim, zDim, distance3DPtr, temp3dPtr, dim2D, h);
						}

						if (j == 0) //because if j is unsigned type, then j(0) - 1 == type maximum value
							break;
					}
				}

				if (k == 0) //because if k is unsigned type, then k(0) - 1 == type maximum value
					break;
			}
			sweepNumber++;
			break;
		case 5:
			for (k = (zDim - 1); k >= 0; k--)
			{
				for (i = (xDim - 1); i >= 0; i--)
				{
					for (j = 0; j < yDim; j++)
					{
						i_n = i + j * xDim;
						if (tmpcurve3DPtr[k][i_n] != fgroundValue)
						{
							compute3dDistance(i_n, k, xDim, yDim, zDim, distance3DPtr, temp3dPtr, dim2D, h);
						}
					}

					if (i == 0) //because if i is unsigned type, then i(0) - 1 == type maximum value
						break;
				}

				if (k == 0) //because if k is unsigned type, then k(0) - 1 == type maximum value
					break;
			}
			sweepNumber++;
			break;
		case 6:
			for (k = (zDim - 1); k >= 0; k--)
			{
				for (i = (xDim - 1); i >= 0; i--)
				{
					for (j = (yDim - 1); j >= 0; j--)
					{
						i_n = i + j * xDim;
						if (tmpcurve3DPtr[k][i_n] != fgroundValue)
						{
							compute3dDistance(i_n, k, xDim, yDim, zDim, distance3DPtr, temp3dPtr, dim2D, h);
						}

						if (j == 0) //because if j is unsigned type, then j(0) - 1 == type maximum value
							break;
					}

					if (i == 0) //because if i is unsigned type, then i(0) - 1 == type maximum value
						break;
				}

				if (k == 0) //because if k is unsigned type, then k(0) - 1 == type maximum value
					break;
			}
			sweepNumber++;
			break;
		default:
			for (k = 0; k < zDim; k++)
			{
				for (i = 0; i < xDim; i++)
				{
					for (j = 0; j < yDim; j++)
					{
						i_n = i + j * xDim;
						if (tmpcurve3DPtr[k][i_n] != fgroundValue)
						{
							compute3dDistance(i_n, k, xDim, yDim, zDim, distance3DPtr, temp3dPtr, dim2D, h);
						}
					}
				}
			}
			sweepNumber++;
			break;
		}
	}
	for (i = 0; i < zDim; i++)
	{
		free(temp3dPtr[i]);
		free(tmpcurve3DPtr[i]);
	}

	free(temp3dPtr);
	free(tmpcurve3DPtr);
	temp3dPtr = NULL;
	tmpcurve3DPtr = NULL;
	return true;
}
//==============================================================================
void fSweeping3D(dataType ** distance3DPtr, dataType ** curve3DPtr, const size_t xDim, const size_t yDim, const size_t zDim, const dataType h, const dataType largeValue, const dataType fgroundValue, const ClipBox bestFit)
{
	size_t j, i, k, i_n;//i, j and k are loop counters. i_n is also a loop counter given by i_n = i + j * xDim
	size_t sweepNumber = 0, sweepDirection;
	// Within the best fit box
	//==============================================================================
	size_t k_min = bestFit.k_min, k_max = bestFit.k_max;
	size_t i_min = bestFit.i_min, i_max = bestFit.i_max;
	size_t j_min = bestFit.j_min, j_max = bestFit.j_max;
	//==============================================================================
	//const size_t dim2D = (i_max + 1) * (j_max + 1);
	const size_t dim2D = (xDim) * (yDim);
	dataType ** temp3dPtr = (dataType **)malloc(sizeof(dataType*) * (zDim));// temporary array used in the computation of distance
	dataType ** tmpcurve3DPtr = (dataType **)malloc(sizeof(dataType*) * zDim);// temporary array to copy curve3DPtr
	const size_t memory_alloc_size = sizeof(dataType) * dim2D;
	for (i = 0; i < zDim; i++) {
		temp3dPtr[i] = (dataType *)malloc(memory_alloc_size);
		tmpcurve3DPtr[i] = (dataType *)malloc(memory_alloc_size);
	}
	//==============================================================================
	// curve3DPtr is in the range 0-1 binary
	const double threshold = 0.5, range = 0.5 / 255.0;
	//Initialization stage
	for (k = 0; k < zDim; k++)
	{
		for (i = 0; i < xDim; i++)
		{
			for (j = 0; j < yDim; j++)
			{
				i_n = i + j * xDim;
				if (curve3DPtr[k][i_n] < threshold + range)
				{
					tmpcurve3DPtr[k][i_n] = distance3DPtr[k][i_n] = fgroundValue;
				}
				else
				{
					distance3DPtr[k][i_n] = largeValue;
					tmpcurve3DPtr[k][i_n] = curve3DPtr[k][i_n];
				}
			}
		}
	}
	//==============================================================================
	//
	int band = 1;

	i_max = min(i_max + band, xDim - 1);
	j_max = min(j_max + band, yDim - 1);
	k_max = min(k_max + band, zDim - 1);

	if (i_min > 0)
	{
		i_min = max(i_min - band, 0);
	}
	if (j_min > 0)
	{
		j_min = max(j_min - band, 0);
	}
	if (k_min > 0)
	{
		k_min = max(k_min - band, 0);
	}
	//

	//==============================================================================
	while (sweepNumber < 8)
	{
		sweepDirection = sweepNumber;
		switch (sweepDirection)
		{
		case 0:
			for (k = k_min; k < k_max + 1; k++)
			{
				for (i = i_min; i < i_max + 1; i++)
				{
					for (j = j_max; j >= j_min; j--)
					{
						i_n = i + j * xDim;
						if (tmpcurve3DPtr[k][i_n] != fgroundValue)
						{
							compute3dDistance(i_n, k, xDim, yDim, zDim, distance3DPtr, temp3dPtr, dim2D, h);
						}

						if (j == j_min) //because if j is unsigned type, then j(0) - 1 == type maximum value
							break;
					}
				}
			}
			sweepNumber++;
			break;

		case 1:
			for (k = k_min; k < k_max + 1; k++)
			{
				for (i = i_max; i >= i_min; i--)
				{
					for (j = j_max; j >= j_min; j--)
					{
						i_n = i + j * xDim;
						if (tmpcurve3DPtr[k][i_n] != fgroundValue)
						{
							compute3dDistance(i_n, k, xDim, yDim, zDim, distance3DPtr, temp3dPtr, dim2D, h);
						}

						if (j == j_min) //because if j is unsigned type, then j(0) - 1 == type maximum value
							break;
					}

					if (i == i_min) //because if i is unsigned type, then i(0) - 1 == type maximum value
						break;
				}
			}
			sweepNumber++;
			break;

		case 2:
			for (k = k_min; k < k_max + 1; k++)
			{
				for (i = i_max; i >= i_min; i--)
				{
					for (j = j_min; j < j_max + 1; j++)
					{
						i_n = i + j * xDim;
						if (tmpcurve3DPtr[k][i_n] != fgroundValue)
						{
							compute3dDistance(i_n, k, xDim, yDim, zDim, distance3DPtr, temp3dPtr, dim2D, h);
						}
					}

					if (i == i_min) //because if i is unsigned type, then i(0) - 1 == type maximum value
						break;
				}
			}
			sweepNumber++;
			break;

		case 3:
			for (k = k_max; k >= k_min; k--)
			{
				for (i = i_min; i < i_max + 1; i++)
				{
					for (j = j_min; j < j_max + 1; j++)
					{
						i_n = i + j * xDim;
						//i_n = i + j * xD;
						if (tmpcurve3DPtr[k][i_n] != fgroundValue)
						{
							compute3dDistance(i_n, k, xDim, yDim, zDim, distance3DPtr, temp3dPtr, dim2D, h);
						}
					}
				}
				if (k == k_min) //because if k is unsigned type, then k(0) - 1 == type maximum value
					break;
			}
			sweepNumber++;
			break;
		case 4:
			for (k = k_max; k >= k_min; k--)
			{
				for (i = i_min; i < i_max + 1; i++)
				{
					for (j = j_max; j >= j_min; j--)
					{
						i_n = i + j * xDim;
						if (tmpcurve3DPtr[k][i_n] != fgroundValue)
						{
							compute3dDistance(i_n, k, xDim, yDim, zDim, distance3DPtr, temp3dPtr, dim2D, h);
						}

						if (j == j_min) //because if j is unsigned type, then j(0) - 1 == type maximum value
							break;
					}
				}

				if (k == k_min) //because if k is unsigned type, then k(0) - 1 == type maximum value
					break;
			}
			sweepNumber++;
			break;
		case 5:
			for (k = k_max; k >= k_min; k--)
			{
				for (i = i_max; i >= i_min; i--)
				{
					for (j = j_min; j < j_max + 1; j++)
					{
						i_n = i + j * xDim;
						if (tmpcurve3DPtr[k][i_n] != fgroundValue)
						{
							compute3dDistance(i_n, k, xDim, yDim, zDim, distance3DPtr, temp3dPtr, dim2D, h);
						}
					}
					if (i == i_min) //because if i is unsigned type, then i(0) - 1 == type maximum value
						break;
				}
				if (k == k_min) //because if k is unsigned type, then k(0) - 1 == type maximum value
					break;
			}
			sweepNumber++;
			break;
		case 6:
			for (k = k_max; k >= k_min; k--)
			{
				for (i = i_max; i >= i_min; i--)
				{
					for (j = j_max; j >= j_min; j--)
					{
						i_n = i + j * xDim;
						//i_n = i + j * xD;
						if (tmpcurve3DPtr[k][i_n] != fgroundValue)
						{
							//compute3dDistance(i_n, k, xDim, yDim, zDim, distance3DPtr, temp3dPtr, dim2D, h);
							compute3dDistance(i_n, k, xDim, yDim, zDim, distance3DPtr, temp3dPtr, dim2D, h);
						}

						if (j == j_min) //because if j is unsigned type, then j(0) - 1 == type maximum value
							break;
					}

					if (i == i_min) //because if i is unsigned type, then i(0) - 1 == type maximum value
						break;
				}

				if (k == k_min) //because if k is unsigned type, then k(0) - 1 == type maximum value
					break;
			}
			sweepNumber++;
			break;
		default:
			for (k = k_min; k < k_max + 1; k++)
			{
				for (i = i_min; i < i_max + 1; i++)
				{
					for (j = j_min; j < j_max + 1; j++)
					{
						i_n = i + j * xDim;
						if (tmpcurve3DPtr[k][i_n] != fgroundValue)
						{
							compute3dDistance(i_n, k, xDim, yDim, zDim, distance3DPtr, temp3dPtr, dim2D, h);
						}
					}
				}
			}
			sweepNumber++;
			break;
		}
	}
	//==============================================================================
	for (i = 0; i < (k_max + 1); i++)
	{
		free(temp3dPtr[i]);
		free(tmpcurve3DPtr[i]);
	}
	free(temp3dPtr);
	free(tmpcurve3DPtr);
	temp3dPtr = NULL;
	tmpcurve3DPtr = NULL;
}
//==============================================================================
bool fastSweepingFunction_2D(dataType * distance2DPtr, dataType * curve2DPtr, const size_t xDim, const size_t yDim,
	const dataType h, const dataType largeValue, const dataType fgroundValue)
{
	size_t j, i, i_n;//i and k are loop counters. i_n is also a loop counter given by i_n = i + j * xDim
	size_t sweepNumber = 0, sweepDirection;
	const size_t dim2D = xDim * yDim;
	dataType * temp2dPtr = (dataType *)malloc(sizeof(dataType) * dim2D); //temporary array used in computation of distance

	//checks if the memory was allocated
	if (temp2dPtr == NULL)
		return false;

	if (distance2DPtr == NULL)
		return false;

	if (curve2DPtr == NULL)
		return false;

	//Initialization stage
	for (i = 0; i < xDim; i++)
	{
		for (j = 0; j < yDim; j++)
		{
			i_n = i + j * xDim;
			if (curve2DPtr[i_n] == fgroundValue)
			{
				distance2DPtr[i_n] = fgroundValue;
			}
			else
			{
				distance2DPtr[i_n] = largeValue;
			}
		}
	}

	//Main loop for four sweeps
	while (sweepNumber < 4)
	{
		sweepDirection = sweepNumber;
		switch (sweepDirection)
		{
		case 0:
			for (i = 0; i < xDim; i++)
			{
				for (j = 0; j < yDim; j++)
				{
					i_n = i + j * xDim;
					if (curve2DPtr[i_n] != fgroundValue)
					{
						compute2dDistance(i_n, xDim, yDim, distance2DPtr, temp2dPtr, dim2D, h);
					}
				}
			}
			sweepNumber++;
			break;

		case 1:
			for (i = 0; i < xDim; i++)
			{
				for (j = (yDim - 1); j >= 0; --j)
				{
					i_n = i + j * xDim;
					if (curve2DPtr[i_n] != fgroundValue)
					{
						compute2dDistance(i_n, xDim, yDim, distance2DPtr, temp2dPtr, dim2D, h);
					}
					if (j == 0) //because if i is unsigned type, then j(0) - 1 == type maximum value
						break;
				}
			}
			sweepNumber++;
			break;

		case 2:
			for (i = (xDim - 1); i >= 0; --i)
			{
				for (j = (yDim - 1); j >= 0; --j)
				{
					i_n = i + j * xDim;
					if (curve2DPtr[i_n] != fgroundValue)
					{
						compute2dDistance(i_n, xDim, yDim, distance2DPtr, temp2dPtr, dim2D, h);
					}
					if (j == 0) //because if i is unsigned type, then j(0) - 1 == type maximum value
						break;
				}
				if (i == 0) //because if i is unsigned type, then i(0) - 1 == type maximum value
					break;
			}
			sweepNumber++;
			break;

		default:
			for (i = (xDim - 1); i >= 0; --i)
			{
				for (j = 0; j < yDim; j++)
				{
					i_n = i + j * xDim;
					if (curve2DPtr[i_n] != fgroundValue)
					{
						compute2dDistance(i_n, xDim, yDim, distance2DPtr, temp2dPtr, dim2D, h);
					}
				}
				if (i == 0) //because if i is unsigned type, then i(0) - 1 == type maximum value
					break;
			}
			sweepNumber++;
			break;
		}
	}

	free(temp2dPtr);
	return true;
}

//Computation of distance in 3D
dataType compute3dDistance(size_t i_n, size_t k, const size_t xDim, const size_t yDim, const size_t zDim, dataType **distance3DPtr,
	dataType **temp3dPtr, const size_t dim2D, const dataType h)
{
	dataType a, b, c, a_1, a_2, a_3, lastCheck;
	if (k == 0)
	{
		if (i_n == 0)//if at the top left corner,
		{//use one sided difference.
			a = min(distance3DPtr[k][i_n], distance3DPtr[k][i_n + 1]);
			b = min(distance3DPtr[k][i_n], distance3DPtr[k][i_n + xDim]);
			c = min(distance3DPtr[k][i_n], distance3DPtr[k + 1][i_n]);
		}
		else if (i_n == (xDim - 1)) //if at the top right coner,
		{//use one sided difference.
			a = min(distance3DPtr[k][i_n - 1], distance3DPtr[k][i_n]);
			b = min(distance3DPtr[k][i_n], distance3DPtr[k][i_n + xDim]);
			c = min(distance3DPtr[k][i_n], distance3DPtr[k + 1][i_n]);
		}
		else if (i_n == ((yDim - 1) * xDim)) //if at the bottom left corner,
		{//use one sided difference.
			a = min(distance3DPtr[k][i_n], distance3DPtr[k][i_n + 1]);
			b = min(distance3DPtr[k][i_n - xDim], distance3DPtr[k][i_n]);
			c = min(distance3DPtr[k][i_n], distance3DPtr[k + 1][i_n]);
		}
		else if (i_n == (dim2D - 1)) //if at the bottom right corner,
		{// use one sided difference.
			a = min(distance3DPtr[k][i_n - 1], distance3DPtr[k][i_n]);
			b = min(distance3DPtr[k][i_n - xDim], distance3DPtr[k][i_n]);
			c = min(distance3DPtr[k][i_n], distance3DPtr[k + 1][i_n]);
		}
		else if ((i_n > 0) && (i_n < (xDim - 1)))//if at first row (excluding left and right corner),
		{//use one sided difference.
			a = min(distance3DPtr[k][i_n - 1], distance3DPtr[k][i_n + 1]);
			b = min(distance3DPtr[k][i_n], distance3DPtr[k][i_n + xDim]);
			c = min(distance3DPtr[k][i_n], distance3DPtr[k + 1][i_n]);
		}
		else if ((i_n > ((yDim - 1) * xDim)) && (i_n < (dim2D - 1))) //if at last row
		{//(excluding left and right corner), use one sided difference.
			a = min(distance3DPtr[k][i_n - 1], distance3DPtr[k][i_n + 1]);
			b = min(distance3DPtr[k][i_n - xDim], distance3DPtr[k][i_n]);
			c = min(distance3DPtr[k][i_n], distance3DPtr[k + 1][i_n]);
		}
		else if ((i_n != 0) && (i_n != ((yDim - 1) * xDim)) && ((i_n % xDim) == 0))//if at first column
		{//(excluding top and bottom left corner), use one sided difference.
			a = min(distance3DPtr[k][i_n], distance3DPtr[k][i_n + 1]);
			b = min(distance3DPtr[k][i_n - xDim], distance3DPtr[k][i_n + xDim]);
			c = min(distance3DPtr[k][i_n], distance3DPtr[k + 1][i_n]);
		}
		else if ((i_n != (xDim - 1)) && (i_n != (dim2D - 1)) && ((i_n % xDim) == (xDim - 1)))//if at last column
		{//(excluding top and bottom right corner), use one sided difference.
			a = min(distance3DPtr[k][i_n - 1], distance3DPtr[k][i_n]);
			b = min(distance3DPtr[k][i_n - xDim], distance3DPtr[k][i_n + xDim]);
			c = min(distance3DPtr[k][i_n], distance3DPtr[k + 1][i_n]);
		}
		else {//otherwise, use one sided difference.
			a = min(distance3DPtr[k][i_n - 1], distance3DPtr[k][i_n + 1]);
			b = min(distance3DPtr[k][i_n - xDim], distance3DPtr[k][i_n + xDim]);
			c = min(distance3DPtr[k][i_n], distance3DPtr[k + 1][i_n]);
		}
	}

	else if (k == (zDim - 1))
	{
		if (i_n == 0)//if at the top left corner,
		{//use one sided difference.
			a = min(distance3DPtr[k][i_n], distance3DPtr[k][i_n + 1]);
			b = min(distance3DPtr[k][i_n], distance3DPtr[k][i_n + xDim]);
			c = min(distance3DPtr[k - 1][i_n], distance3DPtr[k][i_n]);
		}
		else if (i_n == (xDim - 1)) //if at the top right coner,
		{//use one sided difference.
			a = min(distance3DPtr[k][i_n - 1], distance3DPtr[k][i_n]);
			b = min(distance3DPtr[k][i_n], distance3DPtr[k][i_n + xDim]);
			c = min(distance3DPtr[k - 1][i_n], distance3DPtr[k][i_n]);
		}
		else if (i_n == ((yDim - 1) * xDim)) //if at the bottom left corner,
		{//use one sided difference.
			a = min(distance3DPtr[k][i_n], distance3DPtr[k][i_n + 1]);
			b = min(distance3DPtr[k][i_n - xDim], distance3DPtr[k][i_n]);
			c = min(distance3DPtr[k - 1][i_n], distance3DPtr[k][i_n]);
		}
		else if (i_n == (dim2D - 1)) //if at the bottom right corner,
		{// use one sided difference.
			a = min(distance3DPtr[k][i_n - 1], distance3DPtr[k][i_n]);
			b = min(distance3DPtr[k][i_n - xDim], distance3DPtr[k][i_n]);
			c = min(distance3DPtr[k - 1][i_n], distance3DPtr[k][i_n]);
		}
		else if ((i_n > 0) && (i_n < (xDim - 1)))//if at first row (excluding left and right corner),
		{//use one sided difference.
			a = min(distance3DPtr[k][i_n - 1], distance3DPtr[k][i_n + 1]);
			b = min(distance3DPtr[k][i_n], distance3DPtr[k][i_n + xDim]);
			c = min(distance3DPtr[k - 1][i_n], distance3DPtr[k][i_n]);
		}
		else if ((i_n > ((yDim - 1) * xDim)) && (i_n < (dim2D - 1))) //if at last row
		{//(excluding left and right corner), use one sided difference.
			a = min(distance3DPtr[k][i_n - 1], distance3DPtr[k][i_n + 1]);
			b = min(distance3DPtr[k][i_n - xDim], distance3DPtr[k][i_n]);
			c = min(distance3DPtr[k - 1][i_n], distance3DPtr[k][i_n]);
		}
		else if ((i_n != 0) && (i_n != ((yDim - 1) * xDim)) && ((i_n % xDim) == 0))//if at first column
		{//(excluding top and bottom left corner), use one sided difference.
			a = min(distance3DPtr[k][i_n], distance3DPtr[k][i_n + 1]);
			b = min(distance3DPtr[k][i_n - xDim], distance3DPtr[k][i_n + xDim]);
			c = min(distance3DPtr[k - 1][i_n], distance3DPtr[k][i_n]);
		}
		else if ((i_n != (xDim - 1)) && (i_n != (dim2D - 1)) && ((i_n % xDim) == (xDim - 1)))//if at last column
		{//(excluding top and bottom right corner), use one sided difference.
			a = min(distance3DPtr[k][i_n - 1], distance3DPtr[k][i_n]);
			b = min(distance3DPtr[k][i_n - xDim], distance3DPtr[k][i_n + xDim]);
			c = min(distance3DPtr[k - 1][i_n], distance3DPtr[k][i_n]);
		}
		else {//otherwise, use one sided difference.
			a = min(distance3DPtr[k][i_n - 1], distance3DPtr[k][i_n + 1]);
			b = min(distance3DPtr[k][i_n - xDim], distance3DPtr[k][i_n + xDim]);
			c = min(distance3DPtr[k - 1][i_n], distance3DPtr[k][i_n]);
		}
	}

	else
	{
		if (i_n == 0)//if at the top left corner,
		{//use one sided difference.
			a = min(distance3DPtr[k][i_n], distance3DPtr[k][i_n + 1]);
			b = min(distance3DPtr[k][i_n], distance3DPtr[k][i_n + xDim]);
			c = min(distance3DPtr[k - 1][i_n], distance3DPtr[k + 1][i_n]);
		}
		else if (i_n == (xDim - 1)) //if at the top right coner,
		{//use one sided difference.
			a = min(distance3DPtr[k][i_n - 1], distance3DPtr[k][i_n]);
			b = min(distance3DPtr[k][i_n], distance3DPtr[k][i_n + xDim]);
			c = min(distance3DPtr[k - 1][i_n], distance3DPtr[k + 1][i_n]);
		}
		else if (i_n == ((yDim - 1) * xDim)) //if at the bottom left corner,
		{//use one sided difference.
			a = min(distance3DPtr[k][i_n], distance3DPtr[k][i_n + 1]);
			b = min(distance3DPtr[k][i_n - xDim], distance3DPtr[k][i_n]);
			c = min(distance3DPtr[k - 1][i_n], distance3DPtr[k + 1][i_n]);
		}
		else if (i_n == (dim2D - 1)) //if at the bottom right corner,
		{// use one sided difference.
			a = min(distance3DPtr[k][i_n - 1], distance3DPtr[k][i_n]);
			b = min(distance3DPtr[k][i_n - xDim], distance3DPtr[k][i_n]);
			c = min(distance3DPtr[k - 1][i_n], distance3DPtr[k + 1][i_n]);
		}
		else if ((i_n > 0) && (i_n < (xDim - 1)))//if at first row (excluding left and right corner),
		{//use one sided difference.
			a = min(distance3DPtr[k][i_n - 1], distance3DPtr[k][i_n + 1]);
			b = min(distance3DPtr[k][i_n], distance3DPtr[k][i_n + xDim]);
			c = min(distance3DPtr[k - 1][i_n], distance3DPtr[k + 1][i_n]);
		}
		else if ((i_n > ((yDim - 1) * xDim)) && (i_n < (dim2D - 1))) //if at last row
		{//(excluding left and right corner), use one sided difference.
			a = min(distance3DPtr[k][i_n - 1], distance3DPtr[k][i_n + 1]);
			b = min(distance3DPtr[k][i_n - xDim], distance3DPtr[k][i_n]);
			c = min(distance3DPtr[k - 1][i_n], distance3DPtr[k + 1][i_n]);
		}
		else if ((i_n != 0) && (i_n != ((yDim - 1) * xDim)) && ((i_n % xDim) == 0))//if at first column
		{//(excluding top and bottom left corner), use one sided difference.
			a = min(distance3DPtr[k][i_n], distance3DPtr[k][i_n + 1]);
			b = min(distance3DPtr[k][i_n - xDim], distance3DPtr[k][i_n + xDim]);
			c = min(distance3DPtr[k - 1][i_n], distance3DPtr[k + 1][i_n]);
		}
		else if ((i_n != (xDim - 1)) && (i_n != (dim2D - 1)) && ((i_n % xDim) == (xDim - 1)))//if at last column
		{//(excluding top and bottom right corner), use one sided difference.
			a = min(distance3DPtr[k][i_n - 1], distance3DPtr[k][i_n]);
			b = min(distance3DPtr[k][i_n - xDim], distance3DPtr[k][i_n + xDim]);
			c = min(distance3DPtr[k - 1][i_n], distance3DPtr[k + 1][i_n]);
		}
		else {//otherwise, use one sided difference.
			a = min(distance3DPtr[k][i_n - 1], distance3DPtr[k][i_n + 1]);
			b = min(distance3DPtr[k][i_n - xDim], distance3DPtr[k][i_n + xDim]);
			c = min(distance3DPtr[k - 1][i_n], distance3DPtr[k + 1][i_n]);
		}
	}

	a_1 = -max(-a, max(-b, -c)); // min
	a_3 = max(a, max(b, c)); // max
	a_2 = (a + b + c) - (a_3 + a_1);//mid
	temp3dPtr[k][i_n] = a_1 + h;

	if (temp3dPtr[k][i_n] < a_2)
	{
		distance3DPtr[k][i_n] = temp3dPtr[k][i_n];
	}
	else
	{
		if (fabs(a_1 - a_2) >= h)
		{
			temp3dPtr[k][i_n] = min(a_1, a_2) + h;
		}
		else
		{
			temp3dPtr[k][i_n] = (dataType)(a_1 + a_2 + sqrt((2 * h * h) - (pow((a_1 - a_2), 2)))) / 2;
		}
		if (temp3dPtr[k][i_n] < a_3)
		{
			distance3DPtr[k][i_n] = temp3dPtr[k][i_n];
		}
		else
		{
			lastCheck = 2 * ((a_1 * a_1) + (a_2 * a_2) + (a_3 * a_3) - ((a_1 * a_2) + (a_1 * a_3) + (a_2 * a_3)));
			if (lastCheck <= (3 * h * h))
			{
				temp3dPtr[k][i_n] = (dataType)(a_1 + a_2 + a_3 + sqrt((3 * h * h) - lastCheck)) / 3;
			}
			distance3DPtr[k][i_n] = temp3dPtr[k][i_n];
		}
	}
	return distance3DPtr[k][i_n];
}

//Computation of distance in 2D
dataType compute2dDistance(size_t i_n, const size_t xDim, const size_t yDim, dataType *distance2DPtr, dataType *temp2dPtr, const size_t dim2D,
	const dataType h)
{
	dataType a, b;
	if (i_n == 0)//if at the top left corner,
	{//use one sided difference.
		a = min(distance2DPtr[i_n], distance2DPtr[i_n + 1]);
		b = min(distance2DPtr[i_n], distance2DPtr[i_n + xDim]);
	}
	else if (i_n == (xDim - 1)) //if at the top right corner,
	{//use one sided difference.
		a = min(distance2DPtr[i_n - 1], distance2DPtr[i_n]);
		b = min(distance2DPtr[i_n], distance2DPtr[i_n + xDim]);
	}
	else if (i_n == ((yDim - 1) * xDim)) //if at the bottom left corner,
	{//use one sided difference.
		a = min(distance2DPtr[i_n], distance2DPtr[i_n + 1]);
		b = min(distance2DPtr[i_n - xDim], distance2DPtr[i_n]);
	}
	else if (i_n == (dim2D - 1)) //if at the bottom right corner,
	{// use one sided difference.
		a = min(distance2DPtr[i_n - 1], distance2DPtr[i_n]);
		b = min(distance2DPtr[i_n - xDim], distance2DPtr[i_n]);
	}
	else if ((i_n > 0) && (i_n < (xDim - 1)))//if at first row (excluding left and right corner),
	{//use one sided difference.
		a = min(distance2DPtr[i_n - 1], distance2DPtr[i_n + 1]);
		b = min(distance2DPtr[i_n], distance2DPtr[i_n + xDim]);
	}
	else if ((i_n > ((yDim - 1) * xDim)) && (i_n < (dim2D - 1))) //if at last row
	{//(excluding left and right corner), use one sided difference.
		a = min(distance2DPtr[i_n - 1], distance2DPtr[i_n + 1]);
		b = min(distance2DPtr[i_n - xDim], distance2DPtr[i_n]);
	}
	else if ((i_n != 0) && (i_n != ((yDim - 1) * xDim)) && ((i_n % xDim) == 0))//if at first column
	{//(excluding top and bottom left corner), use one sided difference.
		a = min(distance2DPtr[i_n], distance2DPtr[i_n + 1]);
		b = min(distance2DPtr[i_n - xDim], distance2DPtr[i_n + xDim]);
	}
	else if ((i_n != (xDim - 1)) && (i_n != (dim2D - 1)) && ((i_n % xDim) == (xDim - 1)))//if at last column
	{//(excluding top and bottom right corner), use one sided difference.
		a = min(distance2DPtr[i_n - 1], distance2DPtr[i_n]);
		b = min(distance2DPtr[i_n - xDim], distance2DPtr[i_n + xDim]);
	}
	else {//otherwise, use one sided difference.
		a = min(distance2DPtr[i_n - 1], distance2DPtr[i_n + 1]);
		b = min(distance2DPtr[i_n - xDim], distance2DPtr[i_n + xDim]);
	}

	// calculation of distance
	if (fabs(a - b) >= h)
	{
		temp2dPtr[i_n] = min(a, b) + h;
	}
	else
	{
		temp2dPtr[i_n] = (dataType)(a + b + sqrt((2 * h * h) - (pow((a - b), 2)))) / 2;
	}

	distance2DPtr[i_n] = min(distance2DPtr[i_n], temp2dPtr[i_n]);

	return distance2DPtr[i_n];
}