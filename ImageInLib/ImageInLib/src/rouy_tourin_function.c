/*
* Author: Markjoe Olunna UBA
* Purpose: ImageInLife project - 4D Image Segmentation Methods
* Language:  C
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "data_initialization.h"
#include "distance_function.h"
#include "common_functions.h"


bool rouyTourinFunction_3D(dataType ** distance3DPtr, dataType ** image3DPtr, dataType tolerance,
	const size_t xDim, const size_t yDim, const size_t zDim, dataType tau, const dataType h)
{
	if (distance3DPtr == NULL || image3DPtr == NULL)
		return false;

	size_t k, i, i_n, iter;//i and k are loop counters. i_n is also a loop counter given by 

	size_t rowDim_ext = xDim + 2;
	size_t columnDim_ext = yDim + 2;
	size_t sliceDim_ext = zDim + 2;

	size_t origRowDim = xDim;
	size_t origColumnDim = yDim;

	size_t sliceBound = (rowDim_ext - 1)* columnDim_ext;
	size_t k_o, i_o, j_o, j;
	const size_t noIteration = 1000;//(size_t)(sqrt((double)(xDim * xDim + yDim * yDim + zDim * zDim)) / tau);

	const size_t dim2D_ext = rowDim_ext * columnDim_ext;
	dataType tauu = tau / h, mass = 10.0;
	dataType current_dist;

	dataType ** zPtrtemp = (dataType **)malloc(sizeof(dataType*) * sliceDim_ext);

	//checks if the memory was allocated
	if (zPtrtemp == NULL)
		return false;

	for (i = 0; i < sliceDim_ext; i++)
	{
		zPtrtemp[i] = (dataType *)malloc(sizeof(dataType) * dim2D_ext);

		//checks if the memory was allocated
		if (zPtrtemp[i] == NULL)
			return false;
	}

	// filling the 3D array with number=0
	initialize3dArrayD(distance3DPtr, xDim, yDim, zDim, 0);

	iter = 0;

	// implementation of Rouy Tourin scheme in 3D
	while (mass > tolerance && iter < noIteration)
	{
		mass = 0;
		iter++;

		copyDataToExtendedArea(distance3DPtr, zPtrtemp, zDim, origRowDim, origColumnDim);
		//reflection of zPtrtemp 
		reflection3D(zPtrtemp, sliceDim_ext, rowDim_ext, columnDim_ext);

		//computation of 3D distance
		k_o = 0;
		for (k = 1; k <= zDim; k++, k_o++)
		{
			i_o = 0;
			for (i = 1; i <= xDim; i++, i_o++)//row loop
			{
				j_o = 0;
				for (j = 1; j <= yDim; j++, j_o++)// column loop i_n = i(row) + j(column) * rowDim
				{
					i_n = i + j * rowDim_ext;
					if (image3DPtr[k_o][(i_o + j_o * xDim)])
					{
						current_dist = zPtrtemp[k][i_n];

						distance3DPtr[k_o][(i_o + j_o * xDim)] = zPtrtemp[k][i_n] + tau - tauu *
							((dataType)sqrt(
								max(pow(min(zPtrtemp[k][i_n - 1] - current_dist, 0), 2),
									pow(min(zPtrtemp[k][i_n + 1] - current_dist, 0), 2)) +
								max(pow(min(zPtrtemp[k][i_n - rowDim_ext] - current_dist, 0), 2),
									pow(min(zPtrtemp[k][i_n + rowDim_ext] - current_dist, 0), 2)) +
								max(pow(min(zPtrtemp[k - 1][i_n] - current_dist, 0), 2),
									pow(min(zPtrtemp[k + 1][i_n] - current_dist, 0), 2))));

						mass += (dataType)pow(distance3DPtr[k_o][(i_o + j_o * xDim)] - zPtrtemp[k][i_n], 2);
					}
				}
			}
		}

		mass = (dataType)sqrt(mass);
	}

	for (i = 0; i < sliceDim_ext; i++)
		free(zPtrtemp[i]);

	free(zPtrtemp);
	return true;
}

