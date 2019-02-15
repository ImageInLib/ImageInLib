/*
* Author: Markjoe Olunna UBA
* Purpose: ImageInLife project - 4D Image Segmentation Methods
* Language:  C
*/
#include <stdlib.h>
#include <math.h>
#include "data_initialization.h"
#include "distance_function.h"
#include "common_functions.h"


bool bruteForceFunction_3D(dataType ** distance3DPtr, dataType ** curve3DPtr,
	const size_t xDim, const size_t yDim, const size_t zDim, const dataType largeValue,
	const dataType fgroundValue)
{
	size_t k_1, i_1, k_2, i_2;//loop counter for z dimension
	const size_t dim2D = xDim * yDim;
	dataType dx, dy, dz;
	dataType dist;

	//checks if the memory was allocated
	if (distance3DPtr == NULL)
		return false;

	//checks if the memory was allocated
	if (curve3DPtr == NULL)
		return false;

	// filling the 3D array with number=value
	initialize3dArrayD(distance3DPtr, xDim, yDim, zDim, largeValue);

	Point3D * surface_points = malloc(sizeof(Point3D) * zDim * dim2D);
	size_t ptsNum = 0;

	for (k_2 = 0; k_2 < zDim; k_2++) // z axis of the input surface or image
	{
		for (i_2 = 0; i_2 < dim2D; i_2++)// x-y axis of the input surface or image
		{
			if (curve3DPtr[k_2][i_2] == fgroundValue)
			{
				surface_points[ptsNum].x = (int)(i_2 / xDim);
				surface_points[ptsNum].y = (dataType)(i_2 % xDim);
				surface_points[ptsNum].z = (dataType)k_2;
				ptsNum++;
			}
		}
	}


	for (k_1 = 0; k_1 < zDim; k_1++) // z axis of the original image
	{
		for (i_1 = 0; i_1 < dim2D; i_1++)// x-y axis of the original image 
		{
			for (size_t i = 0; i < ptsNum; i++)
			{
				dz = k_1 - surface_points[i].z; // difference between z axes of both images

				dx = (int)(i_1 / xDim) - surface_points[i].x;// difference between x axes of both images
				dy = (i_1 % xDim) - surface_points[i].y;// difference between y axes of both images
				dist = dx * dx + dy * dy + dz * dz;

				if (dist <= distance3DPtr[k_1][i_1]) {
					distance3DPtr[k_1][i_1] = dist;
				}
			}
		}
	}

	for (k_1 = 0; k_1 < zDim; k_1++) // z axis of the original image
	{
		for (i_1 = 0; i_1 < dim2D; i_1++)// x-y axis of the original image 
		{
			distance3DPtr[k_1][i_1] = sqrt(distance3DPtr[k_1][i_1]);
		}
	}

	free(surface_points);

	return true;
}
