/*
* Author: Markjoe Olunna UBA
* Purpose: ImageInLife project - 4D Image Segmentation Methods
* Language:  C
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "setting_boundary_values.h"
#include "common_functions.h"

bool setBoundaryExactValues3D(dataType **inputDataArrayPtr, size_t length, size_t width, size_t height, dataType sphereRadius, dataType t, dataType h)
{
	size_t i, j, k;//loop counter for z dimension
	const size_t dim2D = length * width;
	dataType x, y, z;
	dataType hh = (h / 2.) + 1.25;
	size_t dim2Dminus = dim2D - 1;
	size_t widthminus = width - 1;
	size_t lengthminus = length - 1;
	size_t heightminus = height - 1;

	//checks if the memory was allocated
	if (inputDataArrayPtr == NULL)
		return false;

	//Setting of boundary values boundary values to exact values.
	for (k = 0; k < height; k++)// zDim
	{
		for (i = 0; i < length; i++)//xDim
		{
			for (j = 0; j < width; j++) //yDim
			{
				// 1D representation
				size_t x_n = x_new(i, j, length);
				x = (i * h) - hh;
				y = (j * h) - hh;
				z = (k * h) - hh;
				/*
				* (x_n == 0) - top left corner
				* (x_n == (length - 1)) - top right coner
				* (x_n == ((width - 1) * length)) - bottom left corner
				* (x_n == (dim2D - 1)) - bottom right corner
				* ((x_n > 0) && (x_n < (length - 1))) - first row (excluding left and right corner)
				* ((x_n >((width - 1) * length)) && (x_n < (dim2D - 1))) - last row (excluding left and right corner)
				* ((x_n != 0) && (x_n != ((width - 1) * length)) && ((i % length) == 0)) - first column (excluding top and bottom left corner)
				* ((x_n != (length - 1)) && (x_n != (dim2D - 1)) && ((x_n % length) == (length - 1))) - last column (excluding top and bottom right corner)
				*/
				if (((x_n == 0) || (x_n == lengthminus) || (x_n == (widthminus * length)) || (x_n == dim2Dminus) || ((x_n > 0) && (x_n < lengthminus))
						|| ((x_n >(widthminus * length)) && (x_n < dim2Dminus)) || ((x_n != 0) && (x_n != (widthminus * length)) && ((i % length) == 0))
						|| ((x_n != lengthminus) && (x_n != dim2Dminus) && ((x_n % length) == lengthminus)))
					|| k == 0 || k == heightminus)
				{
					inputDataArrayPtr[k][x_n] = ((pow(x, 2) + pow(y, 2)
						+ pow(z, 2) - pow(sphereRadius, 2)) / 4.) + t;
				}
			}
		}
	}
	return true;
}
	
bool setBoundaryToZeroDirichletBC(dataType **inputDataArrayPtr, size_t length, size_t width, size_t height)
{
	size_t i, j, k;//loop counter for z dimension
	const size_t dim2D = length * width;
	size_t dim2Dminus = dim2D - 1;
	size_t widthminus = width - 1;
	size_t lengthminus = length - 1;
	size_t heightminus = height - 1;

	//checks if the memory was allocated
	if (inputDataArrayPtr == NULL)
		return false;

	// Setting of boundary values boundary values to zero.
	for (k = 0; k < height; k++)// zDim
	{
		for (i = 0; i < length; i++)//xDim
		{
			for (j = 0; j < width; j++) //yDim
			{
				// 1D representation
				size_t x_n = x_new(i, j, length);
				/*
				* (x_n == 0) - top left corner
				* (x_n == (length - 1)) - top right coner
				* (x_n == ((width - 1) * length)) - bottom left corner
				* (x_n == (dim2D - 1)) - bottom right corner
				* ((x_n > 0) && (x_n < (length - 1))) - first row (excluding left and right corner)
				* ((x_n >((width - 1) * length)) && (x_n < (dim2D - 1))) - last row (excluding left and right corner)
				* ((x_n != 0) && (x_n != ((width - 1) * length)) && ((i % length) == 0)) - first column (excluding top and bottom left corner)
				* ((x_n != (length - 1)) && (x_n != (dim2D - 1)) && ((x_n % length) == (length - 1))) - last column (excluding top and bottom right corner)
				*/
				if (((x_n == 0) || (x_n == lengthminus) || (x_n == (widthminus * length)) || (x_n == dim2Dminus) || ((x_n > 0) && (x_n < lengthminus))
					|| ((x_n >(widthminus * length)) && (x_n < dim2Dminus)) || ((x_n != 0) && (x_n != (widthminus * length)) && ((i % length) == 0))
					|| ((x_n != lengthminus) && (x_n != dim2Dminus) && ((x_n % length) == lengthminus)))
					|| k == 0 || k == heightminus)
				{
					inputDataArrayPtr[k][x_n] = 0.;
				}
			}
		}
	}
	return true;
}
