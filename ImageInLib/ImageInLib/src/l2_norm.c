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

bool l2norm3dDataArrayD(dataType ** dataArray3DPtr1, dataType ** dataArray3DPtr2, unsigned char * pathPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType h,
	size_t step)
{
	size_t i, k;//loop counter for z dimension
	const size_t dim2D = xDim * yDim;
	dataType norm, sumPower = 0;
	dataType tau = h * h;
	dataType hhh = h * h * h;

	FILE * outputfile; //file stream

					   //checks if the memory was allocated
	if (dataArray3DPtr1 == NULL || dataArray3DPtr2 == NULL)
		return false;

	/*Computation of norm*/
	for (k = 0; k < zDim; k++)
	{
		for (i = 0; i < dim2D; i++)
		{
			sumPower += (pow(dataArray3DPtr1[k][i] - dataArray3DPtr2[k][i], 2) * hhh);
		}
	}
	norm = sqrt(sumPower);

	//checks if the file was sucessfully opened
	if ((fopen_s(&outputfile, pathPtr, "a+")) != 0) {
		return false;
	}
	else
	{
		fprintf(outputfile, "L2 norm in time %f with h = %f is %f\n", step * tau, h, norm);
	}
	fclose(outputfile);
	return true;
}

dataType l2normD(double ** dataArray3DPtr1, double ** dataArray3DPtr2, const size_t xDim, const size_t yDim, const size_t zDim, double h)
{
	size_t i, k;//loop counter for z dimension
	const size_t dim2D = xDim * yDim;
	double l2norm, sumPower = 0;
	double hhh = h * h * h;

	//checks if the memory was allocated
	if (dataArray3DPtr1 == NULL || dataArray3DPtr2 == NULL)
		return false;

	/*Computation of norm*/
	for (k = 0; k < zDim; k++)
	{
		for (i = 0; i < dim2D; i++)
		{
			sumPower += (pow(dataArray3DPtr1[k][i] - dataArray3DPtr2[k][i], 2) * hhh);
		}
	}
	l2norm = sqrt(sumPower);

	return l2norm;
}