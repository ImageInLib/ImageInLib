/*
* Author : Konan Allaly
* Purpose : Updates for the INFLANET Project
*/

#include <iostream>
#include <climits>
#include <crtdbg.h>
#include <corecrt_malloc.h>
#include <cmath>

#include "segmentation2d.h"

dataType min(dataType a, dataType b) {
	if (a > b) {
		return b;
	}
	else {
		return a;
	}
}

bool rescaleToZeroOne2d(dataType* imageDataPtr, const size_t height, const size_t width)
{
	//check if the memory was allocated successfully
	if (imageDataPtr == NULL)
		return false;

	size_t i;
	dataType max = 0, min = 100000, quotient, offset;

	//Determine minimum and maximum value
	for (i = 0; i < height * width; i++) {
		if (imageDataPtr[i] < min)
			min = imageDataPtr[i];
		if (imageDataPtr[i] > max)
			max = imageDataPtr[i];
	}

	quotient = (dataType)1. / (max - min);
	offset = min * quotient;
	//Rescale values to interval (0, 1)
	for (i = 0; i < height * width; i++) {
		imageDataPtr[i] = (quotient * imageDataPtr[i] - offset);
	}

	return true;
}

bool generateInitialSegmentationFunction(dataType* imageDataPtr, const size_t height, const size_t width, point2D* center, dataType v, dataType R)
{
	size_t i, j;
	int dx, dy;
	dataType norm_of_distance = 0.0, new_value = 0.0;

	if (imageDataPtr == NULL)
		return false;

	for (i = 0; i < height; i++) {
		dx = i - center->x;
		for (j = 0; j < width; j++) {
			dy = j - center->y;
			norm_of_distance = sqrt(dx * dx + dy * dy);
			new_value = (dataType)((1.0 / (norm_of_distance + v)) - (1.0 / (R + v)));
			if (norm_of_distance > R) {
				imageDataPtr[x_new(i, j, height)] = 0;
			}
			else {
				imageDataPtr[x_new(i, j, height)] = new_value;
			}
		}
	}
	rescaleToZeroOne2d(imageDataPtr, height, width);

	return true;
}

bool set2dDirichletBoundaryCondition(dataType* imageDataPtr, const size_t height, const size_t width) {
	size_t i, j;
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			if (i == 0 || i == height - 1 || j == 0 || j == width - 1) {
				imageDataPtr[x_new(i, j, height)] = 0.0;
			}
		}
	}
	return true;
}

bool computeNormOfGradientDiamondCells(dataType* arrayPtr, neighPtrs neigbours, const size_t height, const size_t width, dataType h) {

	size_t i, j, i_ext, j_ext, currentIndx;
	const size_t height_ext = height + 2, width_ext = width + 2;
	size_t dim2D_ext = height_ext * width_ext;

	dataType* extendedArray = new dataType[dim2D_ext];
	if (extendedArray == NULL)
		return false;

	copyDataTo2dExtendedArea(arrayPtr, extendedArray, height, width);
	reflection2D(extendedArray, height_ext, width_ext);

	dataType uP, uN, uNW, uNE, uS, uSW, uSE, uW, uE;
	dataType ux, uy;

	for (i = 0, i_ext = 1; i < height; i++, i_ext++) {
		for (j = 0, j_ext = 1; j < width; j++, j_ext++) {

			size_t iplus = i_ext + 1;
			size_t iminus = i_ext - 1;
			size_t jplus = j_ext + 1;
			size_t jminus = j_ext - 1;

			currentIndx = x_new(i, j, height);
			uP = extendedArray[x_new(i_ext, j_ext, height_ext)];
			uE = extendedArray[x_new(iplus, j_ext, height_ext)];
			uW = extendedArray[x_new(iminus, j_ext, height_ext)];
			uN = extendedArray[x_new(i_ext, jminus, height_ext)];
			uS = extendedArray[x_new(i_ext, jplus, height_ext)];
			uNE = extendedArray[x_new(iplus, jminus, height_ext)];
			uNW = extendedArray[x_new(iminus, jminus, height_ext)];
			uSE = extendedArray[x_new(iplus, jplus, height_ext)];
			uSW = extendedArray[x_new(iminus, jplus, height_ext)];

			//East
			ux = (uE - uP) / h;
			uy = (uNE + uN - uS - uSE) / (4.0 * h);
			neigbours.East[currentIndx] = sqrt(ux * ux + uy * uy);

			//West
			ux = (uP - uW) / h;
			uy = (uNW + uN - uSW - uS) / (4.0 * h);
			neigbours.West[currentIndx] = sqrt(ux * ux + uy * uy);

			//North
			ux = (uNE + uE - uNW - uW) / (4.0 * h);
			uy = (uN - uP) / h;
			neigbours.North[currentIndx] = sqrt(ux * ux + uy * uy);

			//South
			ux = (uSE + uE - uSW - uW) / (4.0 * h);
			uy = (uP - uS) / h;
			neigbours.South[currentIndx] = sqrt(ux * ux + uy * uy);
		}
	}

	delete[] extendedArray;

	return true;
}

bool epsilonRegularization(neighPtrs neighbours, const size_t height, const size_t width, dataType epsilon) {
	size_t i, j, currentIndx;
	dataType current = 0.0;
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			currentIndx = x_new(i, j, height);

			current = neighbours.East[currentIndx];
			neighbours.East[currentIndx] = (dataType)(sqrt(current * current + epsilon));

			current = neighbours.West[currentIndx];
			neighbours.West[currentIndx] = (dataType)(sqrt(current * current + epsilon));

			current = neighbours.North[currentIndx];
			neighbours.North[currentIndx] = (dataType)(sqrt(current * current + epsilon));

			current = neighbours.South[currentIndx];
			neighbours.South[currentIndx] = (dataType)(sqrt(current * current + epsilon));
		}
	}
	return true;
}