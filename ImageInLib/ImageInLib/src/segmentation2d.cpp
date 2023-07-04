/*
* Author : Konan Allaly
* Purpose : Updates for the INFLANET Project
*/

#include <string>
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
		return;

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
}