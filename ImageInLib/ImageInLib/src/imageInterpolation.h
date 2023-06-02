/*
* Author: Konan ALLALY
* Purpose: INFLANET project - Image Processing in Nuclear Medicine (2D/3D)
* Language:  C
*/
#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdbool.h>
#include "common_functions.h"
#include "transformation.h"


	bool nearestNeighborInterpolation(dataType ** originalImage, dataType ** newImage, size_t imageLength, size_t imageWidth, size_t imageHeight,
		                  dataType originalSpacing, dataType newSpacing);

	bool linear2dInterpolation(dataType** originalImage, dataType** newImage, size_t imageLength, size_t imageWidth, size_t imageHeight,
						  dataType originalSpacing, dataType newSpacing);

	bool downSampling(dataType** originalImage, dataType** newImage, size_t length, size_t width, size_t height);

	bool upSampling(dataType** originalImage, dataType** newImage, size_t length, size_t width, size_t height);

#ifdef __cplusplus
}
#endif