#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdbool.h>
#include "common_functions.h"
#include "transformation.h"

	typedef struct {
		dataType sx, sy, sz;
	} Spacing;

	typedef struct {
		Point3D origin;
		size_t height, width, length;
		dataType** dataPtr;
		Spacing toRealCoordinates;
	} patientImageData;

	//Interpolation regarding z-Direction
	bool nearestNeighborInterpolation(dataType ** originalImage, dataType ** newImage, size_t imageLength, size_t imageWidth, size_t imageHeight,
		                  dataType originalSpacing, dataType newSpacing);

	bool linear2dInterpolation(dataType** originalImage, dataType** newImage, size_t imageLength, size_t imageWidth, size_t imageHeight,
						  dataType originalSpacing, dataType newSpacing);
	//====================

	bool downSampling(dataType** originalImage, dataType** newImage, size_t length, size_t width, size_t height);

	bool upSampling(dataType** originalImage, dataType** newImage, size_t length, size_t width, size_t height);

	//Interpolation in all directions
	bool interpolateToRealDimension(patientImageData imageSrc, const char * outputPathPtr);

#ifdef __cplusplus
}

#endif