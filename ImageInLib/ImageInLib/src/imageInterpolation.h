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

	//Change coordinate system
	bool imageCoordToRealCoord(Point3D srcPoint, Spacing imageSpacing, Point3D destPoint);

	bool IJKtoRAS(Point3D srcPoint, Spacing imageSpacing, Point3D destPoint);

	bool IJKtoLPS(Point3D srcPoint, Spacing imageSpacing, Point3D destPoint);

	//2D images for test

	typedef struct {
		dataType sx, sy;
	} Spacing2D;

	typedef struct {
		Point2D origin;
		size_t height, width;
		dataType* dataPtr;
		Spacing2D toRealCoordinates;
	} patientImageData2D;

	bool interpolateToRealDimension2D(patientImageData2D imageSrc, Spacing2D newSpacing, const char* outputPathPtr);

#ifdef __cplusplus
}

#endif