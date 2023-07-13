#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdio.h>
#include <stdbool.h>
#include "common_functions.h"
#include "transformation.h"

	//Interpolation regarding z-Direction
	bool nearestNeighborInterpolation(dataType ** originalImage, dataType ** newImage, size_t imageLength, size_t imageWidth, size_t imageHeight,
		                  dataType originalSpacing, dataType newSpacing);

	bool linear2dInterpolation(dataType** originalImage, dataType** newImage, size_t imageLength, size_t imageWidth, size_t imageHeight,
						  dataType originalSpacing, dataType newSpacing);

	bool downSampling(dataType** originalImage, dataType** newImage, size_t length, size_t width, size_t height);

	bool upSampling(dataType** originalImage, dataType** newImage, size_t length, size_t width, size_t height);

	//=========================

	//Structure to handle image spacing
	typedef struct {
		dataType sx, sy, sz;
	} Spacing3D;

	//Structure to handle image coordinates
	//We need new structure because coordinates are array indexes (integer type is needed)
	typedef struct {
		size_t x, y, z;
	} imgPoint3D;

	/*
	* Get real world cordinate from image coordinate
	* srcPoint : contains the voxel indexes
	* realOrigin : image origin in real world
	* imageSpacing : - distance between voxels in x and y direction
	*                - distance between slices
	* orientation : IJK(1,1,1), RAS(-1,-1,1) or LPS(1,-1,-1)
	*/
	Point3D imageCoordToRealCoord(imgPoint3D srcPoint, Point3D realOrigin, Spacing3D imageSpacing, Point3D orientation);

	//Get image coordinate from real coordinate
	imgPoint3D realCoordToImageCoord(Point3D srcPoint, Point3D realOrigin, Spacing3D imageSpacing, Point3D orientation);

	bool resizeImage(dataType* oldImage, dataType* newImage);

#ifdef __cplusplus
}

#endif