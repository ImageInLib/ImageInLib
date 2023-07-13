#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdio.h>
#include <stdbool.h>
#include "common_functions.h"
#include "transformation.h"

	//Interpolation regarding z-Direction
	bool nearestNeighborInterpolation(dataType** originalImage, dataType** newImage, size_t imageLength, size_t imageWidth, size_t imageHeight,
		dataType originalSpacing, dataType newSpacing);

	bool linear2dInterpolation(dataType** originalImage, dataType** newImage, size_t imageLength, size_t imageWidth, size_t imageHeight,
		dataType originalSpacing, dataType newSpacing);

	bool downSampling(dataType** originalImage, dataType** newImage, size_t length, size_t width, size_t height);

	bool upSampling(dataType** originalImage, dataType** newImage, size_t length, size_t width, size_t height);

	//=========================

	//Structure to handle image spacing
	typedef struct {
		dataType sx, sy, sz;
	} voxelSpacing;

	typedef struct {
		Point3D v1, v2, v3;
	}orientationMatrix;

	/*
	* Get real world cordinate from image coordinate
	* srcPoint : contains the voxel indexes
	* realOrigin : image origin in real world
	* imageSpacing : - distance between voxels in x and y direction
	*                - distance between slices
	* orientation : IJK((1, 0, 0),(0,1,0),(0,0,1)), RAS((-1,0,0),(0,-1,0),(0,0,1))
	* or LPS((0,0,1),(0,-1,0),(0,0,-))
	*/
	Point3D imageCoordToRealCoord(Point3D srcPoint, Point3D realOrigin, voxelSpacing imageSpacing, orientationMatrix orientation);

	//Get image coordinate from real coordinate
	//Point3D realCoordToImageCoord(Point3D srcPoint, Point3D realOrigin, voxelSpacing3D imageSpacing, orientationMatrix orientation);

	typedef struct {
		dataType sx, sy;
	} pixelSpacing;

	bool resizeImage(Image_Data2D oldImage, Image_Data2D newImage);

#ifdef __cplusplus
}

#endif