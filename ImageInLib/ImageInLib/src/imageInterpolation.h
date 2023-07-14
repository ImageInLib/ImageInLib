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

	/*
	* srcPoint : source 3d points in image coordinates system
	* ImageOrigin : image origin in real world (ranslation)
	* VoxelSpacing : - distance between voxels in x and y direction
	*                - distance between slices (scaling)
	* orientation : matrix for orientation/rotation IJK((1, 0, 0),(0,1,0),(0,0,1)), RAS((-1,0,0),(0,-1,0),(0,0,1))
	* or LPS((0,0,1),(0,-1,0),(0,0,-))
	*/
	Point3D imageCoordToRealCoord(Point3D srcPoint, Point3D realOrigin, VoxelSpacing imageSpacing, OrientationMatrix orientation);

	/*
	* srcPoint : source 3d points in real world coordinates system
	* ImageOrigin : image origin in image coordinate system
	* VoxelSpacing : - distance between voxels in x and y direction
	*                - distance between slices (scaling)
	* orientation : matrix for orientation/rotation IJK((1, 0, 0),(0,1,0),(0,0,1)), RAS((-1,0,0),(0,-1,0),(0,0,1))
	* or LPS((0,0,1),(0,-1,0),(0,0,-))
	*/
	Point3D realCoordToImageCoord(Point3D srcPoint, Point3D realOrigin, VoxelSpacing imageSpacing, OrientationMatrix orientation);

	/*
	* oldImage to handle the input image (height, width and pointer for pixel value)
	* newImage to handle the result image (height, width and pointer for pixel value)
	*/
	bool resizeImage(Image_Data2D oldImage, Image_Data2D newImage);

#ifdef __cplusplus
}

#endif