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

	//=====================================================

	/*
	 * Get 2D point in real world Cordinates System from 2D image Coordinates System
	 * Three operations are done here : scaling, translation and rotation
	*/
	Point2D getRealCoordFromImageCoord2D(Point2D srcPoint, Point2D realOrigin, PixelSpacing imageSpacing, OrientationMatrix2D orientation);

	/*
	* Get 2D point in image Coordinates System from 2D point in real Coordinates System
	* Three operations are done here : scaling, translation and rotation
	*/
	Point2D getImageCoordFromRealCoord2D(Point2D srcPoint, Point2D imageOrigin, PixelSpacing imageSpacing, OrientationMatrix2D orientation);

	/*
	* This function find the four neighbors of given point
	* PointNeighbors2D findPointNeighbors(Point2D point, size_t Xmax, size_t Ymax)
	* point : the given point
	* Xmax : maximum in x direction (height)
	* Ymax : maximum in y direction (width)
	*/
	PointNeighbors2D findPointNeighbors(Point2D point, const size_t Xmax, const size_t Ymax);

	/*
	* This functions check if two given points are similar
	* point1 and point2 are points to be tested
	*/
	bool ArePointsTheSame(Point2D point1, Point2D point2);

	/*
	* This function compute 2D euclidian distance between two given points
	* point1 and point2 are given points
	*/
	dataType getDistance2D(Point2D point1, Point2D point2);

	/*
	* This function return the nearest neighbor of given point
	* Point2D findNearestNeighbor2d(Point2D point, const size_t Xmax, const size_t Ymax)
	* point : the given point
	* Xmax : maximum in x direction (height)
	* Ymax : maximum in y direction (width)
	*/
	Point2D findNearestNeighbor2d(Point2D point, const size_t Xmax, const size_t Ymax);

	/*
	* This function perform image interpolation from image coordinates system to real world coordinate system
	* resizeImageFromImageCoordToRealCoord(Image_Data2D src_image, Image_Data2D dest_image, OrientationMatrix2D orientation)
	* Image_Data2D src_image : structure to handle image in "IMAGE Cordinates System"
	*                          - height, width : image dimension
	*                          - imageDataPtr : pointer for pixels value
	*                          - origin : image origin
	*                          - spacing : pixel size
	* Image_Data2D dest_image : structure to handle interpolated image in "Real World Coordinates System"
	*                          - height, width : interpolated image dimension
	*                          - imageDataPtr : pointer for interpolated pixels value
	*                          - origin : interpolated image origin
	*                          - spacing : interpolated pixel size
	* OrientationMatrix2D orientation : orientation of the rotation
	*/
	bool resizeImageFromImageCoordToRealCoord(Image_Data2D src_image, Image_Data2D dest_image, OrientationMatrix2D orientation);

#ifdef __cplusplus
}

#endif