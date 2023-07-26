#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdio.h>
#include <stdbool.h>
#include "common_functions.h"
#include "transformation.h"
#include "interpolations.h"

	typedef enum
	{
		NEAREST_NEIGHBOR = 1,
		BILINEAR, // for 2D image
		TRILINEAR // for 3D image
	} interpolationMethod;

	//====================
	//3D Functions

	/*
	* Point3D getImageCoordFromRealCoord3D(Point3D src_point, Point3D real_origin, VoxelSpacing real_spacing, OrientationMatrix orientation)
	* src_point : source 3D points in image coordinates system
	* real_origin : origin in real world
	* real_spacing : - distance between voxels in x and y direction
	*                 - distance between slices
	* orientation : matrix for orientation/rotation IJK((1, 0, 0),(0,1,0),(0,0,1)), RAS((-1,0,0),(0,-1,0),(0,0,1))
	* or LPS((0,0,1),(0,-1,0),(0,0,-1))
	* For given point in real world cordinates system the function provide
	* the corresponding point in image coordinates system
	*/
	Point3D getRealCoordFromImageCoord3D(Point3D image_point, Point3D real_origin, VoxelSpacing real_spacing, OrientationMatrix orientation);

	/*
	* Point3D getRealCoordFomImageCoord3D(Point3D src_point, Point3D real_origin, VoxelSpacing real_spacing, OrientationMatrix orientation)
	* src_point : source 3D points in image coordinates system
	* real_origin : origin in real world
	* real_spacing : - distance between voxels in x and y direction
	*                 - distance between slices
	* orientation : matrix for orientation/rotation IJK((1, 0, 0),(0,1,0),(0,0,1)), RAS((-1,0,0),(0,-1,0),(0,0,1))
	* or LPS((0,0,1),(0,-1,0),(0,0,-1))
	* For given point in image cordinates system the function provide
	* the corresponding point in real world coordinates system
	*/
	Point3D getImageCoordFromRealCoord3D(Point3D real_point, Point3D real_origin, VoxelSpacing real_spacing, OrientationMatrix orientation);

	/*
	* This function return the interpolated value by nearest neighbor approach
	* dataType getInterpolatedValueNearestNeighbor3D(Image_Data src_image, Point3D point)
	* Image_Data src_image : source image
	* Point3D point : the current point
	*/
	dataType getInterpolatedValueNearestNeighbor3D(Image_Data src_image, Point3D point);

	/*
	* This function return the interpolated value by bilinear interpolation approach
	* dataType getInterpolatedValueTrilinear3D(Image_Data src_image, Point3D point)
	* Image_Data src_image : source image
	* Point3D point : the current point
	*/
	dataType getInterpolatedValueTrilinear3D(Image_Data src_image, Point3D point);

	/*
	* This function perform image interpolation
	* bool imageInterpolation3D(Image_Data src_image, Image_Data dest_image, dataType scale_factor, interpolationMethod method)
	* Image_Data src_image : source image structure
	*                          - height, lenght, width : image dimension
	*                          - imageDataPtr : pointer for pixels value
	*                          - origin : image origin
	*                          - spacing : pixel size
	* Image_Data dest_image : interpolated image structure
	*                          - height, Lenght, width : interpolated image dimension
	*                          - imageDataPtr : pointer for interpolated pixels value
	*                          - origin : interpolated image origin
	*                          - spacing : interpolated pixel size
	* interpolationMethod method : nearest neighbor, bilinear (for 2D image), trilinear
	*/
	bool imageInterpolation3D(Image_Data src_image, Image_Data dest_image, interpolationMethod method);

	//=====================================================
	//2D Functions

	/*
	* Point2D getImageCoordFromRealCoord2D(Point2D src_point, Point2D real_origin, PixelSpacing real_spacing, OrientationMatrix2D orientation)
	* Point2D src_point : source 2D points in image coordinates system
	* Point2D real_origin : origin in real world coordinates system
	* PixelSpacing real_spacing : distance between voxels in x and y direction
	* orientation : matrix for orientation
	* For given point in image cordinates system, the function provides
	* the corresponding point in image real world coordinates system
	*/
	Point2D getRealCoordFromImageCoord2D(Point2D src_point, Point2D real_origin, PixelSpacing real_spacing, OrientationMatrix2D orientation);

	/*
	* Point2D getImageCoordFromRealCoord2D(Point2D src_point, Point2D real_origin, PixelSpacing real_spacing, OrientationMatrix2D orientation)
	* Point2D src_point : source 2D points in image coordinates system
	* Point2D real_origin : origin in real world coordinates system
	* PixelSpacing real_spacing : distance between voxels in x and y direction
	* OrientationMatrix2D orientation : matrix for orientation
	* For given point in real world cordinates system, the function provides
	* the corresponding point in image coordinates system
	*/
	Point2D getImageCoordFromRealCoord2D(Point2D src_point, Point2D image_origin, PixelSpacing real_spacing, OrientationMatrix2D orientation);

	/*
	* This function return the interpolated value by nearest neighbor approach
	* dataType getInterpolatedValueNearestNeighbor2D(Image_Data2D src_image, Point2D point)
	* Image_Data2D src_image : source image
	* Point2D point : the current point
	*/
	dataType getInterpolatedValueNearestNeighbor2D(Image_Data2D src_image, Point2D point);

	/*
	* This function return the interpolated value by bilinear interpolation approach
	* dataType getInterpolatedValueNearestNeighbor2D(Image_Data2D src_image, Point2D point)
	* Image_Data2D src_image : source image
	* Point2D point : the current point
	*/
	dataType getInterpolatedValueBilinear2D(Image_Data2D src_image, Point2D point);

	/*
	* This function perform image interpolation
	* bool imageInterpolation2D(Image_Data2D src_image, Image_Data2D dest_image, dataType scale_factor, interpolationMethod method)
	* Image_Data2D src_image : source image structure
	*                          - height, width : image dimension
	*                          - imageDataPtr : pointer for pixels value
	*                          - origin : image origin
	*                          - spacing : pixel size
	* Image_Data2D dest_image : interpolated image structure
	*                          - height, width : interpolated image dimension
	*                          - imageDataPtr : pointer for interpolated pixels value
	*                          - origin : interpolated image origin
	*                          - spacing : interpolated pixel size
	* interpolationMethod method : nearest neighbor, bilinear, trilinear(for 3D images)
	*/
	bool imageInterpolation2D(Image_Data2D src_image, Image_Data2D dest_image, interpolationMethod method);

#ifdef __cplusplus
}

#endif