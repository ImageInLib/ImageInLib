#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdio.h>
#include <stdbool.h>
#include "common_functions.h"
#include "transformation.h"
#include "interpolations.h"

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
	Point3D getImageCoordFromRealCoord3D(Point3D src_point, Point3D real_origin, VoxelSpacing real_spacing, OrientationMatrix orientation);

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
	Point3D getRealCoordFomImageCoord3D(Point3D src_point, Point3D real_origin, VoxelSpacing real_spacing, OrientationMatrix orientation);

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
	* This function get the neighbors of 2D given point
	* point : the given point
	* spacing : distance between pixels
	*/
	PointNeighbors2D getPointNeighbors2D(Point2D point, PixelSpacing spacing);

	/*
	* This function return the nearest neighbor of given point
	* Point2D getNearestNeighbor2D(Point2D point, PixelSpacing spacing)
	* point : the given point
	*/
	Point2D getNearestNeighbor2D(Point2D point, PixelSpacing spacing);

	/*
	* This function perform image interpolation from image coordinates system to real world coordinate system
	* bool interpolateImageFromImageCSToRealWorldCS(Image_Data2D src_image, Image_Data2D interp_image, dataType scale_factor)
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
	* dataType scale_factor : the scaling. if 0 < scale_factor < 1 ---> shrink
	*                                      if 1 < scale_factor < infty ---> magnify
	*/
	bool interpolateImageFromImageCStoRealWorldCS(Image_Data2D src_image, Image_Data2D interp_image, dataType scale_factor);

#ifdef __cplusplus
}

#endif