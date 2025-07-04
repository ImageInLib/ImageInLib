#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef GENERATE3DSHAPES_H
#define GENERATE3DSHAPES_H

	// INCLUDEs
#include <math.h>
#include <stdbool.h> // Boolean function bool
#include "common_functions.h" // Point3D struct
// Structures

// Macros

// Function Prototypes
/*
* Function to fill the points of Cuboid (block)
* inputHeight - inputDataArrayPtr height
* inputLength - inputDataArrayPtr length
* inputWidth - inputDataArrayPtr width
* Block Dimension to fill within - {Length, Width, Height}
* Fill value is set as inputDataArrayPtr[x][y][z] = fillValue
*/
	void fillBlock3D(dataType **inputDataArrayPtr, size_t inputHeight, size_t inputLength, size_t inputWidth, Point3D blockCorner, dataType *fillBlockDimension, dataType fillValue);
	/*
	* Function to fill 3D Ball Shape Points
	* inputHeight - inputDataArrayPtr height
	* inputLength - inputDataArrayPtr length
	* inputWidth - inputDataArrayPtr width
	*/
	void fillBall3D(dataType **inputDataArrayPtr, size_t inputHeight, size_t inputLength, size_t inputWidth, dataType sphereRadius, Point3D sphereCenter, dataType fillValue);

	bool generateSphereWithSixHoles(dataType ** dataArray3D, Point3D center, size_t length, size_t width, size_t height,
		dataType sphereRadius, dataType smallRadius, dataType fillValue, const char * outputPathPtr);

	bool generateSphere(dataType ** dataArray3D, Point3D center, size_t length, size_t width, size_t height,
		dataType sphereRadius, dataType fillValue, const char * outputPathPtr);

	bool generate3DSpiralTube(Image_Data image, Point3D center, double radius, double length, const size_t turns, const size_t segment_per_turn);

#endif // !GENERATE3DSHAPES_H

#ifdef __cplusplus
}
#endif