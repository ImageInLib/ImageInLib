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
	void fillBlock3D(double **inputDataArrayPtr, size_t inputHeight, size_t inputLength, size_t inputWidth, Point3D blockCorner, double *fillBlockDimension, double fillValue);
	/*
	* Function to fill 3D Ball Shape Points
	* inputHeight - inputDataArrayPtr height
	* inputLength - inputDataArrayPtr length
	* inputWidth - inputDataArrayPtr width
	*/
	void fillBall3D(double **inputDataArrayPtr, size_t inputHeight, size_t inputLength, size_t inputWidth, double sphereRadius, Point3D sphereCenter, double fillValue);

	bool generateSphereWithSixHoles(dataType ** dataArray3D, Point3D center, size_t length, size_t width, size_t height,
		dataType sphereRadius, dataType smallRadius, dataType fillValue, unsigned char * outputPathPtr);

	bool generateSphere(dataType ** dataArray3D, Point3D center, size_t length, size_t width, size_t height,
		dataType sphereRadius, dataType fillValue, unsigned char * outputPathPtr);

#endif // !GENERATE3DSHAPES_H

#ifdef __cplusplus
}
#endif