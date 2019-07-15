#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H
	//==============================================================================
	/*
	* Header file to contain general functions, variables, structs
	* These can be used in transformation algorithms
	* Avoids redefinitions elsewhere
	*/
	//==============================================================================
	// INCLUDES
#include "common_math.h" // PI value
#include "interpolations.h" // Interpolation Function
#include <stdbool.h>
// MACROs
//==============================================================================
// STRUCTS
//==============================================================================
	/*
* Rotated indices
*/
	void coordinate_rotate(dataType z, dataType x, dataType y, dataType theta, dataType psi, dataType phi, dataType* k_t, dataType* i_t, dataType* j_t);
	dataType x_rotate(dataType z, dataType x, dataType y, dataType theta, dataType psi);
	dataType y_rotate(dataType z, dataType x, dataType y, dataType theta, dataType psi, dataType phi);
	dataType z_rotate(dataType z, dataType x, dataType y, dataType theta, dataType psi, dataType phi);
	//==============================================================================
	// Inverse
	dataType x_rotateInv(dataType z, dataType x, dataType y, dataType theta, dataType psi, dataType phi);
	dataType y_rotateInv(dataType z, dataType x, dataType y, dataType theta, dataType psi, dataType phi);
	dataType z_rotateInv(dataType z, dataType x, dataType y, dataType theta, dataType psi, dataType phi);
	//==============================================================================
//==============================================================================
// FUNCTION PROTOTYPES
/*
* Transform Function for imageDataPtr
* Point3D is structure for 3 3D points - x, y,z
* translation - tx,ty,tz
* scaling - sx, sy, sz
* rotation - theta, psi, phi - Are angles in degrees!
* imageHeight is the actual Z dimension
* imageLength is the actual X dimension
* imageWidth is the actual Y dimension
* bgValue is the background value to set if the point is not transformed
*/
	void transform3DImage(dataType ** sourceDataPtr, dataType **transformPointsPtr, Point3D translation, Point3D scaling, Point3D rotation, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType bgValue, dataType centroid[3], dataType imageForeground, bool parallelize);
	void transformInverse3DImage(dataType **sourceDataPtr, dataType **imageDataPtr, Point3D translation, Point3D scaling, Point3D rotation, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType bgValue, dataType centroid[3]);
	//==============================================================================
#endif // !TRANSFORMATION_H

#ifdef __cplusplus
}
#endif