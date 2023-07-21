#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef INTERPOLATION
#define INTERPOLATION
	//==============================================================================
	/*
	* Header file to contain functions, variables, structs
	* These can be used in different implementing interpolation algorithms
	* Avoids redefinitions elsewhere
	*/
	//==============================================================================
	// Includes
#include "common_functions.h"
//==============================================================================
// Macro's
//==============================================================================
// STRUCTs
//==============================================================================
// ENUMs
//==============================================================================
// Global Variables
//==============================================================================

// Local Function Prototypes
	/*
	* Function to calculate linear values R1, R2 inputs for Bilinear
	* 1D
	*/
	dataType linearInterpolation(dataType x, dataType x1, dataType x2, dataType q00, dataType q01);

	/*
	* Function to calculate Bilinear Interpolation
	* 2D
	* Assume a Rectangular shaper
	* Out rectangle
	* x - bottom center, x1 bottom left , x2 -bottom right
	* y - vertical center, y1 vertical bottom, y2 vertical top
	* q00, q01 slopes
	* P = (y,x) - Point to be interpolated/Returned
	* Inner rectangle
	* Intersections of above points
	* q11 = (y1,x1) - Bottom left, q12 = (y1,x2) - Bottom right
	* q21 = (y2,x1) - Top left, q22 = (y2,x2) - Top Right
	* R1 = (y1,x) - Between q11 and q21
	* R2 = (y1,x) - Between q12 and q22
	* P = (y,x) - Point to be interpolated/Returned. Between R1 and R2
	*
	*/
	dataType bilinearInterpolation(dataType x, dataType x1, dataType x2, dataType q11, dataType q12, dataType q21, dataType q22, dataType y, dataType y1, dataType y2);

	/*
	* This function performs Tr-Linear interpolation
	* x, y, are the point coordinates where we want to find the function at
	* x1, x2, y1, y2 are the 4 additional points used in rectilinear 2D grid surrounding the interpolation point:bilinear interpolation
	* Cxx are the Cuboidd 8 corner points surrounding the interpolation point
	*/
	dataType trilinearInterpolation(dataType x, dataType x1, dataType x2, dataType y, dataType y1, dataType y2, dataType c000, dataType c001, dataType c010, dataType c011, dataType c100, dataType c101, dataType c110, dataType c111, dataType z, dataType z1, dataType z2);
	//==============================================================================
#endif // !INTERPOLATION

#ifdef __cplusplus
}
#endif
