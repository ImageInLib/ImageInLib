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
// Function Prototypes
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
