#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef SEGMENTATION3D_COMMON
#define SEGMENTATION3D_COMMON
	//==============================================================================
	/*
	* Header file to contain general functions, variables, structs
	* These can be used in different segmentation algorithms
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
	// Function to calclate the gamma parameter - eq 91
	dataType chooseGamma(dataType value, dataType dist1, dataType dist2);
//==============================================================================
#endif // !SEGMENTATION3D_COMMON

#ifdef __cplusplus
}
#endif