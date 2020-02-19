#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef SEGMENTATIO3D_COMMON
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
// Function to approximate value n + 1 U_pq, 1st order accurate, implicit in time, upwind method
	// Inputs - step n + 1 Up, step n + 1 Uq, step n Apq 
	// Eq - 82
	dataType approx_U_pq(dataType U_p, dataType U_q, dataType A_pq);
//==============================================================================
	// Function to calclate the gamma parameter - eq 91
	dataType chooseGamma(dataType value, dataType dist1, dataType dist2);
//==============================================================================
#endif // !SEGMENTATIO3D_COMMON

#ifdef __cplusplus
}
#endif