#ifdef __cplusplus
extern "C" {
#endif
	//==============================================================================
#pragma once
#ifndef MHD
#define MHD
	//==============================================================================
	/*
	* Header file to contain general functions, variables, structs
	* These can be used to calculate the mean hausdorff distance
	* Avoids redefinitions elsewhere
	*/
	//==============================================================================
	// Includes
#include "common_functions.h"
	//==============================================================================
	dataType mean_hausdorff(dataType ** curveA_Pointer, dataType curveB_Pointer, dataType curveIntensity, size_t height, size_t length, size_t width);
	//==============================================================================
#endif // !MHD


#ifdef __cplusplus
}
#endif