#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef FAST_MARCHING
#define FAST_MARCHING
	//==============================================================================
	/*
	* Header file to contain general functions, variables, structs
	* These can be used in different segmentation algorithms
	* Avoids redefinitions elsewhere
	*/
	//==============================================================================
	// Includes
#include "linked_list.h"
//==============================================================================
// Macro's
//==============================================================================
// STRUCTs
// Arrival time (s)
	typedef struct { dataType T; } Arrival_Time;
	//==============================================================================
	// ENUMs
	//==============================================================================
	// Global Variables
	//==============================================================================
	// Function Prototypes
	// Fast Marching 2D &3D
	void fastMarching3D(struct Node * band, Obj_Structure ** object, Point3D points[], Arrival_Time *known, size_t imageHeight, size_t imageLength, size_t imageWidth, size_t countPoints);
	//==============================================================================
	// Calls Fmm Method
	void fastMarching(dataType **distancePtr, dataType **dataSourcePtr, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType objPixel);
	//==============================================================================
#endif // !FAST_MARCHING

#ifdef __cplusplus
}
#endif