#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef COMMON_FUNCTIONS
#define COMMON_FUNCTIONS
	//==============================================================================
	//debug constants
#ifndef MEASURE_TIME
#define MEASURE_TIME
#endif
// Show/Hide Output in Console
#ifndef CONSOLE_OUTPUT
#define CONSOLE_OUTPUT
#endif
#undef CONSOLE_OUTPUT
//==============================================================================
// Typedefs
// Data Type

#ifndef USE_DOUBLE_PRECISION
	typedef float dataType;
#else
	typedef double dataType;
#endif // USE_DOUBLE_PRECISION

	//==============================================================================
	// Includes
#include <stddef.h>
//==============================================================================
// MACROs
//==============================================================================
// STRUCTs
// Common 3D Points - {x,y,z}
	typedef struct {
		dataType x, y, z;
	} Point3D;

	// Image Container and Properties
	typedef struct {
		// Image Dimensions
		size_t height, length, width; // Absolute Dimension

		dataType ** imageDataPtr; // Image Data Containers
	} Image_Data;
	//==============================================================================
	// Shapes Container
	typedef struct {
		// Shapes
		dataType **shp;
		int num; // shape number
	} Shapes;
	//==============================================================================
	// Enumeration
	enum FiniteDifference
	{
		FINITE_FORWARD = 1, FINITE_BACKWARD, FINITE_CENTRAL
	};
	//==============================================================================
	// FUNCTION PROTOTYPES
	/*
	* Calculate to 1D representation from 2D Array
	* x from column and row indices
	*/
	size_t x_new(const size_t rowIndex, const size_t columnIndex, const size_t rowLength);
	//==============================================================================
	/* Flattens 3D array to 1D array
	*/
	size_t x_flat(const size_t rowIndex, const size_t columnIndex, const size_t heightIndex, const size_t rowLength, const size_t columnLength);
	//==============================================================================
	/*
	* Reflection Function Zuska
	* toReflectImage is the data to perform reflection on
	* imageHeight is the actual Z dimension
	* imageLength is the actual X dimension
	* imageWidth is the actual Y dimension
	* P is the boundary constant/distance
	* Reflection does not include the current point but only its adjacent neighbor
	* Preserves Neumann Boundary condition - Target is for those methods
	*/
	void reflection3DB(dataType ** toReflectImage, size_t imageHeight, size_t imageLength, size_t imageWidth, size_t p);
	//==============================================================================
	/*
	* Reflection Function Jozef
	* Reflection includes the current point together with its adjacent neighbor
	* Targets Distance Function types
	*/
	void reflection3D(dataType ** toReflectImage, size_t imageHeight, size_t imageLength, size_t imageWidth);
	//==============================================================================
	/*
	* Gradient Calculation function
	*/
	dataType gradientFunction(dataType value, dataType coef);
	//==============================================================================
	/*
	void copyDataToExtendedArea(const dataType ** originalDataPtr, dataType ** extendedDataPtr, const size_t originalHeight, const size_t originalLength, const size_t originalWidth);

	The function copies data from source given by originalDataPtr
	representing data of dimensions given by originalHeight, originalLength and originalWidth
	to an array extendedDataPtr representing enlarged data (by 1 vx in each direction)
	of dimensions originalHeight + 2, originalLength + 2 and originalWidth + 2
	*/
	void copyDataToExtendedArea(const dataType ** originalDataPtr, dataType ** extendedDataPtr, const size_t originalHeight, const size_t originalLength, const size_t originalWidth);
	//==============================================================================
	void copyDataToReducedArea(dataType ** originalDataPtr, const dataType ** extendedDataPtr, const size_t originalHeight, const size_t originalLength, const size_t originalWidth);
	//==============================================================================
#endif // !COMMON_FUNCTIONS

#ifdef __cplusplus
}
#endif