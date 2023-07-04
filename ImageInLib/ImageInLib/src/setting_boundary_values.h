#ifdef __cplusplus
extern "C" {
#endif

	/*
	* Author: Markjoe Olunna UBA
	* Purpose: ImageInLife project - 4D Image Segmentation Methods
	* Language:  C
	*/
#pragma once
#include <stdbool.h>
#include <stddef.h>
#include "common_functions.h"

	//t = step * tau
	bool setBoundaryExactValues3D(dataType **inputDataArrayPtr, size_t length, size_t width, size_t height, dataType sphereRadius, dataType t, dataType h);

	bool setBoundaryToZeroDirichletBC(dataType **inputDataArrayPtr, size_t length, size_t width, size_t height);

#ifdef __cplusplus
}
#endif