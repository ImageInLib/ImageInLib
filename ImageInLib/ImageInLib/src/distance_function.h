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
	/*bruteForceFunction_3D is a function that implements brute force algorithm in 3D. distance3DPtr is a
	pointer to the 3D array that will hold the so called distance function, curve3DPtr is a pointer to the 3D
	array that will hold the input boundary image, xDim, yDim, zDim are the x, y, and z dimensions respectively,
	largeValue is a very relatively large number which according to the algorithm should be set to infinity,
	fgroundValue is the value at the object which is different from background value*/
	bool bruteForceFunction_3D(dataType ** distance3DPtr, dataType ** curve3DPtr,
		const size_t xDim, const size_t yDim, const size_t zDim, const dataType largeValue,
		const dataType fgroundValue);

	/*rouyTourinFunction_3D is a function that implements Rouy - Tourin scheme in 3D. distance3DPtr is a
	pointer to the 3D array that will hold the so called distance function, image3DPtr is a pointer to the 3D
	array that will hold the original input image, xDim, yDim, zDim are the x, y, and z dimensions respectively,
	tau is a scale value and h is the lenght of the voxel*/
	bool rouyTourinFunction_3D(dataType ** distance3DPtr, dataType ** image3DPtr, dataType tolerance,
		const size_t xDim, const size_t yDim, const size_t zDim, dataType tau, const dataType h);

	bool fastSweepingFunction_3D(dataType ** distance3DPtr, dataType ** curve3DPtr, const size_t xDim, const size_t yDim, const size_t zDim,
		const dataType h, const dataType largeValue, const dataType fgroundValue);

	bool fastSweepingFunction_2D(dataType * distance2DPtr, dataType * curve2DPtr, const size_t xDim, const size_t yDim,
		const dataType h, const dataType largeValue, const dataType fgroundValue);

	dataType compute2dDistance(size_t i_n, const size_t xDim, const size_t yDim, dataType *distance2DPtr, dataType *temp2dPtr, const size_t dim2D,
		const dataType h);

	dataType compute3dDistance(size_t i_n, size_t k, const size_t xDim, const size_t yDim, const size_t zDim, dataType **distance3DPtr,
		dataType **temp3dPtr, const size_t dim2D, const dataType h);

#ifdef __cplusplus
}
#endif