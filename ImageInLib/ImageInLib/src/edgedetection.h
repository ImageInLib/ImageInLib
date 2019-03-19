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

	//edge detection for double
	/* image3DPtr is pointer to the 3D input image, xDim, yDim and zDim are dimensions of x, y, and z respectively.
	* edge3DPtr is pointer to the 3D output image containing edge values.
	* fgroundvalue is foreground input value
	* bgroundvalue is background input value
	*/
	bool edgeDetection3dFunctionD(dataType ** image3DPtr, dataType ** edge3DPtr, const size_t xDim,
		const size_t yDim, const size_t zDim, const dataType bgroundvalue, const dataType fgroundvalue);

	/* image3DPtr is pointer to the 3D input image, xDim, yDim and zDim are dimensions of x, y, and z respectively.
	* edge3DPtr is pointer to the 3D output image containing edge values.
	* fgroundvalue is foreground input value
	* bgroundvalue is background input value
	*/
	bool edgeDetection3dFunctionUC(unsigned char ** image3DPtr, unsigned char ** edge3DPtr, const size_t xDim,
		const size_t yDim, const size_t zDim, const unsigned char bgroundvalue, const unsigned char fgroundvalue);
	bool edgeDetection2dFunctionUC(unsigned char * image2DPtr, unsigned char * edge2DPtr, const size_t xDim,
		const size_t yDim, const unsigned char bgroundvalue, const unsigned char fgroundvalue);

#ifdef __cplusplus
}
#endif