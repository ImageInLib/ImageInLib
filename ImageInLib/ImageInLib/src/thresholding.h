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

	/* image3DPtr is pointer to the 3D image, xDim, yDim and zDim are dimensions of x, y, and z respectively.
	* inputbgvalue is the value the voxel value (image3DPtr[k][i]) will be set to if the voxel value is
	less or equal to the threshold value
	* inputfgvalue is the value the voxel value (image3DPtr[k][i]) will be set to if the voxel value is
	greater than the threshold value
	*/
	bool thresholding3dFunctionUC(unsigned char ** image3DPtr, const size_t xDim, const size_t yDim,
		const size_t zDim, const unsigned char thresholdvalue, const unsigned char inputbgvalue,
		const unsigned char inputfgvalue);
	bool thresholding3dFunctionD(dataType ** image3DPtr, const size_t xDim, const size_t yDim,
		const size_t zDim, const dataType thresholdvalue, const dataType inputbgvalue,
		const dataType inputfgvalue);

#ifdef __cplusplus
}
#endif
