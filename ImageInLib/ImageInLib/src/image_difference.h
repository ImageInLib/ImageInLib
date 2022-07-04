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

	/* origImage3DPtr is pointer to the original 3D image, xDim, yDim and zDim are dimensions of x, y, and z respectively.
	* processedImage3DPtr is pointer to the processed 3D image, xDim, yDim and zDim are dimensions of x, y, and z respectively.
	* outputImage3DPtr is pointer to the output 3D image, xDim, yDim and zDim are dimensions of x, y, and z respectively.
	*/
	bool imageDifference3dDataArrayD(dataType ** origImage3DPtr, dataType ** processedImage3DPtr, dataType ** outputImage3DPtr,
		const size_t xDim, const size_t yDim, const size_t zDim);

#ifdef __cplusplus
}
#endif
