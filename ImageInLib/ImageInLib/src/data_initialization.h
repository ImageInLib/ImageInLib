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
#include "common_functions.h"

	bool initialize3dArrayUC(unsigned char ** array3DPtr, const size_t xDim, const size_t yDim,
		const size_t zDim, const unsigned char value);
	bool initialize2dArrayUC(unsigned char * array2DPtr, const size_t xDim, const size_t yDim, const unsigned char value);

	bool initialize3dArrayD(dataType ** array3DPtr, const size_t xDim, const size_t yDim,
		const size_t zDim, dataType value);
	bool initialize2dArrayD(dataType * array2DPtr, const size_t xDim, const size_t yDim, dataType value);

#ifdef __cplusplus
}
#endif