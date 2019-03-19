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

	double norm3dDataArrayD(dataType ** dataArray3DPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType h);

	bool l2norm3dDataArrayD(dataType ** dataArray3DPtr1, dataType ** dataArray3DPtr2, unsigned char * pathPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType h,
		size_t step);

	dataType timespacel2norm3dDataArrayD(dataType ** dataArray3DPtr1, dataType ** dataArray3DPtr2, const size_t xDim, const size_t yDim, const size_t zDim, dataType h);

	dataType l2normD(double ** dataArray3DPtr1, double ** dataArray3DPtr2, const size_t xDim, const size_t yDim, const size_t zDim, double h);

#ifdef __cplusplus
}
#endif
