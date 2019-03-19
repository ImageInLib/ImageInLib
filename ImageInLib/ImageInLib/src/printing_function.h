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

	bool printing_function3D(unsigned char **array3DPtr, const size_t xDim, const size_t yDim, const size_t zDim);
	bool printing_function2D(unsigned char *array2DPtr, const size_t xDim, const size_t yDim);

#ifdef __cplusplus
}
#endif