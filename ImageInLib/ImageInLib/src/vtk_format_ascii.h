#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdbool.h>

	//function for initialization of 3D array with some constant value.
	//xDim is the x dimension, yDim is the y dimension and zDim is the z dimension
	//value is the initial constant value
	bool vtkformat_3D(unsigned char ** array3DPtr, const size_t xDim, const size_t yDim,
		const size_t zDim, unsigned char * pathPtr);

	//function for initialization of 2D array with some constant value.
	bool vtkformat_2D(unsigned char * array2DPtr, const size_t yDim, size_t dim2D,
		unsigned char * pathPtr);

#ifdef __cplusplus
}
#endif