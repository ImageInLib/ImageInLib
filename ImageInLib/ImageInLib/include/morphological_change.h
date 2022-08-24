#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"
#include <stdbool.h>

	//Erosion ---> Shrink objects and remove boundaries pixels
	bool erosion2D(int** imageDataPtr, const size_t xDim, const size_t yDim, int background, int object);

	bool erosion3D(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background);

	bool dilatation3D(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background);


#ifdef __cplusplus
}
#endif
