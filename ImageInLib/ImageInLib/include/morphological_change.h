#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"
#include <stdbool.h>
#include "../src/filter_params.h"

	//Erosion removes floating pixels and thin lines so that only substantive objects remain
	bool erosion3D(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background);

	bool erosion3dHeighteenNeigbours(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background);

	//Dilation makes objects more visible and fills in small holes in objects.
	bool dilatation3D(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background);

	bool dilatation3dHeighteenNeigbours(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background);

#ifdef __cplusplus
}
#endif
