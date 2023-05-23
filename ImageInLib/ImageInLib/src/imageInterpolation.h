#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdbool.h>
#include "common_functions.h"
#include "transformation.h"

	bool nearestNeighborInterpolation(dataType ** originalImage, dataType ** newImage, size_t imageLength, size_t imageWidth, size_t imageHeight,
		                  dataType originalSpacing, dataType newSpacing);

	bool linear2dInterpolation(dataType** originalImage, dataType** newImage, size_t imageLength, size_t imageWidth, size_t imageHeight,
						  dataType originalSpacing, dataType newSpacing);

#ifdef __cplusplus
}
#endif