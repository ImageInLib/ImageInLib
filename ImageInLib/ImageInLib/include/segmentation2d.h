#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

#include "common_functions.h"
#include "../src/heat_equation.h"
#include "../src/filter_params.h"
#include "../distanceForPathFinding.h"

	dataType min(dataType a, dataType b);

	bool rescaleToZeroOne2d(dataType* imageDataPtr, const size_t height, const size_t width);

	bool generateInitialSegmentationFunction(dataType* imageDataPtr, const size_t height, const size_t width, point2D* center, dataType v, dataType R);

#ifdef __cplusplus
}
#endif