//#ifdef __cplusplus
//extern "C" {
//#endif

#pragma once
#include<iostream>
#include <stdio.h>
#include <string>

#include "common_functions.h"
#include "../src/heat_equation.h"
#include "../src/filter_params.h"
#include "../distanceForPathFinding.h"

	typedef struct
	{
		dataType* East;
		dataType* West;
		dataType* North;
		dataType* South;
	} neighPtrs;

	dataType min(dataType a, dataType b);

	bool rescaleToZeroOne2d(dataType* imageDataPtr, const size_t height, const size_t width);

	bool generateInitialSegmentationFunction(dataType* imageDataPtr, const size_t height, const size_t width, point2D* center, dataType v, dataType R);

	bool set2dDirichletBoundaryCondition(dataType* imageDataPtr, const size_t height, const size_t width);

	bool computeNormOfGradientDiamondCells(dataType* arrayPtr, neighPtrs neigbours, const size_t height, const size_t width, dataType h);

	bool epsilonRegularization(neighPtrs neighbours, const size_t height, const size_t width, dataType epsilon);

//#ifdef __cplusplus
//}
//#endif