#ifdef __cplusplus
extern "C" {
#endif

#pragma once

#include<iostream>
#include "common_functions.h"
#include "../src/data_load.h"
#include "../src/endianity_bl.h"
#include <stdio.h>
#include <string.h>

	typedef struct {
		size_t x, y;
	} Point2D;

	dataType solve2dQuadratic(dataType X, dataType Y, dataType W);

	dataType selectX(dataType* distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t I, const size_t J);

	dataType selectY(dataType* distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t I, const size_t J);

	bool computeImageGradient(dataType * imageDataPtr, dataType * imageGradient, const size_t height, const size_t width, dataType h);

	dataType computeImageNorm2d(dataType* imageDataPtr, const size_t height, const size_t width);

	bool computePotential(dataType * imageDataPtr, dataType* potentialFuncPtr, const size_t height, const size_t width, Point2D* seedPoints);

	bool fastMarching2d(dataType* imageDataPtr, dataType* distanceFuncPtr, dataType* potentialFuncPtr, const size_t height, const size_t width, Point2D* seedPoints);

	bool shortestPath2d(dataType* distanceFuncPtr, dataType* resultedPath, const size_t height, const size_t width, dataType h, Point2D* seedPoints);

	bool findPath2d(dataType* distanceFuncPtr, dataType* resultedPath, const size_t height, const size_t width, Point2D* seedPoints);

	//bool bruteForce2d(dataType* imageDataPtr, dataType* distanceFuncPtr, const size_t height, const size_t width, dataType backGround);


#ifdef __cplusplus
}
#endif