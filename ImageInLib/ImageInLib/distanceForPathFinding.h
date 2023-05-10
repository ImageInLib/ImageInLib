#ifdef __cplusplus
extern "C" {
#endif

#pragma once

#include<iostream>
#include<vector>
#include "common_functions.h"
#include "../src/data_load.h"
#include "../src/endianity_bl.h"
#include <stdio.h>
#include <string.h>

using namespace std;

	//2D functions

	typedef struct {
		size_t x, y;
	} point2D;

	typedef struct {
		size_t x, y;
		dataType arrival;
	}pointFastMarching2D;

	// heap functions

	void swap2dPoints(pointFastMarching2D* a, pointFastMarching2D* b);

	void heapifyDown2D(vector<pointFastMarching2D>& in_Process, int pos);

	void heapifyUp2D(vector<pointFastMarching2D>& in_Process, int pos);

	void heapifyVector2D(vector<pointFastMarching2D>& in_Process);

	void deleteRootHeap2D(vector<pointFastMarching2D>& in_Process);

	void addPointHeap2D(vector<pointFastMarching2D>& in_Process, pointFastMarching2D point);

	bool fastMarching2D(dataType* imageDataPtr, dataType* distancePtr, dataType* potentialPtr, const size_t height, const size_t width, point2D* seedPoints);

	//fast marching functions

	dataType solve2dQuadratic(dataType X, dataType Y, dataType W);

	dataType selectX(dataType* distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t I, const size_t J);

	dataType selectY(dataType* distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t I, const size_t J);

	bool computeImageGradient(dataType* imageDataPtr, dataType* gradientVectorX, dataType* gradientVectorY, const size_t height, const size_t width, dataType h);

	dataType computeGradientNorm2d(dataType* gradientVectorX, dataType* gradientVectorY, const size_t height, const size_t width);

	bool computePotential(dataType * imageDataPtr, dataType* potentialFuncPtr, const size_t height, const size_t width, point2D* seedPoints);

	bool shortestPath2d(dataType* distanceFuncPtr, dataType* resultedPath, const size_t height, const size_t width, dataType h, point2D* seedPoints);

	//==============================================================
	 
	//3D functions
	typedef struct {
		size_t x, y, z;
	} point3d;

	typedef struct {
		size_t x, y, z;
		dataType arrival;
	}pointFastMarching3D;

	dataType solve3dQuadratic(dataType X, dataType Y, dataType Z, dataType W);

	dataType select3dX(dataType** distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t dimK, const size_t I, const size_t J, const size_t K);

	dataType select3dY(dataType** distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t dimK, const size_t I, const size_t J, const size_t K);

	dataType select3dZ(dataType** distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t dimK, const size_t I, const size_t J, const size_t K);

	dataType computeGradientNorm3d(dataType** gradientVectorX, dataType** gradientVectorY, dataType** gradientVectorZ, const size_t length, const size_t width, const size_t height);

	bool compute3dImageGradient(dataType** imageDataPtr, dataType** gradientVectorX, dataType** gradientVectorY, dataType** gradientVectorZ, const size_t lenght, const size_t width, const size_t height, dataType h);

	bool compute3dPotential(dataType** imageDataPtr, dataType** potentialFuncPtr, const size_t length, const size_t width, const size_t height, point3d* seedPoints);

	//heap functions
	void swap3dPoints(pointFastMarching3D* a, pointFastMarching3D* b);

	void heapifyDown3D(vector<pointFastMarching3D>& in_Process, int i);

	void heapifyVector3D(vector<pointFastMarching3D>& in_Process);

	void heapifyUp3D(vector<pointFastMarching3D>& in_Process, int i);

	void deleteRootHeap3D(vector<pointFastMarching3D>& in_Process);

	void addPointHeap3D(vector<pointFastMarching3D>& in_Process, pointFastMarching3D point);

	bool fastMarching3D_N(dataType** imageDataPtr, dataType** distanceFuncPtr, dataType** potentialFuncPtr, const size_t length, const size_t width, const size_t height, point3d* seedPoints);

	bool shortestPath3d(dataType** distanceFuncPtr, dataType** resultedPath, const size_t length, const size_t width, const size_t height, dataType h, point3d* seedPoints);

	//================================================================

#ifdef __cplusplus
}
#endif