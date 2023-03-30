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

	typedef struct {
		size_t x, y;
	} Point2D;

	typedef struct {
		size_t x, y, z;
		dataType arrival;
	}neighborPoint;

	//class neighborPoint {
	//private:
	//	size_t x, y, z;
	//	dataType arrival;
	//	size_t label;
	//	size_t position;
	//public:
	//	neighborPoint(size_t i, size_t j, size_t k, size_t dimI, size_t  dimJ) {
	//		x = j;
	//		y = i;
	//		z = k;
	//		arrival = INFINITY;
	//		label = 3;
	//		position = x_flat(j, i, k, dimJ, dimJ);
	//	}
	//	size_t getJ() {
	//		return x;
	//	}
	//	size_t getI() {
	//		return y;
	//	}
	//	size_t getZ() {
	//		return z;
	//	}
	//	size_t getLabel() {
	//		return label;
	//	}
	//	size_t getPosition() {
	//		return position;
	//	}
	//	dataType getArrival() {
	//		return arrival;
	//	}
	//	void setArrival(dataType a) {
	//		arrival = a;
	//	}
	//	void setLabel(size_t l) {
	//		label = l;
	//	}
	//};

	dataType solve2dQuadratic(dataType X, dataType Y, dataType W);

	dataType selectX(dataType* distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t I, const size_t J);

	dataType selectY(dataType* distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t I, const size_t J);

	bool computeImageGradient(dataType* imageDataPtr, dataType* gradientVectorX, dataType* gradientVectorY, const size_t height, const size_t width, dataType h);

	dataType computeGradientNorm2d(dataType* gradientVectorX, dataType* gradientVectorY, const size_t height, const size_t width);

	bool computePotential(dataType * imageDataPtr, dataType* potentialFuncPtr, const size_t height, const size_t width, Point2D* seedPoints);

	bool fastMarching2d(dataType* imageDataPtr, dataType* distanceFuncPtr, dataType* potentialFuncPtr, const size_t height, const size_t width, Point2D* seedPoints);

	bool shortestPath2d(dataType* distanceFuncPtr, dataType* resultedPath, const size_t height, const size_t width, dataType h, Point2D* seedPoints);

	//==============================================================
	 
	//3D functions
	typedef struct {
		size_t x, y, z;
	} Point3d;

	dataType solve3dQuadratic(dataType X, dataType Y, dataType Z, dataType W);

	dataType select3dX(dataType** distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t dimK, const size_t I, const size_t J, const size_t K);

	dataType select3dY(dataType** distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t dimK, const size_t I, const size_t J, const size_t K);

	dataType select3dZ(dataType** distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t dimK, const size_t I, const size_t J, const size_t K);

	dataType computeGradientNorm3d(dataType** gradientVectorX, dataType** gradientVectorY, dataType** gradientVectorZ, const size_t length, const size_t width, const size_t height);

	bool compute3dImageGradient(dataType** imageDataPtr, dataType** gradientVectorX, dataType** gradientVectorY, dataType** gradientVectorZ, const size_t lenght, const size_t width, const size_t height, dataType h);

	bool compute3dPotential(dataType** imageDataPtr, dataType** potentialFuncPtr, const size_t length, const size_t width, const size_t height, Point3d* seedPoints);

	bool fastMarching3d_N(dataType** imageDataPtr, dataType** distanceFuncPtr, dataType** potentialFuncPtr, const size_t length, const size_t width, const size_t height, Point3d* seedPoints);

	bool shortestPath3d(dataType** distanceFuncPtr, dataType** resultedPath, const size_t length, const size_t width, const size_t height, dataType h, Point3d* seedPoints);

	bool FastMarching3DNewVersion(dataType** imageDataPtr, dataType** distanceFuncPtr, dataType** potentialFuncPtr, const size_t length, const size_t width, const size_t height, Point3d* seedPoints);

	void swapNeighbor(neighborPoint* a, neighborPoint* b);

	void heapify(vector<neighborPoint> &in_Process, int length_InProcess, int i);

	void createMinHeapStructure(vector<neighborPoint> &in_Process, int length_InProcess);

	void heapifyBottomToUp(vector<neighborPoint> &in_Process, int i);

#ifdef __cplusplus
}
#endif