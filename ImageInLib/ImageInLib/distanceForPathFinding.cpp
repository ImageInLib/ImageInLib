/*
* Author: Konan ALLALY
* Purpose: INFLANET project - Image Processing in Nuclear Medicine (2D/3D)
* Language:  C and C++
*/
#include <iostream>
#include <climits>
#include <crtdbg.h>
#include <corecrt_malloc.h>
#include<cmath>
#include <omp.h>
#include<vector>
#include<tuple>
#include "distanceForPathFinding.h"
#include<template_functions.h>

#define BIG_VALUE INFINITY

using namespace std;

//J.A Sethian, A Fast Marching Level Set method for Monotonically advancing fronts, 1995, page 8 and 10.
//link to article ---> http://ugweb.cs.ualberta.ca/~vis/courses/CompVis/readings/modelrec/sethian95fastlev.pdf

//Functions for 2D images

dataType solve2dQuadratic(dataType X, dataType Y, dataType W) {
	//This fuction is used the solve the following quadratic coming the discretization 
	// by upwind principle in the implementation of the fast marching method 
	// aU^2 -2U(X+Y) + (X^2 + Y^2 - W) = 0
	dataType sol = 0.0, a, b, c, delta;

	a = 2.0; 
	if (X == INFINITY) {
		X = 0; a--;
	}
	if (Y == INFINITY) {
		Y = 0; a--;
	}

	b = -2 * (X + Y); c = pow(X, 2) + pow(Y, 2) - W;
	delta = pow(b, 2) - 4 * a * c;

	if (delta >= 0) {
		sol = (-b + sqrt(delta)) / (2 * a);
	}
	else {
		sol = min(X, Y) + W;
	}

	if (sol < 0) {
		cout << "The solution is negative " << endl;
		return 0;
	}
	else {
		return sol;
	}
	
}

dataType selectX(dataType* distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t I, const size_t J) {
	// this function return the minimum in the upwind principle
	// x--->j and y--->i
	dataType j_minus, j_plus;

	if (J == 0) {
		j_minus = INFINITY;
	}
	else {
		j_minus = distanceFuncPtr[x_new(J - 1, I, dimJ)];
	}

	if (J == dimJ - 1) {
		j_plus = INFINITY;
	}
	else {
		j_plus = distanceFuncPtr[x_new(J + 1, I, dimJ)];
	}

	return min(j_minus, j_plus);
}

dataType selectY(dataType* distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t I, const size_t J) {
	// this function return the minimum in the upwind principle
	// x--->j and y--->i
	dataType i_minus, i_plus;

	if (I == 0) {
		i_minus = BIG_VALUE;
	}
	else {
		i_minus = distanceFuncPtr[x_new(J, I - 1, dimJ)];
	}

	if (I == dimI - 1) {
		i_plus = BIG_VALUE;
	}
	else {
		i_plus = distanceFuncPtr[x_new(J, I + 1, dimJ)];
	}

	return min(i_minus, i_plus);
}

bool computeImageGradient(dataType* imageDataPtr, dataType* gradientVectorX, dataType* gradientVectorY, const size_t height, const size_t width, dataType h) {
	//This function compute the gradient of image by finite difference
	// x--->j and y--->i

	if (imageDataPtr == NULL || gradientVectorX == NULL || gradientVectorY == NULL) {
		return false;
	}

	size_t i, j, xd;
	dataType ux = 0.0, uy = 0.0;

	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			xd = x_new(j, i, width);
			if (j == 0) {
				ux = (imageDataPtr[x_new(j + 1, i, width)] - imageDataPtr[x_new(j, i, width)]) / h;
			}
			else {
				if (j == width - 1) {
					ux = (imageDataPtr[x_new(j, i, width)] - imageDataPtr[x_new(j - 1, i, width)]) / h;
				}
				else {
					ux = (imageDataPtr[x_new(j + 1, i, width)] - imageDataPtr[x_new(j - 1, i, width)]) / (2 * h);
				}
			}
			if (i == 0) {
				uy = (imageDataPtr[x_new(j, i + 1, width)] - imageDataPtr[x_new(j, i, width)]) / h;
			}
			else {
				if (i == height - 1) {
					uy = (imageDataPtr[x_new(j, i, width)] - imageDataPtr[x_new(j, i - 1, width)]) / h;
				}
				else {
					uy = (imageDataPtr[x_new(j, i + 1, width)] - imageDataPtr[x_new(j, i - 1, width)]) / (2 * h);
				}
			}
			gradientVectorX[xd] = ux;
			gradientVectorY[xd] = uy;
		}
	}

	return true;
}

dataType computeGradientNorm2d(dataType* gradientVectorX, dataType* gradientVectorY, const size_t height, const size_t width) {
	//this function is used to compute the norm of the gradient
	size_t i, j, xd;
	dataType norm_array = 0.0;

	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			xd = x_new(j, i, width);
			norm_array = norm_array + pow(gradientVectorX[xd], 2) + pow(gradientVectorY[xd], 2);
		}
	}
	return sqrt(norm_array);
}

bool computePotential(dataType* imageDataPtr, dataType* potentialFuncPtr, const size_t height, const size_t width, point2D * seedPoints)
{
	//This function is used to compute the potential function
	//were epsilon can be used also as parameter
	//epislon = 0.01 was perfect for our experimentations

	if (imageDataPtr == NULL || potentialFuncPtr == NULL || seedPoints == NULL)
		return false;

	size_t i, j;
	size_t i1 = seedPoints[0].y, j1 = seedPoints[0].x;
	size_t i2 = seedPoints[1].y, j2 = seedPoints[1].x;

	dataType seedVal = (imageDataPtr[x_new(j1, i1, width)] + imageDataPtr[x_new(j2, i2, width)]) / 2;
	size_t currentIndx = 0;
	dataType epsylon = 0.01;
	dataType K = 0.00005;

	dataType* gradientVectorX = new dataType[height * width];
	dataType* gradientVectorY = new dataType[height * width];
	computeImageGradient(imageDataPtr, gradientVectorX, gradientVectorY, height, width, 1.0);

	size_t seedIndice = x_new(j1, i1, width);
	dataType ux0 = gradientVectorX[seedIndice], uy0 = gradientVectorY[seedIndice], ux = 0.0, uy = 0.0;

	//Computation of potential function
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			currentIndx = x_new(j, i, width);
			potentialFuncPtr[currentIndx] = abs(seedVal - imageDataPtr[currentIndx]);
		}
	}

	//find the max potential
	dataType max_potential = -1 * INFINITY;
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			currentIndx = x_new(j, i, width);
			if (potentialFuncPtr[currentIndx] > max_potential) {
				max_potential = potentialFuncPtr[currentIndx];
			}
		}
	}

	//Normalization
	dataType weight = 0.0, edgeValue = 0.0, norm_of_gradient = 0.0;
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			currentIndx = x_new(j, i, width);
			ux = gradientVectorX[currentIndx]; uy = gradientVectorY[currentIndx]; norm_of_gradient = sqrt(ux * ux + uy * uy);
			edgeValue = 1 + K * (ux * ux + uy * uy);
			weight = potentialFuncPtr[currentIndx] / max_potential;
			potentialFuncPtr[currentIndx] = (dataType)(epsylon + weight * edgeValue);
		}
	}

	delete[] gradientVectorX;
	delete[] gradientVectorY;

	return true;
}

//function to swap 2D points in the fast marching contest
void swap2dPoints(pointFastMarching2D* a, pointFastMarching2D* b) {
	pointFastMarching2D temp = *a;
	*a = *b;
	*b = temp;
}

//pos is the index of the element to be heapified
void heapifyDown2D(vector<pointFastMarching2D>& in_Process, int pos) {
	//we use type int for indexes because we do operations like pos--
	int length_array = in_Process.size();
	int current = pos;
	int left_child = 2 * pos + 1;
	int right_child = 2 * pos + 2;

	dataType val_current = 0.0, val_left = 0.0, val_right = 0.0;

	if (current >= 0 && current < length_array) {
		val_current = in_Process[current].arrival;
	}

	if (left_child < length_array) {
		val_left = in_Process[left_child].arrival;
		if (val_left < val_current) {
			current = left_child;
			val_current = in_Process[current].arrival;
		}
	}

	if (right_child < length_array) {
		val_right = in_Process[right_child].arrival;
		if (val_right < val_current) {
			current = right_child;
		}
	}

	if (current != pos) {
		swap2dPoints(&in_Process[pos], &in_Process[current]);
		heapifyDown2D(in_Process, current);
	}
}

void heapifyUp2D(vector<pointFastMarching2D>& in_Process, int i) {

	int current = i;

	if (i > 0) {
		int parent = (i - 1) / 2;
		dataType val_current = in_Process[current].arrival;
		dataType val_parent = in_Process[parent].arrival;
		if (val_current < val_parent) {
			current = parent;
		}
	}

	if (current != i) {
		swap2dPoints(&in_Process[current], &in_Process[i]);
		heapifyUp2D(in_Process, current);
	}

}

void heapifyVector2D(vector<pointFastMarching2D>& in_Process) {
	//we use type int for indexes because we do operations like pos--
	int length_array = in_Process.size();
	int indx, start = length_array / 2 - 1;
	for (indx = start; indx >= 0; indx--) {
		heapifyDown2D(in_Process, indx);
	}
}

void deleteRootHeap2D(vector<pointFastMarching2D>& in_Process) {
	//we use type int for indexes because we do operations like pos--
	int l = in_Process.size();
	swap2dPoints(&in_Process[0], &in_Process[l - 1]);
	in_Process.pop_back();
	heapifyDown2D(in_Process, 0);
}

void addPointHeap2D(vector<pointFastMarching2D>& in_Process, pointFastMarching2D point) {
	//we use type int for indexes because we do operations like pos--
	in_Process.push_back(point);
	int l = in_Process.size();
	heapifyUp2D(in_Process, l - 1);
}

bool fastMarching2D(dataType* imageDataPtr, dataType* distancePtr, dataType* potentialPtr, const size_t height, const size_t width, point2D* seedPoints) {

	short* labelArray = new short[height * width];

	if (imageDataPtr == NULL || distancePtr == NULL || potentialPtr == NULL || seedPoints == NULL || labelArray == NULL) {
		return false;
	}

	size_t i = 0, j = 0, k = 0, dim2D = height * width;
	vector<pointFastMarching2D> inProcess;

	dataType x = 0.0, y = 0.0, coefSpeed = 0.0;

	//Compute the potential function
	computePotential(imageDataPtr, potentialPtr, height, width, seedPoints);

	//STEP 1
	//In labelAray we have : 1 ---> already processed, 2 ---> in process and 3 ---> not processed
	for (k = 0; k < dim2D; k++) {
		distancePtr[k] = INFINITY;
		labelArray[k] = 3;
	}

	i = seedPoints[0].y;
	j = seedPoints[0].x;
	size_t currentIndx = x_new(j, i, width);
	distancePtr[currentIndx] = 0.0;
	labelArray[currentIndx] = 1;

	//East
	if (j < width - 1 && i >= 0 && i < height) {
		size_t jplus = j + 1;
		x = selectX(distancePtr, height, width, i, jplus);
		y = selectY(distancePtr, height, width, i, jplus);
		size_t indxEast = x_new(jplus, i, width);
		coefSpeed = potentialPtr[indxEast];
		dataType dEast = solve2dQuadratic(x, y, coefSpeed);
		pointFastMarching2D EastNeighbor = { jplus, i, dEast };
		distancePtr[indxEast] = dEast;
		inProcess.push_back(EastNeighbor);
		labelArray[indxEast] = 2;
	}

	//West
	if (j > 0 && i >= 0 && i < height) {
		size_t jminus = j - 1;
		x = selectX(distancePtr, height, width, i, jminus);
		y = selectY(distancePtr, height, width, i, jminus);
		size_t indxWest = x_new(jminus, i, width);
		coefSpeed = potentialPtr[indxWest];
		dataType dWest = solve2dQuadratic(x, y, coefSpeed);
		pointFastMarching2D WestNeighbor = { jminus, i, dWest };
		distancePtr[indxWest] = dWest;
		inProcess.push_back(WestNeighbor);
		labelArray[indxWest] = 2;
	}

	//North
	if (j >= 0 && j < width && i > 0) {
		size_t iminus = i - 1;
		x = selectX(distancePtr, height, width, iminus, j);
		y = selectY(distancePtr, height, width, iminus, j);
		size_t indxNorth = x_new(j, iminus, width);
		coefSpeed = potentialPtr[indxNorth];
		dataType dNorth = solve2dQuadratic(x, y, coefSpeed);
		pointFastMarching2D NorthNeighbor = { j, iminus, dNorth };
		distancePtr[indxNorth] = dNorth;
		inProcess.push_back(NorthNeighbor);
		labelArray[indxNorth] = 2;
	}

	//South
	if (j >= 0 && j < width && i < height - 1) {
		size_t iplus = i + 1;
		x = selectX(distancePtr, height, width, iplus, j);
		y = selectY(distancePtr, height, width, iplus, j);
		size_t indxSouth = x_new(j, iplus, width);
		coefSpeed = potentialPtr[indxSouth];
		dataType dSouth = solve2dQuadratic(x, y, coefSpeed);
		pointFastMarching2D SouthNeighbor = { j, iplus, dSouth };
		distancePtr[indxSouth] = dSouth;
		inProcess.push_back(SouthNeighbor);
		labelArray[indxSouth] = 2;
	}

	//heapify 2D vector
	heapifyVector2D(inProcess);

	pointFastMarching2D current;
	short label = 0;

	while (inProcess.size() != 0) {

		current = inProcess[0];
		j = current.x;
		i = current.y;
		currentIndx = x_new(j, i, width);
		labelArray[currentIndx] = 1;
		distancePtr[currentIndx] = current.arrival;
		deleteRootHeap2D(inProcess);

		//East
		if (j < width - 1 && i >= 0 && i < height) {
			size_t jplus = j + 1;
			size_t indxEast = x_new(jplus, i, width);
			label = labelArray[indxEast];
			if (label != 1) {
				x = selectX(distancePtr, height, width, i, jplus);
				y = selectY(distancePtr, height, width, i, jplus);
				coefSpeed = potentialPtr[indxEast];
				dataType dEast = solve2dQuadratic(x, y, coefSpeed);
				if (label == 3) {
					distancePtr[indxEast] = dEast;
					labelArray[indxEast] = 2;
					pointFastMarching2D EastNeighbor = { jplus, i, dEast };
					addPointHeap2D(inProcess, EastNeighbor);
				}
				else {
					if (label == 2) {
						if (dEast < distancePtr[indxEast]) {
							distancePtr[indxEast] = dEast;
						}
					}
				}
			}
		}

		//West
		if (j > 0 && i >= 0 && i < height) {
			size_t jminus = j - 1;
			size_t indxWest = x_new(jminus, i, width);
			label = labelArray[indxWest];
			if (label != 1) {
				x = selectX(distancePtr, height, width, i, jminus);
				y = selectY(distancePtr, height, width, i, jminus);
				coefSpeed = potentialPtr[indxWest];
				dataType dWest = solve2dQuadratic(x, y, coefSpeed);
				if (label == 3) {
					distancePtr[indxWest] = dWest;
					labelArray[indxWest] = 2;
					pointFastMarching2D WestNeighbor = { jminus, i, dWest };
					addPointHeap2D(inProcess, WestNeighbor);
				}
				else {
					if (label == 2) {
						if (dWest < distancePtr[indxWest]) {
							distancePtr[indxWest] = dWest;
						}
					}
				}
			}
		}

		//North
		if (j >= 0 && j < width && i > 0) {
			size_t iminus = i - 1;
			size_t indxNorth = x_new(j, iminus, width);
			label = labelArray[indxNorth];
			if (label != 1) {
				x = selectX(distancePtr, height, width, iminus, j);
				y = selectY(distancePtr, height, width, iminus, j);
				coefSpeed = potentialPtr[indxNorth];
				dataType dNorth = solve2dQuadratic(x, y, coefSpeed);
				if (label == 3) {
					distancePtr[indxNorth] = dNorth;
					labelArray[indxNorth] = 2;
					pointFastMarching2D NorthNeighbor = { j, iminus, dNorth };
					addPointHeap2D(inProcess, NorthNeighbor);
				}
				else {
					if (label == 2) {
						if (dNorth < distancePtr[indxNorth]) {
							distancePtr[indxNorth] = dNorth;
						}
					}
				}
			}


		}

		//South
		if (j >= 0 && j < width && i < height - 1) {
			size_t iplus = i + 1;
			size_t indxSouth = x_new(j, iplus, width);
			label = labelArray[indxSouth];
			if (label != 1) {
				x = selectX(distancePtr, height, width, iplus, j);
				y = selectY(distancePtr, height, width, iplus, j);
				coefSpeed = potentialPtr[indxSouth];
				dataType dSouth = solve2dQuadratic(x, y, coefSpeed);
				if (label == 3) {
					distancePtr[indxSouth] = dSouth;
					labelArray[indxSouth] = 2;
					pointFastMarching2D SouthNeighbor = { j, iplus, dSouth };
					addPointHeap2D(inProcess, SouthNeighbor);
				}
				else {
					if (label == 2) {
						if (dSouth < distancePtr[indxSouth]) {
							distancePtr[indxSouth] = dSouth;
						}
					}
				}
			}
		}

	}
}

bool shortestPath2d(dataType* distanceFuncPtr, dataType* resultedPath, const size_t height, const size_t width, dataType h, point2D* seedPoints) {

	if (distanceFuncPtr == NULL || resultedPath == NULL || seedPoints == NULL)
		return false;

	size_t i, j, xd, dim2d = height * width;
	dataType tau = 0.8, dist_min = INFINITY, tol = 1.0;
	size_t i_init = seedPoints[0].y, j_init = seedPoints[0].x, i_end = seedPoints[1].y, j_end = seedPoints[1].x;
	size_t i_current, j_current;
	dataType iNew = 0.0;
	dataType jNew = 0.0;

	dataType * gradientVectorX = new dataType[dim2d];
	dataType * gradientVectorY = new dataType[dim2d];

	//Normalization of the gradient
	computeImageGradient(distanceFuncPtr, gradientVectorX, gradientVectorY, height, width, 1.0);
	dataType norm_of_gradient = computeGradientNorm2d(gradientVectorX, gradientVectorY, height, width);
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			xd = x_new(j, i, width);
			gradientVectorX[xd] = gradientVectorX[xd] / norm_of_gradient;
			gradientVectorY[xd] = gradientVectorY[xd] / norm_of_gradient;
		}
	}


	//Find the closest point till the last point
	size_t cpt = 1;
	i_current = i_end; j_current = j_end;
	resultedPath[x_new(j_current, i_current, width)] = 1;

	// Make the end point visible on the result
	//===============
	resultedPath[x_new(j_init, i_init, width)] = 1; 
	resultedPath[x_new(j_init, i_init - 1, width)] = 1;
	resultedPath[x_new(j_init, i_init + 1, width)] = 1;
	resultedPath[x_new(j_init - 1, i_init, width)] = 1; 
	resultedPath[x_new(j_init - 1, i_init - 1, width)] = 1;
	resultedPath[x_new(j_init - 1, i_init + 1, width)] = 1; 
	resultedPath[x_new(j_init + 1, i_init, width)] = 1; 
	resultedPath[x_new(j_init + 1, i_init - 1, width)] = 1;
	resultedPath[x_new(j_init + 1, i_init + 1, width)] = 1;
	//===============

	iNew = i_current; jNew = j_current;
	dataType currentDist = 0;

	do{

		iNew = iNew - tau * gradientVectorY[x_new(j_current, i_current, width)];
		jNew = jNew - tau * gradientVectorX[x_new(j_current, i_current, width)];

		dist_min = sqrt((iNew - i_init) * (iNew - i_init) + (jNew - j_init) * (jNew - j_init));

		i_current = round(iNew); j_current = round(jNew);
		resultedPath[x_new(j_current, i_current, width)] = 1;
		
		currentDist = distanceFuncPtr[x_new(j_current, i_current, width)];

		cpt++;
	}
	while(dist_min > tol && cpt < 10000000);

	cout << "\nDistance to the end point : " << dist_min << endl;
	cout << "\nNumber of iterations : " << cpt << endl;

	delete[] gradientVectorX;
	delete[] gradientVectorY;

	return true;
}

//===========================================================

//Functions for 3D images
// 3U^2 - 2U(X+Y+Z) + (X^2 + Y^2 + Z^2 - W) = 0 ---> aU + 2bU + c = 0
dataType solve3dQuadratic(dataType X, dataType Y, dataType Z, dataType W) {

	dataType sol = 0.0, a = 0.0, b = 0.0, c = 0.0, delta = 0.0;

	if (X == INFINITY && Y != INFINITY && Z != INFINITY) {
		a = 2;
		b = -2 * (Y + Z);
		c = pow(Y, 2) + pow(Z, 2) - W;
		delta = b * b - 4 * a * c;
		if (delta >= 0) {
			return (-b + sqrt(delta)) / (2 * a);
		}
		else {
			return min(Y, Z) + W;
		}
	}

	if (Y == INFINITY && X != INFINITY && Z != INFINITY) {
		a = 2;
		b = -2 * (X + Z);
		c = pow(X, 2) + pow(Z, 2) - W;
		delta = b * b - 4 * a * c;
		if (delta >= 0) {
			return (-b + sqrt(delta)) / (2 * a);
		}
		else {
			return min(X, Z) + W;
		}
	}

	if (Z == INFINITY && X != INFINITY && Y != INFINITY) {
		a = 2;
		b = -2 * (X + Y);
		c = pow(X, 2) + pow(Y, 2) - W;
		delta = b * b - 4 * a * c;
		if (delta >= 0) {
			return (-b + sqrt(delta)) / (2 * a);
		}
		else {
			return min(X, Y) + W;
		}
	}

	if (X == INFINITY && Y == INFINITY && Z != INFINITY) {
		a = 1;
		b = -2 * Z;
		c = pow(Z, 2) - W;
		delta = b * b - 4 * a * c;
		if (delta >= 0) {
			return (-b + sqrt(delta)) / (2 * a);
		}
		else {
			return Z + W;
		}
	}

	if (X == INFINITY && Z == INFINITY && Y != INFINITY) {
		a = 1;
		b = -2 * Y;
		c = pow(Y, 2) - W;
		delta = b * b - 4 * a * c;
		if (delta >= 0) {
			return (-b + sqrt(delta)) / (2 * a);
		}
		else {
			return Y + W;
		}
	}

	if (Y == INFINITY && Z == INFINITY && X != INFINITY) {
		a = 1;
		b = -2 * X;
		c = pow(X, 2) - W;
		delta = b * b - 4 * a * c;
		if (delta >= 0) {
			return (-b + sqrt(delta)) / (2 * a);
		}
		else {
			return X + W;
		}
	}

	if (X != INFINITY && Y != INFINITY && Z != INFINITY) {
		a = 3;
		b = -2 * (X + Y + Z);
		c = pow(X, 2) + pow(Y, 2) + pow(Z, 2) - W;
		delta = b * b - 4 * a * c;
		if (delta >= 0) {
			return (-b + sqrt(delta)) / (2 * a);
		}
		else {
			return min(X, min(Y, Z)) + W;
		}
	}

}

dataType select3dX(dataType** distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t dimK, const size_t I, const size_t J, const size_t K) {

	dataType j_minus, j_plus;

	if (J == 0) {
		j_minus = BIG_VALUE;
	}
	else {
		j_minus = distanceFuncPtr[K][x_new(J - 1, I, dimJ)];
	}

	if (J == dimJ - 1) {
		j_plus = BIG_VALUE;
	}
	else {
		j_plus = distanceFuncPtr[K][x_new(J + 1, I, dimJ)];
	}

	return min(j_minus, j_plus);
}

dataType select3dY(dataType** distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t dimK, const size_t I, const size_t J, const size_t K) {

	dataType i_minus, i_plus;

	if (I == 0) {
		i_minus = BIG_VALUE;
	}
	else {
		i_minus = distanceFuncPtr[K][x_new(J, I - 1, dimJ)];
	}

	if (I == dimI - 1) {
		i_plus = BIG_VALUE;
	}
	else {
		i_plus = distanceFuncPtr[K][x_new(J, I + 1, dimJ)];
	}

	return min(i_minus, i_plus);
}

dataType select3dZ(dataType** distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t dimK, const size_t I, const size_t J, const size_t K) {

	dataType k_minus, k_plus;

	size_t xd = x_new(J, I, dimJ);

	if (K == 0) {
		k_minus = BIG_VALUE;
	}
	else {
		k_minus = distanceFuncPtr[K - 1][xd];
	}

	if (K == dimK - 1) {
		k_plus = BIG_VALUE;
	}
	else {
		k_plus = distanceFuncPtr[K + 1][xd];
	}

	return min(k_minus, k_plus);
}

dataType computeGradientNorm3d(dataType** gradientVectorX, dataType** gradientVectorY, dataType** gradientVectorZ, const size_t length, const size_t width, const size_t height) {
	
	size_t i, j, k, xd;
	dataType norm_array = 0.0;

	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				xd = x_new(j, i, width);
				norm_array = norm_array + pow(gradientVectorX[k][xd], 2) + pow(gradientVectorY[k][xd], 2) + pow(gradientVectorZ[k][xd], 2);
			}
		}
	}
	
	return sqrt(norm_array);
}

bool compute3dImageGradient(dataType** imageDataPtr, dataType** gradientVectorX, dataType** gradientVectorY, dataType** gradientVectorZ, const size_t lenght, const size_t width, const size_t height, dataType h) {

	if (imageDataPtr == NULL || gradientVectorX == NULL || gradientVectorY == NULL || gradientVectorZ == NULL) {
		return false;
	}

	size_t i, j, k, currentInd;
	dataType ux = 0.0, uy = 0.0, uz = 0.0;

	for (k = 0; k < height; k++) {
		for (i = 0; i < lenght; i++) {
			for (j = 0; j < width; j++) {

				currentInd = x_new(j, i, width);

				if (k == 0) {
					uz = (imageDataPtr[k + 1][currentInd] - imageDataPtr[k][currentInd]) / h;
				}
				else {
					if (k == (height - 1)) {
						uz = (imageDataPtr[k][currentInd] - imageDataPtr[k - 1][currentInd]) / h;
					}
					else {
						uz = (imageDataPtr[k + 1][currentInd] - imageDataPtr[k - 1][currentInd]) / (2 * h);
					}
				}

				if (i == 0) {
					uy = (imageDataPtr[k][x_new(j, i + 1, width)] - imageDataPtr[k][x_new(j, i, width)]) / h;
				}
				else {
					if (i == (lenght - 1)) {
						uy = (imageDataPtr[k][x_new(j, i, width)] - imageDataPtr[k][x_new(j, i - 1, width)]) / h;
					}
					else {
						uy = (imageDataPtr[k][x_new(j, i + 1, width)] - imageDataPtr[k][x_new(j, i - 1, width)]) / (2 * h);
					}
				}

				if (j == 0) {
					ux = (imageDataPtr[k][x_new(j + 1, i, width)] - imageDataPtr[k][x_new(j, i, width)]) / h;
				}
				else {
					if (j == (width - 1)) {
						ux = (imageDataPtr[k][x_new(j, i, width)] - imageDataPtr[k][x_new(j - 1, i, width)]) / h;
					}
					else {
						ux = (imageDataPtr[k][x_new(j + 1, i, width)] - imageDataPtr[k][x_new(j - 1, i, width)]) / (2 * h);
					}
				}

				gradientVectorX[k][currentInd] = ux; // j
				gradientVectorY[k][currentInd] = uy; // i
				gradientVectorZ[k][currentInd] = uz; // k

			}
		}
	}
	return true;
}

bool compute3dPotential(dataType** imageDataPtr, dataType** potentialFuncPtr, const size_t length, const size_t width, const size_t height, point3d* seedPoints) {

	if (imageDataPtr == NULL || potentialFuncPtr == NULL || seedPoints == NULL)
		return false;

	size_t i, j, k;
	const size_t dim2D = length * width;
	size_t i0 = seedPoints[0].y, j0 = seedPoints[0].x, k0 = seedPoints[0].z;
	size_t i1 = seedPoints[1].y, j1 = seedPoints[1].x, k1 = seedPoints[1].z;

	dataType** gradientVectorX = new dataType* [height];
	dataType** gradientVectorY = new dataType* [height];
	dataType** gradientVectorZ = new dataType* [height];
	for (k = 0; k < height; k++) {
		gradientVectorX[k] = new dataType [dim2D];
		gradientVectorY[k] = new dataType [dim2D];
		gradientVectorZ[k] = new dataType [dim2D];
	}
	if (gradientVectorX == NULL || gradientVectorY == NULL || gradientVectorZ == NULL)
		return false;
	
	compute3dImageGradient(imageDataPtr, gradientVectorX, gradientVectorY, gradientVectorZ, length, width, height, 1.0);

	size_t seedIndice = x_new(j0, i0, width), currentIndx = 0;
	dataType seedVal = (imageDataPtr[k0][x_new(j0, i0, width)] + imageDataPtr[k1][x_new(j1, i1, width)]) / 2;
	dataType ux = 0.0, uy = 0.0, uz = 0.0;
	dataType epsilon = 0.01, K = 0.00005;

	//Computation of potential function
	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				currentIndx = x_new(j, i, width);
				potentialFuncPtr[k][currentIndx] = abs(seedVal - imageDataPtr[k][currentIndx]);
			}
		}
	}

	//Find max
	dataType max_potential = -1 * INFINITY;
	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				currentIndx = x_new(j, i, width);
				if (potentialFuncPtr[k][currentIndx] > max_potential) {
					max_potential = potentialFuncPtr[k][currentIndx];
				}
			}
		}
	}

	//Normalization
	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				currentIndx = x_new(j, i, width);
				ux = gradientVectorX[k][currentIndx];
				uy = gradientVectorY[k][currentIndx];
				uz = gradientVectorZ[k][currentIndx];
				potentialFuncPtr[k][currentIndx] = epsilon + (potentialFuncPtr[k][currentIndx] / max_potential) * (1 + K * (ux * ux + uy * uy + uz * uz));
				//potentialFuncPtr[k][currentIndx] = (dataType)(1./ (1 + sqrt(ux * ux + uy * uy + uz * uz)));
				//potentialFuncPtr[k][currentIndx] = (dataType)(1.0 / (10 + abs(seedVal - imageDataPtr[k][currentIndx])));
			}
		}
	}

	for (k = 0; k < height; k++) {
		delete[] gradientVectorX[k];
		delete[] gradientVectorY[k];
		delete[] gradientVectorZ[k];
	}
	delete[] gradientVectorX; 
	delete[] gradientVectorY; 
	delete[] gradientVectorZ;

	return true;
}

//function to swap 3D points in the fast marching contest
void swap3dPoints(pointFastMarching3D* a, pointFastMarching3D* b) {
	pointFastMarching3D temp = *a;
	*a = *b;
	*b = temp;
}

void heapifyDown3D(vector<pointFastMarching3D>& in_Process, int i) {

	int length_array = in_Process.size();
	int current = i;
	int left_child = 2 * i + 1;
	int right_child = 2 * i + 2;

	dataType val_current = 0.0, val_left = 0.0, val_right = 0.0;

	if (current >= 0 && current < length_array) {
		val_current = in_Process[current].arrival;
	}

	if (left_child < length_array) {
		val_left = in_Process[left_child].arrival;
		if (val_left < val_current) {
			current = left_child;
			val_current = in_Process[current].arrival;
		}
	}

	if (right_child < length_array) {
		val_right = in_Process[right_child].arrival;
		if (val_right < val_current) {
			current = right_child;
		}
	}

	if (current != i) {
		swap3dPoints(&in_Process[i], &in_Process[current]);
		heapifyDown3D(in_Process, current);
	}

}

void heapifyVector3D(vector<pointFastMarching3D>& in_Process) {
	int length_array = in_Process.size();
	int ind, start = length_array / 2 - 1;
	for (ind = start; ind >= 0; ind--) {
		heapifyDown3D(in_Process, ind);
	}
}

void deleteRootHeap3D(vector<pointFastMarching3D>& in_Process) {
	//we use type int for indexes because we do operations like pos--
	int l = in_Process.size();
	swap3dPoints(&in_Process[0], &in_Process[l - 1]);
	in_Process.pop_back();
	heapifyDown3D(in_Process, 0);
}

void addPointHeap3D(vector<pointFastMarching3D>& in_Process, pointFastMarching3D point) {
	//we use type int for indexes because we do operations like pos--
	in_Process.push_back(point);
	int l = in_Process.size();
	heapifyUp3D(in_Process, l - 1);
}

void heapifyUp3D(vector<pointFastMarching3D>& in_Process, int i) {

	int current = i;

	if (i > 0) {
		int parent = (i - 1) / 2;
		dataType val_current = in_Process[current].arrival;
		dataType val_parent = in_Process[parent].arrival;
		if (val_current < val_parent) {
			current = parent;
		}
	}

	if (current != i) {
		swap3dPoints(&in_Process[current], &in_Process[i]);
		heapifyUp3D(in_Process, current);
	}

}

bool fastMarching3D_N(dataType** imageDataPtr, dataType** distanceFuncPtr, dataType** potentialFuncPtr, const size_t length, const size_t width, const size_t height, point3d* seedPoints) {

	if (imageDataPtr == NULL || distanceFuncPtr == NULL || potentialFuncPtr == NULL || seedPoints == NULL) {
		return false;
	}

	vector <pointFastMarching3D> inProcess;
	size_t i = 0, j = 0, k = 0;

	size_t** labelArray = new size_t * [height];
	for (k = 0; k < height; k++) {
		labelArray[k] = new size_t[width * length];
	}
	if (labelArray == NULL)
		return false;

	dataType x = 0.0, y = 0.0, z = 0.0;
	size_t currentIndx = 0;


	//Initialization
	//All the points are notProcessed ---> label = 3
	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				currentIndx = x_new(j, i, width);
				distanceFuncPtr[k][currentIndx] = INFINITY;
				labelArray[k][currentIndx] = 3;
			}
		}
	}

	compute3dPotential(imageDataPtr, potentialFuncPtr, length, width, height, seedPoints);

	//Processed the initial point
	pointFastMarching3D current;
	j = seedPoints[0].x;
	i = seedPoints[0].y;
	k = seedPoints[0].z;
	currentIndx = x_new(j, i, width);
	distanceFuncPtr[k][currentIndx] = 0.0;
	labelArray[k][currentIndx] = 1;

	//find the neighbours of the initial point add add them to inProcess
	size_t height_minus = height - 1, length_minus = length - 1, width_minus = width - 1;
	size_t kminus = 0, kplus = 0, iminus = 0, iplus = 0, jminus = 0, jplus = 0;
	size_t indxNorth = 0, indxSouth = 0, indxWest = 0, indxEast = 0;
	dataType dTop = 0.0, dBottom = 0.0, dNorth = 0.0, dSouth = 0.0, dWest = 0.0, dEast = 0.0, coefSpeed = 0.0;

	//Top
	if (k > 0 && j >= 0 && j < width && i >= 0 && i < length) {
		kminus = k - 1;
		x = select3dX(distanceFuncPtr, length, width, height, i, j, kminus);
		y = select3dY(distanceFuncPtr, length, width, height, i, j, kminus);
		z = select3dZ(distanceFuncPtr, length, width, height, i, j, kminus);
		coefSpeed = potentialFuncPtr[kminus][currentIndx];
		dTop = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D TopNeighbor = { j, i, kminus, dTop };
		distanceFuncPtr[kminus][currentIndx] = dTop;
		inProcess.push_back(TopNeighbor);
		labelArray[kminus][currentIndx] = 2;
	}

	//Bottom
	if (k < height_minus && j >= 0 && j < width && i >= 0 && i < length) {
		kplus = k + 1;
		x = select3dX(distanceFuncPtr, length, width, height, i, j, kplus);
		y = select3dY(distanceFuncPtr, length, width, height, i, j, kplus);
		z = select3dZ(distanceFuncPtr, length, width, height, i, j, kplus);
		coefSpeed = potentialFuncPtr[kplus][currentIndx];
		dBottom = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D BottomNeighbor = { j, i, kplus, dBottom };
		distanceFuncPtr[kplus][currentIndx] = dBottom;
		inProcess.push_back(BottomNeighbor);
		labelArray[kplus][currentIndx] = 2;
	}

	//North
	if (k >= 0 && k < height && j > 0 && i >= 0 && i < length) {
		jminus = j - 1;
		indxNorth = x_new(jminus, i, width);
		x = select3dX(distanceFuncPtr, length, width, height, i, jminus, k);
		y = select3dY(distanceFuncPtr, length, width, height, i, jminus, k);
		z = select3dZ(distanceFuncPtr, length, width, height, i, jminus, k);
		coefSpeed = potentialFuncPtr[k][indxNorth];
		dNorth = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D NorthNeighbor = { jminus, i, k, dNorth };
		distanceFuncPtr[k][indxNorth] = dNorth;
		inProcess.push_back(NorthNeighbor);
		labelArray[k][indxNorth] = 2;
	}

	//South
	if (k >= 0 && k < height && j < width_minus && i >= 0 && i < length) {
		jplus = j + 1;
		indxSouth = x_new(jplus, i, width);
		x = select3dX(distanceFuncPtr, length, width, height, i, jplus, k);
		y = select3dY(distanceFuncPtr, length, width, height, i, jplus, k);
		z = select3dZ(distanceFuncPtr, length, width, height, i, jplus, k);
		coefSpeed = potentialFuncPtr[k][indxSouth];
		dSouth = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D SouthNeighbor = { jplus, i, k, dSouth };
		distanceFuncPtr[k][indxSouth] = dSouth;
		inProcess.push_back(SouthNeighbor);
		labelArray[k][indxSouth] = 2;
	}

	//East
	if (k >= 0 && k < height && j >= 0 && j < width && i < length_minus) {
		iplus = i + 1;
		indxEast = x_new(j, iplus, width);
		x = select3dX(distanceFuncPtr, length, width, height, iplus, j, k);
		y = select3dY(distanceFuncPtr, length, width, height, iplus, j, k);
		z = select3dZ(distanceFuncPtr, length, width, height, iplus, j, k);
		coefSpeed = potentialFuncPtr[k][indxEast];
		dEast = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D EastNeighbor = { j, iplus, k, dEast };
		distanceFuncPtr[k][indxEast] = dEast;
		inProcess.push_back(EastNeighbor);
		labelArray[k][indxEast] = 2;
	}

	//West
	if (k >= 0 && k < height && j >= 0 && j < width && i > 0) {
		iminus = i - 1;
		indxWest = x_new(j, iminus, width);
		x = select3dX(distanceFuncPtr, length, width, height, iminus, j, k);
		y = select3dY(distanceFuncPtr, length, width, height, iminus, j, k);
		z = select3dZ(distanceFuncPtr, length, width, height, iminus, j, k);
		coefSpeed = potentialFuncPtr[k][indxWest];
		dWest = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D WestNeighbor = { j, iminus, k, dWest };
		distanceFuncPtr[k][indxWest] = dWest;
		inProcess.push_back(WestNeighbor);
		labelArray[k][indxWest] = 2;
	}

	//int n = inProcess.size();

	heapifyVector3D(inProcess);
	size_t label = 0;

	int l = 0, m = 0;

	while (inProcess.size() != 0) {

		//processed the minimum
		current = inProcess[0]; // index 0 exist because, if not we will be out of the while loop
		j = current.x;
		i = current.y;
		k = current.z;
		currentIndx = x_new(j, i, width);
		labelArray[k][currentIndx] = 1;
		distanceFuncPtr[k][currentIndx] = current.arrival;
		//l = inProcess.size();
		//swap3dPoints(&inProcess[0], &inProcess[l - 1]);
		//inProcess.pop_back();
		//heapifyDown3D(inProcess, 0);
		deleteRootHeap3D(inProcess);

		//Treat neighbors of the minimum in the narrow band
		//Bottom
		if (k < height_minus && i >= 0 && i < length && j >= 0 && j < width) {
			kplus = k + 1;
			label = labelArray[kplus][currentIndx];
			if (label != 1) {
				x = select3dX(distanceFuncPtr, length, width, height, i, j, kplus);
				y = select3dY(distanceFuncPtr, length, width, height, i, j, kplus);
				z = select3dZ(distanceFuncPtr, length, width, height, i, j, kplus);
				coefSpeed = potentialFuncPtr[kplus][currentIndx];
				dBottom = solve3dQuadratic(x, y, z, coefSpeed);
				if (label == 3) {
					distanceFuncPtr[kplus][currentIndx] = dBottom;
					labelArray[kplus][currentIndx] = 2;
					pointFastMarching3D BottomNeighbor = { j, i, kplus, dBottom };
					//inProcess.push_back(BottomNeighbor);
					//l = inProcess.size();
					//heapifyUp3D(inProcess, l - 1);
					addPointHeap3D(inProcess, BottomNeighbor);
				}
				else {
					if (label == 2) {
						if (dBottom < distanceFuncPtr[kplus][currentIndx]) {
							distanceFuncPtr[kplus][currentIndx] = dBottom;
						}
					}
				}
			}
		}

		//Top
		if (k > 0 && j >= 0 && j < width && i >= 0 && i < length) {
			kminus = k - 1;
			label = labelArray[kminus][currentIndx];
			if (label != 1) {
				x = select3dX(distanceFuncPtr, length, width, height, i, j, kminus);
				y = select3dY(distanceFuncPtr, length, width, height, i, j, kminus);
				z = select3dZ(distanceFuncPtr, length, width, height, i, j, kminus);
				coefSpeed = potentialFuncPtr[kminus][currentIndx];
				dTop = solve3dQuadratic(x, y, z, coefSpeed);
				if (label == 3) {
					distanceFuncPtr[kminus][currentIndx] = dTop;
					labelArray[kminus][currentIndx] = 2;
					pointFastMarching3D TopNeighbor = { j, i, kminus, dTop };
					//inProcess.push_back(TopNeighbor);
					//l = inProcess.size();
					//heapifyUp3D(inProcess, l - 1);
					addPointHeap3D(inProcess, TopNeighbor);
				}
				else {
					if (label == 2) {
						if (dTop < distanceFuncPtr[kminus][currentIndx]) {
							distanceFuncPtr[kminus][currentIndx] = dTop;
						}
					}
				}
			}
		}

		//East
		if (i < length_minus && j >= 0 && j < width && k >= 0 && k < height) {
			iplus = i + 1;
			indxEast = x_new(j, iplus, width);
			label = labelArray[k][indxEast];
			if (label != 1) {
				x = select3dX(distanceFuncPtr, length, width, height, iplus, j, k);
				y = select3dY(distanceFuncPtr, length, width, height, iplus, j, k);
				z = select3dZ(distanceFuncPtr, length, width, height, iplus, j, k);
				coefSpeed = potentialFuncPtr[k][indxEast];
				dEast = solve3dQuadratic(x, y, z, coefSpeed);
				if (label == 3) {
					distanceFuncPtr[k][indxEast] = dEast;
					labelArray[k][indxEast] = 2;
					pointFastMarching3D EastNeighbor = { j, iplus, k, dEast };
					//inProcess.push_back(EastNeighbor);
					//l = inProcess.size();
					//heapifyUp3D(inProcess, l - 1);
					addPointHeap3D(inProcess, EastNeighbor);
				}
				else {
					if (label == 2) {
						if (dEast < distanceFuncPtr[k][indxEast]) {
							distanceFuncPtr[k][indxEast] = dEast;
						}
					}
				}
			}
		}

		//West
		if (i > 0 && j >= 0 && j < width && k >= 0 && k < height) {
			iminus = i - 1;
			indxWest = x_new(j, iminus, width);
			label = labelArray[k][indxWest];
			if (label != 1) {
				x = select3dX(distanceFuncPtr, length, width, height, iminus, j, k);
				y = select3dY(distanceFuncPtr, length, width, height, iminus, j, k);
				z = select3dZ(distanceFuncPtr, length, width, height, iminus, j, k);
				coefSpeed = potentialFuncPtr[k][indxWest];
				dWest = solve3dQuadratic(x, y, z, coefSpeed);
				if (label == 3) {
					distanceFuncPtr[k][indxWest] = dWest;
					labelArray[k][indxWest] = 2;
					pointFastMarching3D WestNeighbor = { j, iminus, k, dWest };
					//inProcess.push_back(WestNeighbor);
					//l = inProcess.size();
					//heapifyUp3D(inProcess, l - 1);
					addPointHeap3D(inProcess, WestNeighbor);
				}
				else {
					if (label == 2) {
						if (dWest < distanceFuncPtr[k][indxWest]) {
							distanceFuncPtr[k][indxWest] = dWest;
						}
					}
				}
			}
		}

		//North
		if (j > 0 && i >= 0 && i < length && k >= 0 && k < height) {
			jminus = j - 1;
			indxNorth = x_new(jminus, i, width);
			label = labelArray[k][indxNorth];
			if (label != 1) {
				x = select3dX(distanceFuncPtr, length, width, height, i, jminus, k);
				y = select3dY(distanceFuncPtr, length, width, height, i, jminus, k);
				z = select3dZ(distanceFuncPtr, length, width, height, i, jminus, k);
				coefSpeed = potentialFuncPtr[k][indxNorth];
				dNorth = solve3dQuadratic(x, y, z, coefSpeed);
				if (label == 3) {
					distanceFuncPtr[k][indxNorth] = dNorth;
					labelArray[k][indxNorth] = 2;
					pointFastMarching3D NorthNeighbor = { jminus, i, k, dNorth };
					//inProcess.push_back(NorthNeighbor);
					//l = inProcess.size();
					//heapifyUp3D(inProcess, l - 1);
					addPointHeap3D(inProcess, NorthNeighbor);
				}
				else {
					if (label == 2) {
						if (dNorth < distanceFuncPtr[k][indxNorth]) {
							distanceFuncPtr[k][indxNorth] = dNorth;
						}
					}
				}
			}
		}

		//South
		if (j < width_minus && i >= 0 && i < length && k >= 0 && k < height) {
			jplus = j + 1;
			indxSouth = x_new(jplus, i, width);
			label = labelArray[k][indxSouth];
			if (label != 1) {
				x = select3dX(distanceFuncPtr, length, width, height, i, jplus, k);
				y = select3dY(distanceFuncPtr, length, width, height, i, jplus, k);
				z = select3dZ(distanceFuncPtr, length, width, height, i, jplus, k);
				coefSpeed = potentialFuncPtr[k][indxSouth];
				dSouth = solve3dQuadratic(x, y, z, coefSpeed);
				if (label == 3) {
					distanceFuncPtr[k][indxSouth] = dSouth;
					labelArray[k][indxSouth] = 2;
					pointFastMarching3D SouthNeighbor = { jplus, i, k, dSouth };
					//inProcess.push_back(SouthNeighbor);
					//l = inProcess.size();
					//heapifyUp3D(inProcess, l - 1);
					addPointHeap3D(inProcess, SouthNeighbor);
				}
				else {
					if (label == 2) {
						if (dSouth < distanceFuncPtr[k][indxSouth]) {
							distanceFuncPtr[k][indxSouth] = dSouth;
						}
					}
				}
			}
		}
	}

	//size_t cpt1 = 0, cpt2 = 0,cpt3 = 0 ;
	//for (k = 0; k < height; k++) {
	//	for (i = 0; i < length; i++) {
	//		for (j = 0; j < width; j++) {
	//			currentIndx = x_new(j, i, width);
	//			if (distanceFuncPtr[k][currentIndx] == INFINITY) {
	//				cpt1++;
	//			}
	//			if (labelArray[k][currentIndx] == 2) {
	//				cpt2++;
	//			}
	//			if (labelArray[k][currentIndx] == 3) {
	//				cpt3++;
	//			}
	//		}
	//	}
	//}
	//cout << cpt1 << " point(s) with distance = INFINITY" << endl;;
	//cout << cpt2 << " point(s) with label = 2" << endl;
	//cout << cpt3 << " point(s) with label = 3" << endl;

	for (k = 0; k < height; k++) {
		delete[] labelArray[k];
	}
	delete[] labelArray;

	return true;
}

bool shortestPath3d(dataType** distanceFuncPtr, dataType** resultedPath, const size_t length, const size_t width, const size_t height, dataType h, point3d* seedPoints) {

	if (distanceFuncPtr == NULL || resultedPath == NULL || seedPoints == NULL)
		return false;

	size_t i, j, k, xd, dim2d = length * width, max_iter = 1000000000; //width* length* height;
	dataType tau = 0.8, tol = 1.0;
	size_t i_init = seedPoints[0].y, j_init = seedPoints[0].x, k_init = seedPoints[0].z;
	size_t i_end = seedPoints[1].y, j_end = seedPoints[1].x, k_end = seedPoints[1].z;

	dataType** gradientVectorX = new dataType * [height];
	dataType** gradientVectorY = new dataType * [height];
	dataType** gradientVectorZ = new dataType * [height];
	for (k = 0; k < height; k++) {
		gradientVectorX[k] = new dataType [dim2d];
		gradientVectorY[k] = new dataType [dim2d];
		gradientVectorZ[k] = new dataType [dim2d];
	}
	if (gradientVectorX == NULL || gradientVectorY == NULL || gradientVectorZ == NULL)
		return false;

	//Normalization of the gradient
	compute3dImageGradient(distanceFuncPtr, gradientVectorX, gradientVectorY, gradientVectorZ, length, width, height, 1.0);
	dataType norm_of_gradient = computeGradientNorm3d(gradientVectorX, gradientVectorY, gradientVectorZ, length, width, height);

	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				xd = x_new(j, i, width);
				gradientVectorX[k][xd] = gradientVectorX[k][xd] / norm_of_gradient;
				gradientVectorY[k][xd] = gradientVectorY[k][xd] / norm_of_gradient;
				gradientVectorZ[k][xd] = gradientVectorZ[k][xd] / norm_of_gradient;
			}
		}
	}
	
	//Find the closest point till the last point
	size_t cpt = 0;
	size_t i_current = i_end;
	size_t j_current = j_end;
	size_t k_current = k_end;
	size_t currentIndx = x_new(j_current, i_current, width);
	resultedPath[k_current][currentIndx] = 1;

	dataType iNew = i_current;
	dataType jNew = j_current;
	dataType kNew = k_current;
	dataType currentDist = 0.0;
	dataType dist_min = 0.0;

	do {

		currentIndx = x_new(j_current, i_current, width);
		iNew = iNew - tau * gradientVectorY[k_current][currentIndx];
		jNew = jNew - tau * gradientVectorX[k_current][currentIndx];
		kNew = kNew - tau * gradientVectorZ[k_current][currentIndx];

		dist_min = sqrt((iNew - i_init) * (iNew - i_init) + (jNew - j_init) * (jNew - j_init) + (kNew - k_init) * (kNew - k_init));

		i_current = round(iNew); 
		j_current = round(jNew);
		k_current = round(kNew);
		resultedPath[k_current][x_new(j_current, i_current, width)] = 1;

		//currentDist = distanceFuncPtr[k_current][x_new(j_current, i_current, width)];
		cpt++;

	} while (dist_min > tol && cpt < max_iter);

	cout << "\nDistance to the end point : " << dist_min << endl;
	cout << "\nNumber of iterations : " << cpt << endl;

	for (k = 0; k < height; k++) {
		delete[] gradientVectorX[k];
		delete[] gradientVectorY[k];
		delete[] gradientVectorZ[k];
	}
	delete[] gradientVectorX;
	delete[] gradientVectorY;
	delete[] gradientVectorZ;

	return true;
}

