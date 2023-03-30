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
// aU^2 -2U(X+Y) + (X^2 + Y^2 - W) = 0
dataType solve2dQuadratic(dataType X, dataType Y, dataType W) {

	dataType sol = 0.0, a, b, c, delta;

	a = 2.0; 
	if (X == BIG_VALUE) {
		X = 0; a--;
	}
	if (Y == BIG_VALUE) {
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

	dataType j_minus, j_plus;

	if (J == 0) {
		j_minus = BIG_VALUE;
	}
	else {
		j_minus = distanceFuncPtr[x_new(J - 1, I, dimJ)];
	}

	if (J == dimJ - 1) {
		j_plus = BIG_VALUE;
	}
	else {
		j_plus = distanceFuncPtr[x_new(J + 1, I, dimJ)];
	}

	return min(j_minus, j_plus);
}

dataType selectY(dataType* distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t I, const size_t J) {

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

bool computePotential(dataType* imageDataPtr, dataType* potentialFuncPtr, const size_t height, const size_t width, Point2D * seedPoints)
{

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
			edgeValue = 1;//+K * (ux * ux + uy * uy);
			weight = potentialFuncPtr[currentIndx] / max_potential;
			potentialFuncPtr[currentIndx] = (dataType)(epsylon + weight * edgeValue);
		}
	}

	delete[] gradientVectorX;
	delete[] gradientVectorY;

	return true;
}

bool fastMarching2d(dataType* imageDataPtr, dataType* distanceFuncPtr, dataType* potentialFuncPtr, const size_t height, const size_t width, Point2D * seedPoints)
{
	short * labelArray = new short[height * width];

	if (imageDataPtr == NULL || distanceFuncPtr == NULL || potentialFuncPtr == NULL || seedPoints == NULL || labelArray == NULL) {
		return false;
	}

	size_t i = 0, j = 0, k = 0, dim2D = height * width, cpt = 0;
	vector<size_t> i_Processed, i_inProcess;
	vector<size_t> j_Processed, j_inProcess;
	dataType x = 0.0, y = 0.0, minSolution = 0.0, coef = 0.0, dist = 0.0;
	size_t nbNeighborsFound = 0;
	size_t h_n = 0, iSol = 0, jSol = 0;
	vector<dataType> tempTimeFunc;
	dataType dNorth = 0.0, dSouth = 0.0, dEast = 0.0, dWest = 0.0;

	//Compute the potential function
	computePotential(imageDataPtr, potentialFuncPtr, height, width, seedPoints);

	//STEP 1
	//In labelAray we have : 1 ---> already processed, 2 ---> in process and 3 ---> not processed
	for (k = 0; k < dim2D; k++) {
		distanceFuncPtr[k] = BIG_VALUE;
		labelArray[k] = 3;
	}
	//--------------------End of STEP 1 -----------------------------------

	//STEP 2
	//Add the neighbors of the starting point in the vector of pixels to be processed
	i = seedPoints[0].y; 
	j = seedPoints[0].x;
	distanceFuncPtr[x_new(j, i, width)] = 0;
	i_Processed.push_back(i); 
	j_Processed.push_back(j);

	size_t iminus = 0, iplus = 0, jminus = 0, jplus = 0;

	//East
	if (j < width - 1 && i >= 0 && i < height) {

		jplus = j + 1;

		if (labelArray[x_new(jplus, i, width)] == 3) {
			i_inProcess.push_back(i);
			j_inProcess.push_back(jplus);
			labelArray[x_new(jplus, i, width)] = 2;
			nbNeighborsFound++;
		}
	}

	//West
	if (j > 0 && i >= 0 && i < height) {

		jminus = j - 1;

		if (labelArray[x_new(jminus, i, width)] == 3) {
			i_inProcess.push_back(i); 
			j_inProcess.push_back(jminus);
			labelArray[x_new(jminus, i, width)] = 2;
			nbNeighborsFound++;
		}
	}

	//North
	if (j >= 0 && j < width && i > 0) {

		iminus = i - 1;

		if (labelArray[x_new(j, iminus, width)] == 3) {
			i_inProcess.push_back(iminus); 
			j_inProcess.push_back(j);
			labelArray[x_new(j, iminus, width)] = 2;
			nbNeighborsFound++;
		}
	}

	//South
	if (j >= 0 && j < width && i < height - 1) {

		iplus = i + 1;

		if (labelArray[x_new(j, iplus, width)] == 3) {
			i_inProcess.push_back(iplus); 
			j_inProcess.push_back(j);
			labelArray[x_new(j, iplus, width)] = 2;
			nbNeighborsFound++;
		}
	}
	

	//Compute the solution for neighbors in the stack
	h_n = i_inProcess.size();

	for (k = 0; k < h_n; k++) {
		i = i_inProcess[k];
		j = j_inProcess[k];
		x = selectX(distanceFuncPtr, height, width, i, j);
		y = selectY(distanceFuncPtr, height, width, i, j);
		coef = potentialFuncPtr[x_new(j, i, width)];
		dist = solve2dQuadratic(x, y, coef);
		tempTimeFunc.push_back(dist);
	}

	//Update the label of seed point as processed
	labelArray[x_new(j, i, width)] = 1;

	//int L = tempTimeFunc.size();
	//heapSortMaxMin(tempTimeFunc, L);

	//---------------------End of STEP 2 -------------------------------------

	//STEP 3
	while (i_inProcess.size() != 0) {

		h_n = i_inProcess.size();
		// i_inProcess and j_inProcess have the same size

		//Find the minimal solution
		minSolution = INFINITY;
		for (k = 0; k < h_n; k++) {
			if (minSolution > tempTimeFunc[k]) {
				minSolution = tempTimeFunc[k];
				iSol = k; 
				i = i_inProcess[k];
				jSol = k; 
				j = j_inProcess[k];
			}
		}

		//cout << "\ndistance  : " << minSolution << endl;
		
		//Set the distance to the processed pixel
		distanceFuncPtr[x_new(j, i, width)] = minSolution;
		labelArray[x_new(j, i, width)] = 1;
		i_Processed.push_back(i); 
		j_Processed.push_back(j);

		//Remove the processed pixel to the stack
		if (iSol == 0 || jSol == 0) {
			i_inProcess.erase(i_inProcess.begin());
			j_inProcess.erase(j_inProcess.begin());
			tempTimeFunc.erase(tempTimeFunc.begin());
		}
		else {
			i_inProcess.erase(i_inProcess.begin() + iSol);
			j_inProcess.erase(j_inProcess.begin() + jSol);
			tempTimeFunc.erase(tempTimeFunc.begin() + jSol);
		}

		//Compute solution for the neigbors of the selected point

		//STEP 4
		//Find the neighbors of the processed pixel and compute they time function

		//East
		if (j < width - 1 && i >= 0 && i < height) {

			jplus = j + 1;

			x = selectX(distanceFuncPtr, height, width, i, jplus);
			y = selectY(distanceFuncPtr, height, width, i, jplus);
			coef = potentialFuncPtr[x_new(jplus, i, width)];
			dEast = solve2dQuadratic(x, y, coef);
			if (labelArray[x_new(jplus, i, width)] == 3) {
				distanceFuncPtr[x_new(jplus, i, width)] = dEast;
				i_inProcess.push_back(i); 
				j_inProcess.push_back(jplus);
				tempTimeFunc.push_back(dEast);
				labelArray[x_new(jplus, i, width)] = 2;
			}
			else {
				if (labelArray[x_new(jplus, i, width)] == 2) {
					if (dEast < distanceFuncPtr[x_new(jplus, i, width)]) {
						distanceFuncPtr[x_new(jplus, i, width)] = dEast;
					}
				}
			}
		}

		//West
		if (j > 0 && i >= 0 && i < height) {

			jminus = j - 1;

			x = selectX(distanceFuncPtr, height, width, i, jminus);
			y = selectY(distanceFuncPtr, height, width, i, jminus);
			coef = potentialFuncPtr[x_new(jminus, i, width)];
			dWest = solve2dQuadratic(x, y, coef);
			if (labelArray[x_new(jminus, i, width)] == 3) {
				distanceFuncPtr[x_new(jminus, i, width)] = dWest;
				i_inProcess.push_back(i); 
				j_inProcess.push_back(jminus);
				tempTimeFunc.push_back(dWest);
				labelArray[x_new(jminus, i, width)] = 2;
			}
			else {
				if (labelArray[x_new(jminus, i, width)] == 2) {
					if (dWest < distanceFuncPtr[x_new(jminus, i, width)]) {
						distanceFuncPtr[x_new(jminus, i, width)] = dWest;
					}
				}
			}
		}

		//North
		if (j >= 0 && j < width && i > 0) {

			iminus = i - 1;

			x = selectX(distanceFuncPtr, height, width, iminus, j);
			y = selectY(distanceFuncPtr, height, width, iminus, j);
			coef = potentialFuncPtr[x_new(j, iminus, width)];
			dNorth = solve2dQuadratic(x, y, coef);
			if (labelArray[x_new(j, iminus, width)] == 3) {
				distanceFuncPtr[x_new(j, iminus, width)] = dNorth;
				i_inProcess.push_back(iminus); 
				j_inProcess.push_back(j);
				tempTimeFunc.push_back(dNorth);
				labelArray[x_new(j, iminus, width)] = 2;
			}
			else {
				if (labelArray[x_new(j, iminus, width)] == 2) {
					if (dNorth < distanceFuncPtr[x_new(j, iminus, width)]) {
						distanceFuncPtr[x_new(j, iminus, width)] = dNorth;
					}
				}
			}
		}

		//South
		if (j >= 0 && j < width && i < height - 1) {

			iplus = i + 1;

			x = selectX(distanceFuncPtr, height, width, iplus, j);
			y = selectY(distanceFuncPtr, height, width, iplus, j);
			coef = potentialFuncPtr[x_new(j, iplus, width)];
			dSouth = solve2dQuadratic(x, y, coef);
			if (labelArray[x_new(j, iplus, width)] == 3) {
				distanceFuncPtr[x_new(j, iplus, width)] = dSouth;
				i_inProcess.push_back(iplus); 
				j_inProcess.push_back(j);
				tempTimeFunc.push_back(dSouth);
				labelArray[x_new(j, iplus, width)] = 2;
			}
			else {
				if (labelArray[x_new(j, iplus, width)] == 2) {
					if (dSouth < distanceFuncPtr[x_new(j, iplus, width)]) {
						distanceFuncPtr[x_new(j, iplus, width)] = dSouth;
					}
				}
			}
		}

	}

	for (k = 0; k < dim2D; k++) {
		if (labelArray[k] == 3) {
			cpt++;
		}
	}
	if (cpt == 0) {
		cout << "All the pixels have been visited" << endl;
	}
	else {
		cout << "\n" << cpt << " pixels haven't been visited" << endl;
	}

	//cout << "\n" << cpt << " pixels with distance < 0 : " << endl;
	//cout << "\nNumber of processed Point :" << i_Processed.size() << endl;

	//free(labelArray);
	delete[] labelArray;

	return true;
}

bool shortestPath2d(dataType* distanceFuncPtr, dataType* resultedPath, const size_t height, const size_t width, dataType h, Point2D* seedPoints) {

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

//Functions for 3D images
// 3U^2 - 2U(X+Y+Z) + (X^2 + Y^2 + Z^2 - W) = 0 ---> aU + 2bU + c = 0
dataType solve3dQuadratic(dataType X, dataType Y, dataType Z, dataType W) {

	dataType sol = 0.0, a = 0.0, b = 0.0, c = 0.0, delta = 0.0;

	//a = 3;
	//if (X == BIG_VALUE) {
	//	a--; X = 0;
	//}
	//if (Y == BIG_VALUE) {
	//	a--; Y = 0;
	//}
	//if (Z == BIG_VALUE) {
	//	a--; Z = 0;
	//}
	//b = -2 * (X + Y + Z); 
	//c = pow(X, 2) + pow(Y, 2) + pow(Z, 2) - W;
	//delta = b * b  - 4 * a * c;
	//if (a == 0) {
	//	sol = 0;
	//}
	//else {
	//	if (delta >= 0) {
	//		sol = (-b + sqrt(delta)) / (2 * a);
	//	}
	//	else {
	//		sol = min(X, min(Y, Z)) + W;
	//	}
	//}

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

	//else {
	//	a = 3;
	//	b = -2 * (X + Y + Z);
	//	c = pow(X, 2) + pow(Y, 2) + pow(Z, 2) - W;
	//	delta = b * b - 4 * a * c;
	//	if (delta >= 0) {
	//		sol = (-b + sqrt(delta)) / (2 * a);
	//	}
	//	else {
	//		sol = min(X, min(Y, Z)) + W;
	//	}
	//}

	//if (Y == BIG_VALUE) {
	//	a = 2;
	//	b = -2 * (X + Z);
	//	c = pow(X, 2) + pow(Z, 2) - W;
	//	delta = b * b - 4 * a * c;
	//	if (delta >= 0) {
	//		sol = (-b + sqrt(delta)) / (2 * a);
	//	}
	//	else {
	//		sol = min(X, Z) + W;
	//	}
	//}
	//if (Z == BIG_VALUE) {
	//	a = 2;
	//	b = -2 * (X + Y);
	//	c = pow(X, 2) + pow(Y, 2) - W;
	//	delta = b * b - 4 * a * c;
	//	if (delta >= 0) {
	//		sol = (-b + sqrt(delta)) / (2 * a);
	//	}
	//	else {
	//		sol = min(X, Y) + W;
	//	}
	//}
	//if (X == BIG_VALUE && Y == BIG_VALUE) {
	//	a = 1;
	//	b = -2 * Z;
	//	c = pow(Z, 2) - W;
	//	delta = b * b - 4 * a * c;
	//	if (delta >= 0) {
	//		sol = (-b + sqrt(delta)) / (2 * a);
	//	}
	//	else {
	//		sol = Z + W;
	//	}
	//}
	//if (X == BIG_VALUE && Z == BIG_VALUE) {
	//	a = 1;
	//	b = -2 * Y;
	//	c = pow(Y, 2) - W;
	//	delta = b * b - 4 * a * c;
	//	if (delta >= 0) {
	//		sol = (-b + sqrt(delta)) / (2 * a);
	//	}
	//	else {
	//		sol = Y + W;
	//	}
	//}
	//if (Y == BIG_VALUE && Z == BIG_VALUE) {
	//	a = 1;
	//	b = -2 * X;
	//	c = pow(X, 2) - W;
	//	delta = b * b - 4 * a * c;
	//	if (delta >= 0) {
	//		sol = (-b + sqrt(delta)) / (2 * a);
	//	}
	//	else {
	//		sol = X + W;
	//	}
	//}
	
	//if (sol < 0) {
	//	cout << "The solution is negative " << endl; //If everything is OK, it never happen
	//	return 0;
	//}
	//else {
	//	return sol;
	//}

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

bool compute3dPotential(dataType** imageDataPtr, dataType** potentialFuncPtr, const size_t length, const size_t width, const size_t height, Point3d* seedPoints) {

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
				//potentialFuncPtr[k][currentIndx] = 1 / (1 + K * (ux * ux + uy * uy + uz * uz));
				//if (k >= h && i >= l && j >= w) {
				//	potentialFuncPtr[k][currentIndx] = 1;
				//}
				//else {
				//	potentialFuncPtr[k][currentIndx] = 2;
				//}
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

bool fastMarching3d_N(dataType** imageDataPtr, dataType** distanceFuncPtr, dataType** potentialFuncPtr, const size_t length, const size_t width, const size_t height, Point3d* seedPoints)
{
	
	if (imageDataPtr == NULL || distanceFuncPtr == NULL || potentialFuncPtr == NULL || seedPoints == NULL) {
		return false;
	}

	size_t i = 0, j = 0, k = 0, n = 0, xd = 0, dim2D = length * width, cpt = 0;
	vector<size_t> i_inProcess, j_inProcess, k_inProcess;
	dataType x = 0.0, y = 0.0, z = 0.0, minSolution = 0.0, coef = 0.0, dist = 0.0;
	size_t lenghtSol = 0, indxSol = 0;
	vector<dataType> tempTimeFunc;
	dataType dNorth = 0.0, dSouth = 0.0, dEast = 0.0, dWest = 0.0, dTop = 0.0, dBottom = 0.0;

	short** labelArray = new short* [height];
	for (k = 0; k < height; k++) {
		labelArray[k] = new short[dim2D];
	}
	if (labelArray == NULL)
		return false;

	//Compute the potential function
	compute3dPotential(imageDataPtr, potentialFuncPtr, length, width, height, seedPoints);

	//STEP 1
	//In labelAray we have : 1 ---> already processed, 2 ---> in process and 3 ---> not processed
	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				xd = x_new(j, i, width);
				distanceFuncPtr[k][xd] = BIG_VALUE;
				labelArray[k][xd] = 3;
			}
		}
	}
	//--------------------End of STEP 1 -----------------------------------

	//STEP 2
	//Treat the starting point
	i = seedPoints[0].y; 
	j = seedPoints[0].x; 
	k = seedPoints[0].z;

	size_t currentIndx = x_new(j, i, width);

	distanceFuncPtr[k][currentIndx] = 0.0;
	
	size_t kplus = 0, kminus = 0, iplus = 0, iminus = 0, jplus = 0, jminus = 0;
	size_t indxWest = 0, indxEast = 0, indxNorth = 0, indxSouth = 0;
	const size_t height_minus = height - 1, length_minus = length - 1, width_minus = width - 1;

	//STEP 3
	//Find the neigbors of the starting point

	//Top
	if (k > 0 && j >= 0 && j < width && i >= 0 && i < length) {
		kminus = k - 1;
		labelArray[kminus][x_new(j, i, width)] = 2;
		k_inProcess.push_back(kminus);
		j_inProcess.push_back(j);
		i_inProcess.push_back(i);
	}

	//Bottom
	if (k < height_minus && j >= 0 && j < width && i >= 0 && i < length) {
		kplus = k + 1;
		labelArray[kplus][x_new(j, i, width)] = 2;
		k_inProcess.push_back(kplus);
		j_inProcess.push_back(j);
		i_inProcess.push_back(i);
	}

	//North
	if (k >= 0 && k < height && j > 0 && i >= 0 && i < length) {
		jminus = j - 1;
		labelArray[k][x_new(jminus, i, width)] = 2;
		k_inProcess.push_back(k);
		j_inProcess.push_back(jminus);
		i_inProcess.push_back(i);
	}

	//South
	if (k >= 0 && k < height && j < width_minus && i >= 0 && i < length) {
		jplus = j + 1;
		labelArray[k][x_new(jplus, i, width)] = 2;
		k_inProcess.push_back(k);
		j_inProcess.push_back(jplus);
		i_inProcess.push_back(i);
	}

	//East
	if (k >= 0 && k < height && j >= 0 && j < width && i < length_minus) {
		iplus = i + 1;
		labelArray[k][x_new(j, iplus, width)] = 2;
		k_inProcess.push_back(k);
		j_inProcess.push_back(j);
		i_inProcess.push_back(iplus);
	}

	//West
	if (k >= 0 && k < height && j >= 0 && j < width && i > 0) {
		iminus = i - 1;
		labelArray[k][x_new(j, iminus, width)] = 2;
		k_inProcess.push_back(k);
		j_inProcess.push_back(j);
		i_inProcess.push_back(iminus);
	}

	//Compute the solution of neighbors in the stack
	lenghtSol = i_inProcess.size();
	dataType coefSpeed = 0.0;
	for (n = 0; n < lenghtSol; n++) {
		i = i_inProcess[n];
		j = j_inProcess[n];
		k = k_inProcess[n];
		x = select3dX(distanceFuncPtr, length, width, height, i, j, k);
		y = select3dY(distanceFuncPtr, length, width, height, i, j, k);
		z = select3dZ(distanceFuncPtr, length, width, height, i, j, k);
		coefSpeed = potentialFuncPtr[k][x_new(j, i, width)];
		dist = solve3dQuadratic(x, y, z, coefSpeed);
		tempTimeFunc.push_back(dist);
	}

	//---------------------End of STEP 2 -------------------------------------

	//STEP 3
	while (i_inProcess.size() != 0) {

		lenghtSol = i_inProcess.size();
		//kplus = 0, kminus = 0, iplus = 0, iminus = 0, jplus = 0, jminus = 0;
		//indxWest = 0, indxEast = 0, indxNorth = 0, indxSouth = 0;
		//dNorth = 0.0, dSouth = 0.0, dEast = 0.0, dWest = 0.0, dTop = 0.0, dBottom = 0.0;

		//Find the minimal solution
		minSolution = INFINITY;
		for (n = 0; n < lenghtSol; n++) {
			if (minSolution > tempTimeFunc[n]) {
				minSolution = tempTimeFunc[n];
				indxSol = n;
			}
		}

		i = i_inProcess[indxSol];
		j = j_inProcess[indxSol];
		k = k_inProcess[indxSol];
		currentIndx = x_new(j, i, width);

		//Set the distance to the processed pixel

		//if (distanceFuncPtr[k][currentIndx] > minSolution) {
		//	distanceFuncPtr[k][currentIndx] = minSolution;
		//}

		//cout << "min distance = " << minSolution << endl;
		 
		distanceFuncPtr[k][currentIndx] = minSolution;
		labelArray[k][currentIndx] = 1;

		//i_Processed.push_back(i); 
		//j_Processed.push_back(j);
		//k_Processed.push_back(k);

		//Remove the processed pixel from inProcess
		if (indxSol == 0 ) {
			i_inProcess.erase(i_inProcess.begin());
			j_inProcess.erase(j_inProcess.begin());
			k_inProcess.erase(k_inProcess.begin());
			tempTimeFunc.erase(tempTimeFunc.begin());
		}
		else {
			i_inProcess.erase(i_inProcess.begin() + indxSol);
			j_inProcess.erase(j_inProcess.begin() + indxSol);
			k_inProcess.erase(k_inProcess.begin() + indxSol);
			tempTimeFunc.erase(tempTimeFunc.begin() + indxSol);
		}


		//STEP 4
		//Find the neighbors of the processed pixel and compute they time

		//Bottom
		if (k < height_minus && i >= 0 && i < length && j >= 0 && j < width){
			
			kplus = k + 1;

			if (labelArray[kplus][currentIndx] != 1) {
				x = select3dX(distanceFuncPtr, length, width, height, i, j, kplus);
				y = select3dY(distanceFuncPtr, length, width, height, i, j, kplus);
				z = select3dZ(distanceFuncPtr, length, width, height, i, j, kplus);
				coefSpeed = potentialFuncPtr[kplus][currentIndx];
				dBottom = solve3dQuadratic(x, y, z, coefSpeed);
				if (labelArray[kplus][currentIndx] == 3) {
					distanceFuncPtr[kplus][currentIndx] = dBottom;
					i_inProcess.push_back(i);
					j_inProcess.push_back(j);
					k_inProcess.push_back(kplus);
					tempTimeFunc.push_back(dBottom);
					labelArray[kplus][currentIndx] = 2;
				}
				else {
					if (labelArray[kplus][currentIndx] == 2) {
						if (dBottom < distanceFuncPtr[kplus][currentIndx]) {
							distanceFuncPtr[kplus][currentIndx] = dBottom;
						}
					}
				}
			}
		}
		
		//Top
		if (k > 0 && i >= 0 && i < length && j >= 0 && j < width) {

			kminus = k - 1;

			if (labelArray[kminus][currentIndx] != 1) {
				x = select3dX(distanceFuncPtr, length, width, height, i, j, kminus);
				y = select3dY(distanceFuncPtr, length, width, height, i, j, kminus);
				z = select3dZ(distanceFuncPtr, length, width, height, i, j, kminus);
				coefSpeed = potentialFuncPtr[kminus][currentIndx];
				dTop = solve3dQuadratic(x, y, z, coefSpeed);
				if (labelArray[kminus][currentIndx] == 3) {
					distanceFuncPtr[kminus][currentIndx] = dTop;
					i_inProcess.push_back(i);
					j_inProcess.push_back(j);
					k_inProcess.push_back(kminus);
					tempTimeFunc.push_back(dTop);
					labelArray[kminus][currentIndx] = 2;
				}
				else {
					if (labelArray[kminus][currentIndx] == 2) {
						if (dTop < distanceFuncPtr[kminus][currentIndx]) {
							distanceFuncPtr[kminus][currentIndx] = dTop;
						}
					}
				}
			}

		}

		//East
		if (i < length_minus && j >= 0 && j < width && k >= 0 && k < height){

			iplus = i + 1;
			indxEast = x_new(j, iplus, width);

			if (labelArray[k][indxEast] != 1) {
				x = select3dX(distanceFuncPtr, length, width, height, iplus, j, k);
				y = select3dY(distanceFuncPtr, length, width, height, iplus, j, k);
				z = select3dZ(distanceFuncPtr, length, width, height, iplus, j, k);
				coefSpeed = potentialFuncPtr[k][indxEast];
				dEast = solve3dQuadratic(x, y, z, coefSpeed);
				if (labelArray[k][indxEast] == 3) {
					distanceFuncPtr[k][indxEast] = dEast;
					i_inProcess.push_back(iplus);
					j_inProcess.push_back(j);
					k_inProcess.push_back(k);
					tempTimeFunc.push_back(dEast);
					labelArray[k][indxEast] = 2;
				}
				else {
					if (labelArray[k][indxEast] == 2) {
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

			if (labelArray[k][indxWest] != 1) {
				x = select3dX(distanceFuncPtr, length, width, height, iminus, j, k);
				y = select3dY(distanceFuncPtr, length, width, height, iminus, j, k);
				z = select3dZ(distanceFuncPtr, length, width, height, iminus, j, k);
				coefSpeed = potentialFuncPtr[k][indxWest];
				dWest = solve3dQuadratic(x, y, z, coefSpeed);
				if (labelArray[k][indxWest] == 3) {
					distanceFuncPtr[k][indxWest] = dWest;
					i_inProcess.push_back(iminus);
					j_inProcess.push_back(j);
					k_inProcess.push_back(k);
					tempTimeFunc.push_back(dWest);
					labelArray[k][indxWest] = 2;
				}
				else {
					if (labelArray[k][indxWest] == 2) {
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

			if (labelArray[k][indxNorth] != 1) {
				x = select3dX(distanceFuncPtr, length, width, height, i, jminus, k);
				y = select3dY(distanceFuncPtr, length, width, height, i, jminus, k);
				z = select3dZ(distanceFuncPtr, length, width, height, i, jminus, k);
				coefSpeed = potentialFuncPtr[k][indxNorth];
				dNorth = solve3dQuadratic(x, y, z, coefSpeed);
				if (labelArray[k][indxNorth] == 3) {
					distanceFuncPtr[k][indxNorth] = dNorth;
					i_inProcess.push_back(i);
					j_inProcess.push_back(jminus);
					k_inProcess.push_back(k);
					tempTimeFunc.push_back(dNorth);
					labelArray[k][indxNorth] = 2;
				}
				else {
					if (labelArray[k][indxNorth] == 2) {
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

			if (labelArray[k][indxSouth] != 1) {
				x = select3dX(distanceFuncPtr, length, width, height, i, jplus, k);
				y = select3dY(distanceFuncPtr, length, width, height, i, jplus, k);
				z = select3dZ(distanceFuncPtr, length, width, height, i, jplus, k);
				coefSpeed = potentialFuncPtr[k][indxSouth];
				dSouth = solve3dQuadratic(x, y, z, coefSpeed);
				if (labelArray[k][indxSouth] == 3) {
					distanceFuncPtr[k][indxSouth] = dSouth;
					i_inProcess.push_back(i);
					j_inProcess.push_back(jplus);
					k_inProcess.push_back(k);
					tempTimeFunc.push_back(dSouth);
					labelArray[k][indxSouth] = 2;
				}
				else {
					if (labelArray[k][indxSouth] == 2) {
						if (dSouth < distanceFuncPtr[k][indxSouth]) {
							distanceFuncPtr[k][indxSouth] = dSouth;
						}
					}
				}
			}

		}

	}

	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				if (labelArray[k][x_new(j, i, width)] == 3) {
					cpt++;
				}
			}
		}
	}
	if (cpt == 0) {
		cout << "\nAll the pixels have been visited" << endl;
	}
	else {
		cout << "\n" << cpt << " pixels haven't been visited" << endl;
	}

	for (k = 0; k < height; k++) {
		delete[] labelArray[k];
	}
	delete[] labelArray;

	return true;
}

bool shortestPath3d(dataType** distanceFuncPtr, dataType** resultedPath, const size_t length, const size_t width, const size_t height, dataType h, Point3d* seedPoints) {

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

//------------------------------------------------------------------

void swapNeighbor(neighborPoint* a, neighborPoint* b) {
	neighborPoint temp = *a;
	*a = *b;
	*b = temp;
}

void heapify(vector<neighborPoint>& in_Process, int length_InProcess, int i) {

	int current = i;
	int left_child = 2 * i + 1;
	int right_child = 2 * i + 2;

	dataType val_current = 0.0, val_left = 0.0, val_right = 0.0;

	if (current >= 0 && current < length_InProcess) {
		val_current = in_Process[current].arrival;
	}
	
	if (left_child < length_InProcess) {
		val_left = in_Process[left_child].arrival;
		if (val_left < val_current) {
			current = left_child;
			val_current = in_Process[current].arrival;
		}
	}

	if (right_child < length_InProcess) {
		val_right = in_Process[right_child].arrival;
		if (val_right < val_current) {
			current = right_child;
		}
	}

	if (current != i) {
		swapNeighbor(&in_Process[i], &in_Process[current]);
		heapify(in_Process, length_InProcess, current);
	}

}

void createMinHeapStructure(vector<neighborPoint>& in_Process, int length_InProcess){
	int ind, start = length_InProcess / 2 - 1; 
	for (ind = start; ind >= 0; ind--) {
		heapify(in_Process, length_InProcess, ind);
	}
}

void heapifyBottomToUp(vector<neighborPoint>& in_Process, int i) {

	//int l = in_Process.size();
	int current = i;

	//if (i > 0 && i % 2 == 0) {
	//	int parent = (i - 1) / 2;
	//	int brother = i - 1;
	//	dataType val_current = in_Process[current].arrival;
	//	dataType val_parent = in_Process[parent].arrival;
	//	dataType val_brother = in_Process[brother].arrival;
	//	if (val_current < val_parent) {
	//		current = parent;
	//		val_current = in_Process[current].arrival;
	//	}
	//	if (val_brother < val_parent) {
	//		current = brother;
	//	}
	//}
	//else {
	//	int parent = (i - 1) / 2;
	//	dataType val_current = in_Process[current].arrival;
	//	dataType val_parent = in_Process[parent].arrival;
	//	if (val_current < val_parent) {
	//		current = parent;
	//	}
	//}

	if (i > 0) {
		int parent = (i - 1) / 2;
		dataType val_current = in_Process[current].arrival;
		dataType val_parent = in_Process[parent].arrival;
		if (val_current < val_parent) {
			current = parent;
		}
	}

	if (current != i) {
		swapNeighbor(&in_Process[current], &in_Process[i]);
		heapifyBottomToUp(in_Process, current);
	}
	
}

bool FastMarching3DNewVersion(dataType** imageDataPtr, dataType** distanceFuncPtr, dataType** potentialFuncPtr, const size_t length, const size_t width, const size_t height, Point3d* seedPoints) {

	if (imageDataPtr == NULL || distanceFuncPtr == NULL || potentialFuncPtr == NULL || seedPoints == NULL) {
		return false;
	}

	vector <neighborPoint> inProcess;
	size_t i = 0, j = 0, k = 0;

	size_t ** labelArray = new size_t*[height];
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
	neighborPoint current;
	j = seedPoints[0].x;
	i = seedPoints[0].y;
	k = seedPoints[0].z;
	currentIndx = x_new(j, i, width);
	//current.label = 1;
	//current.arrival = 0.0;
	//current = {j, i, k, 1, poxy, 0.0};
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
		neighborPoint TopNeighbor = { j, i, kminus, dTop };
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
		neighborPoint BottomNeighbor = { j, i, kplus, dBottom };
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
		neighborPoint NorthNeighbor = { jminus, i, k, dNorth };
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
		neighborPoint SouthNeighbor = { jplus, i, k, dSouth };
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
		neighborPoint EastNeighbor = { j, iplus, k, dEast };
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
		neighborPoint WestNeighbor = {j, iminus, k, dWest};
		distanceFuncPtr[k][indxWest] = dWest;
		inProcess.push_back(WestNeighbor);
		labelArray[k][indxWest] = 2;
	}

	int n = inProcess.size();

	createMinHeapStructure(inProcess, n);
	size_t label = 0;

	int l = 0, m = 0;
	 
	while (inProcess.size() != 0) {

		//processed the minimum
		current = inProcess[0]; // index 0 exist because, if not we will be out of the while loop

		//cout << "min distance find by heap =  " << current.arrival << endl;

		//dataType minS = INFINITY;
		//for (size_t t = 0; t < inProcess.size(); t++) {
		//	if (minS > inProcess[t].arrival) {
		//		minS = inProcess[t].arrival;
		//	}
		//}

		//cout << "min distance find by classical min search =  " << minS << endl;

		j = current.x;
		i = current.y;
		k = current.z;
		currentIndx = x_new(j, i, width);
		labelArray[k][currentIndx] = 1;
		distanceFuncPtr[k][currentIndx] = current.arrival;

		//if (distanceFuncPtr[k][poxy] > current.arrival) {
		//	distanceFuncPtr[k][poxy] = current.arrival;
		//}

		l = inProcess.size();
		//if (l > 1) {
		//	swapNeighbor(&inProcess[0], &inProcess[l - 1]);
		//}

		//dataType minSol = INFINITY;
		//for (int t = 0; i < inProcess.size(); t++) {
		//	if (minSol > inProcess[t].arrival) {
		//		minSol = inProcess[t].arrival;
		//	}
		//}

		swapNeighbor(&inProcess[0], &inProcess[l - 1]);
		inProcess.pop_back();
		m = inProcess.size();
		heapify(inProcess, m, 0);
		
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
					neighborPoint BottomNeighbor = { j, i, kplus, dBottom };
					inProcess.push_back(BottomNeighbor);
					l = inProcess.size();
					heapifyBottomToUp(inProcess, l - 1);
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
					neighborPoint TopNeighbor = { j, i, kminus, dTop };
					inProcess.push_back(TopNeighbor);
					l = inProcess.size();
					heapifyBottomToUp(inProcess, l - 1);
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
					neighborPoint EastNeighbor = { j, iplus, k, dEast };
					inProcess.push_back(EastNeighbor);
					l = inProcess.size();
					heapifyBottomToUp(inProcess, l - 1);
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
					neighborPoint WestNeighbor = { j, iminus, k, dWest };
					inProcess.push_back(WestNeighbor);
					l = inProcess.size();
					heapifyBottomToUp(inProcess, l - 1);
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
					neighborPoint NorthNeighbor = { jminus, i, k, dNorth };
					inProcess.push_back(NorthNeighbor);
					l = inProcess.size();
					heapifyBottomToUp(inProcess, l - 1);
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
					neighborPoint SouthNeighbor = { jplus, i, k, dSouth };
					inProcess.push_back(SouthNeighbor);
					l = inProcess.size();
					heapifyBottomToUp(inProcess, l - 1);
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

	size_t cpt1 = 0, cpt2 = 0,cpt3 = 0 ;
	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				currentIndx = x_new(j, i, width);
				if (distanceFuncPtr[k][currentIndx] == INFINITY) {
					cpt1++;
				}
				if (labelArray[k][currentIndx] == 2) {
					cpt2++;
				}
				if (labelArray[k][currentIndx] == 3) {
					cpt3++;
				}
			}
		}
	}
	cout << cpt1 << " point(s) with distance = INFINITY" << endl;;
	cout << cpt2 << " point(s) with label = 2" << endl;
	cout << cpt3 << " point(s) with label = 3" << endl;

	for (k = 0; k < height; k++) {
		delete[] labelArray[k];
	}
	delete[] labelArray;
	return true;
}


