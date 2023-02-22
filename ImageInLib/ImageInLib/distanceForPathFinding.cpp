#include <iostream>
#include <climits>
#include <crtdbg.h>
#include <corecrt_malloc.h>
#include<cmath>
#include <omp.h>
#include<vector>
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

bool computeImageGradient(dataType* imageDataPtr, dataType* gradientVectorX, dataType* gradientVectorY , const size_t height, const size_t width, dataType h) {
	
	if (imageDataPtr == NULL || gradientVectorX == NULL || gradientVectorY == NULL) {
		return false;
	}

	size_t i, j, xd;
	dataType ux = 0.0, uy = 0.0, norm_gradient = 0.0;

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

			norm_gradient = norm_gradient + ux * ux + uy * uy;
			//norm_gradient = ux * ux + uy * uy + 0.001;
			gradientVectorX[xd] = ux; // / norm_gradient;
			gradientVectorY[xd] = uy; // / norm_gradient;

		}
	}

	norm_gradient = sqrt(norm_gradient);

	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			xd = x_new(j, i, width);
			gradientVectorX[xd] = gradientVectorX[xd] / norm_gradient;
			gradientVectorY[xd] = gradientVectorY[xd] / norm_gradient;
		}
	}

	return true;
}

dataType computeImageNorm2d(dataType* imageDataPtr, const size_t height, const size_t width) {
	size_t i, j, xd;
	dataType norm_array = 0.0;

	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			xd = x_new(j, i, width);
			norm_array = norm_array + imageDataPtr[xd] * imageDataPtr[xd];
		}
	}
	return sqrt(norm_array);
}

bool computePotential(dataType* imageDataPtr, dataType* potentialFuncPtr, const size_t height, const size_t width, Point2D * seedPoints) {

	if (imageDataPtr == NULL || potentialFuncPtr == NULL || seedPoints == NULL)
		return false;

	size_t i, j;
	size_t i1 = seedPoints[0].y, j1 = seedPoints[0].x, i2 = seedPoints[1].y, j2 = seedPoints[1].x;

	dataType seedVal = (imageDataPtr[x_new(j1, i1, width)] + imageDataPtr[x_new(j2, i2, width)]) / 2;
	size_t currentIndx = 0;
	dataType epsylon = 0.00001;

	dataType* gradientVectorX = new dataType[height * width];
	dataType* gradientVectorY = new dataType[height * width];

	computeImageGradient(imageDataPtr, gradientVectorX, gradientVectorY, height, width, 1.0);

	size_t seedIndice = x_new(j1, i1, width);
	dataType ux0 = gradientVectorX[seedIndice], uy0 = gradientVectorY[seedIndice], ux = 0.0, uy = 0.0;

	//Computation of potential function
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			currentIndx = x_new(j, i, width);
			//potentialFuncPtr[currentIndx] = abs(seedVal - imageDataPtr[currentIndx]);
			ux = gradientVectorX[currentIndx]; uy = gradientVectorY[currentIndx];
			potentialFuncPtr[currentIndx] = 1 / (0.001 + sqrt(pow(ux - ux0, 2) + pow(uy - uy0, 2))) ;
		}
	}

	////find the max potential
	//dataType max_potential = -1;
	//for (i = 0; i < height; i++) {
	//	for (j = 0; j < width; j++) {
	//		currentIndx = x_new(j, i, width);
	//		if (potentialFuncPtr[currentIndx] > max_potential) {
	//			max_potential = potentialFuncPtr[currentIndx];
	//		}
	//	}
	//}

	////Normalization
	//for (i = 0; i < height; i++) {
	//	for (j = 0; j < width; j++) {
	//		currentIndx = x_new(j, i, width);
	//		potentialFuncPtr[currentIndx] = (dataType)(epsylon + potentialFuncPtr[currentIndx] / max_potential);
	//	}
	//}

	delete[] gradientVectorX;
	delete[] gradientVectorY;

	return true;
}

bool fastMarching2d(dataType* imageDataPtr, dataType* distanceFuncPtr, dataType* potentialFuncPtr, const size_t height, const size_t width, Point2D * seedPoints)
{
	//short* labelArray = (short*)malloc(height * width * sizeof(short));
	short * labelArray = new short[height * width];

	if (imageDataPtr == NULL || distanceFuncPtr == NULL || potentialFuncPtr == NULL || seedPoints == NULL || labelArray == NULL) {
		return false;
	}

	size_t i = 0, j = 0, k = 0, dim2D = height * width, cpt = 0;
	vector<size_t> i_Processed, i_inProcess;
	vector<size_t> j_Processed, j_inProcess;
	dataType x = 0.0, y = 0.0, minSolution = 0.0, coef = 0.0, dist = 0.0;
	size_t nbNeighborsFound = 0;
	size_t h_n = 0, iSol = 0, jSol = 0, iNew = 0, jNew = 0, iNew_minus = 0, iNew_plus = 0, jNew_minus = 0, jNew_plus = 0;
	vector<dataType> tempTimeFunc;
	dataType dNorth = 0.0, dSouth = 0.0, dEast = 0.0, dWest = 0.0;

	//cout << "\nNumber of pixels : " << height * width  << endl;

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
	//Add the neighbors of the seed point in the vector of pixels to be processed
	i = seedPoints[0].y; j = seedPoints[0].x;
	distanceFuncPtr[x_new(j, i, width)] = 0;
	i_Processed.push_back(i); j_Processed.push_back(j);

	size_t iminus = i - 1, iplus = i + 1, jminus = j - 1, jplus = j + 1;
	if (i == 0) {
		if (j == 0) {

			//East
			if (labelArray[x_new(jplus, i, width)] == 3) {
				i_inProcess.push_back(i); j_inProcess.push_back(jplus);
				labelArray[x_new(jplus, i, width)] = 2;
				nbNeighborsFound++;
			}
			//South
			if (labelArray[x_new(j, iplus, width)] == 3) {
				i_inProcess.push_back(iplus); j_inProcess.push_back(j);
				labelArray[x_new(j, iplus, width)] = 2;
				nbNeighborsFound++;
			}
		}
		else {
			if ( j == (width - 1) ) {
				//West
				if (labelArray[x_new(jminus, i, width)] == 3) {
					i_inProcess.push_back(i); j_inProcess.push_back(jminus);
					labelArray[x_new(jminus, i, width)] = 2;
					nbNeighborsFound++;
				}
				//South
				if (labelArray[x_new(j, iplus, width)] == 3) {
					i_inProcess.push_back(iplus); j_inProcess.push_back(j);
					labelArray[x_new(j, iplus, width)] = 2;
					nbNeighborsFound++;
				}
			}
			else {

				//East
				if (labelArray[x_new(jplus, i, width)] == 3) {
					i_inProcess.push_back(i); j_inProcess.push_back(jplus);
					labelArray[x_new(jplus, i, width)] = 2;
					nbNeighborsFound++;
				}
				//West
				if (labelArray[x_new(jminus, i, width)] == 3) {
					i_inProcess.push_back(i); j_inProcess.push_back(jminus);
					labelArray[x_new(jminus, i, width)] = 2;
					nbNeighborsFound++;
				}
				//South
				if (labelArray[x_new(j, iplus, width)] == 3) {
					i_inProcess.push_back(iplus); j_inProcess.push_back(j);
					labelArray[x_new(j, iplus, width)] = 2;
					nbNeighborsFound++;
				}

			}
		}
	}
	else {
		if (i == (height - 1)) {
			if (j == 0) {

				//East
				if (labelArray[x_new(jplus, i, width)] == 3) {
					i_inProcess.push_back(i); j_inProcess.push_back(jplus);
					labelArray[x_new(jplus, i, width)] = 2;
					nbNeighborsFound++;
				}
				//North
				if (labelArray[x_new(j, iminus, width)] == 3) {
					i_inProcess.push_back(iminus); j_inProcess.push_back(j);
					labelArray[x_new(j, iminus, width)] = 2;
					nbNeighborsFound++;
				}
				
			}
			else {
				if (j == (width - 1)) {

					//North
					if (labelArray[x_new(j, iminus, width)] == 3) {
						i_inProcess.push_back(iminus); j_inProcess.push_back(j);
						labelArray[x_new(j, iminus, width)] = 2;
						nbNeighborsFound++;
					}
					//West
					if (labelArray[x_new(jminus, i, width)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(jminus);
						labelArray[x_new(jminus, i, width)] = 2;
						nbNeighborsFound++;
					}

				}
				else {

					//East
					if (labelArray[x_new(jplus, i, width)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(jplus);
						labelArray[x_new(jplus, i, width)] = 2;
						nbNeighborsFound++;
					}
					//North
					if (labelArray[x_new(j, iminus, width)] == 3) {
						i_inProcess.push_back(iminus); j_inProcess.push_back(j);
						labelArray[x_new(j, iminus, width)] = 2;
						nbNeighborsFound++;
					}
					//West
					if (labelArray[x_new(jminus, i, width)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(jminus);
						labelArray[x_new(jminus, i, width)] = 2;
						nbNeighborsFound++;
					}

				}
			}
		}
		else {
			if (j == 0) {

				//East
				if (labelArray[x_new(jplus, i, width)] == 3) {
					i_inProcess.push_back(i); j_inProcess.push_back(jplus);
					labelArray[x_new(jplus, i, width)] = 2;
					nbNeighborsFound++;
				}
				//North
				if (labelArray[x_new(j, iminus, width)] == 3) {
					i_inProcess.push_back(iminus); j_inProcess.push_back(j);
					labelArray[x_new(j, iminus, width)] = 2;
					nbNeighborsFound++;
				}
				//South
				if (labelArray[x_new(j, iplus, width)] == 3) {
					i_inProcess.push_back(iplus); j_inProcess.push_back(j);
					labelArray[x_new(j, iplus, width)] = 2;
					nbNeighborsFound++;
				}
			}
			else {
				if (j == (width - 1) ) {

					//North
					if (labelArray[x_new(j, iminus, width)] == 3) {
						i_inProcess.push_back(iminus); j_inProcess.push_back(j);
						labelArray[x_new(j, iminus, width)] = 2;
						nbNeighborsFound++;
					}
					//West
					if (labelArray[x_new(jminus, i, width)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(jminus);
						labelArray[x_new(jminus, i, width)] = 2;
						nbNeighborsFound++;
					}
					//South
					if (labelArray[x_new(j, iplus, width)] == 3) {
						i_inProcess.push_back(iplus); j_inProcess.push_back(j);
						labelArray[x_new(j, iplus, width)] = 2;
						nbNeighborsFound++;
					}
				}
				else {

					//East
					if (labelArray[x_new(jplus, i, width)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(jplus);
						labelArray[x_new(jplus, i, width)] = 2;
						nbNeighborsFound++;
					}
					//North
					if (labelArray[x_new(j, iminus, width)] == 3) {
						i_inProcess.push_back(iminus); j_inProcess.push_back(j);
						labelArray[x_new(j, iminus, width)] = 2;
						nbNeighborsFound++;
					}
					//West
					if (labelArray[x_new(jminus, i, width)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(jminus);
						labelArray[x_new(jminus, i, width)] = 2;
						nbNeighborsFound++;
					}
					//South
					if (labelArray[x_new(j, iplus, width)] == 3) {
						i_inProcess.push_back(iplus); j_inProcess.push_back(j);
						labelArray[x_new(j, iplus, width)] = 2;
						nbNeighborsFound++;
					}
				}
			}
		}
	}

	//Compute the solution for neighbors in the stack
	h_n = i_inProcess.size(); // h_n = j_inProcess.size();

	for (k = 0; k < h_n; k++) {
		x = selectX(distanceFuncPtr, height, width, i_inProcess[k], j_inProcess[k]);
		y = selectY(distanceFuncPtr, height, width, i_inProcess[k], j_inProcess[k]);
		coef = (dataType)(1.0 / potentialFuncPtr[x_new(j_inProcess[k], i_inProcess[k], width)] );
		dist = solve2dQuadratic(x, y, coef);
		tempTimeFunc.push_back(dist);
	}

	//Update the label of seed point as processed
	labelArray[x_new(j, i, width)] = 1;

	//cout << "Number of Neigbors found : " << nbNeighborsFound << endl;
	//---------------------End of STEP 2 -------------------------------------

	//STEP 3
	while (i_inProcess.size() != 0 || j_inProcess.size() != 0 ) {

		h_n = i_inProcess.size();
		// i_inProcess and j_inProcess have the same size

		//Find the minimal solution
		minSolution = INFINITY;
		for (k = 0; k < h_n; k++) {
			if (minSolution > tempTimeFunc[k]) {
				minSolution = tempTimeFunc[k];
				iSol = k; iNew = i_inProcess[k];
				jSol = k; jNew = j_inProcess[k];
			}
		}

		//cout << "\ndistance  : " << minSolution << endl;
		
		//Set the distance to the processed pixel
		distanceFuncPtr[x_new(jNew, iNew, width)] = minSolution;
		labelArray[x_new(jNew, iNew, width)] = 1;
		i_Processed.push_back(iNew); j_Processed.push_back(jNew);

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
		if (iNew == 0) {
			if (jNew == 0) {

				jNew_plus = jNew + 1; 
				iNew_plus = iNew + 1;

				//East
				x = selectX(distanceFuncPtr, height, width, iNew, jNew_plus);
				y = selectY(distanceFuncPtr, height, width, iNew, jNew_plus);
				coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_plus, iNew, width)]);
				dEast = solve2dQuadratic(x, y, coef);

				//South
				x = selectX(distanceFuncPtr, height, width, iNew_plus, jNew);
				y = selectY(distanceFuncPtr, height, width, iNew_plus, jNew);
				coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_plus, width)]);
				dSouth = solve2dQuadratic(x, y, coef);

			}
			else {
				if (jNew == (width - 1) ) {

					iNew_plus = iNew + 1;
					jNew_minus = jNew - 1;

					//West
					x = selectX(distanceFuncPtr, height, width, iNew, jNew_minus);
					y = selectY(distanceFuncPtr, height, width, iNew, jNew_minus);
					coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_minus, iNew, width)]);
					dWest = solve2dQuadratic(x, y, coef);

					//South
					x = selectX(distanceFuncPtr, height, width, iNew_plus, jNew);
					y = selectY(distanceFuncPtr, height, width, iNew_plus, jNew);
					coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_plus, width)]);
					dSouth = solve2dQuadratic(x, y, coef);

				}
				else {

					iNew_plus = iNew + 1;
					jNew_plus = jNew + 1;
					jNew_minus = jNew - 1;

					//West
					x = selectX(distanceFuncPtr, height, width, iNew, jNew_minus);
					y = selectY(distanceFuncPtr, height, width, iNew, jNew_minus);
					coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_minus, iNew, width)]);
					dWest = solve2dQuadratic(x, y, coef);

					//East
					x = selectX(distanceFuncPtr, height, width, iNew, jNew_plus);
					y = selectY(distanceFuncPtr, height, width, iNew, jNew_plus);
					coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_plus, iNew, width)]);
					dEast = solve2dQuadratic(x, y, coef);

					//South
					x = selectX(distanceFuncPtr, height, width, iNew_plus, jNew);
					y = selectY(distanceFuncPtr, height, width, iNew_plus, jNew);
					coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_plus, width)]);
					dSouth = solve2dQuadratic(x, y, coef);

				}
			}
		}
		else {
			if (iNew == (height - 1) ) {
				if (jNew == 0) {

					iNew_minus = iNew - 1;
					jNew_plus = jNew + 1;

					//North
					x = selectX(distanceFuncPtr, height, width, iNew_minus, jNew);
					y = selectY(distanceFuncPtr, height, width, iNew_minus, jNew);
					coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_minus, width)]);
					dNorth = solve2dQuadratic(x, y, coef);

					//East
					x = selectX(distanceFuncPtr, height, width, iNew, jNew_plus);
					y = selectY(distanceFuncPtr, height, width, iNew, jNew_plus);
					coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_plus, iNew, width)]);
					dEast = solve2dQuadratic(x, y, coef);

				}
				else {
					if (jNew == (width - 1) ) {

						iNew_minus = iNew - 1;
						jNew_minus = jNew - 1;

						//North
						x = selectX(distanceFuncPtr, height, width, iNew_minus, jNew);
						y = selectY(distanceFuncPtr, height, width, iNew_minus, jNew);
						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_minus, width)]);
						dNorth = solve2dQuadratic(x, y, coef);

						//West
						x = selectX(distanceFuncPtr, height, width, iNew, jNew_minus);
						y = selectY(distanceFuncPtr, height, width, iNew, jNew_minus);
						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_minus, iNew, width)]);
						dWest = solve2dQuadratic(x, y, coef);

					}
					else {

						iNew_minus = iNew - 1;
						jNew_plus = jNew + 1;
						jNew_minus = jNew - 1;

						//West
						x = selectX(distanceFuncPtr, height, width, iNew, jNew_minus);
						y = selectY(distanceFuncPtr, height, width, iNew, jNew_minus);
						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_minus, iNew, width)]);
						dWest = solve2dQuadratic(x, y, coef);

						//North
						x = selectX(distanceFuncPtr, height, width, iNew_minus, jNew);
						y = selectY(distanceFuncPtr, height, width, iNew_minus, jNew);
						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_minus, width)]);
						dNorth = solve2dQuadratic(x, y, coef);

						//East
						x = selectX(distanceFuncPtr, height, width, iNew, jNew_plus);
						y = selectY(distanceFuncPtr, height, width, iNew, jNew_plus);
						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_plus, iNew, width)]);
						dEast = solve2dQuadratic(x, y, coef);

					}
				}
			}
			else {
				if (jNew == 0) {

					iNew_minus = iNew - 1;
					iNew_plus = iNew + 1;
					jNew_plus = jNew + 1;

					//North
					x = selectX(distanceFuncPtr, height, width, iNew_minus, jNew);
					y = selectY(distanceFuncPtr, height, width, iNew_minus, jNew);
					coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_minus, width)]);
					dNorth = solve2dQuadratic(x, y, coef);

					//East
					x = selectX(distanceFuncPtr, height, width, iNew, jNew_plus);
					y = selectY(distanceFuncPtr, height, width, iNew, jNew_plus);
					coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_plus, iNew, width)]);
					dEast = solve2dQuadratic(x, y, coef);

					//South
					x = selectX(distanceFuncPtr, height, width, iNew_plus, jNew);
					y = selectY(distanceFuncPtr, height, width, iNew_plus, jNew);
					coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_plus, width)]);
					dSouth = solve2dQuadratic(x, y, coef);

				}
				else {
					if (jNew == (width - 1) ) {

						iNew_minus = iNew - 1;
						iNew_plus = iNew + 1;
						jNew_minus = jNew - 1;

						//North
						x = selectX(distanceFuncPtr, height, width, iNew_minus, jNew);
						y = selectY(distanceFuncPtr, height, width, iNew_minus, jNew);
						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_minus, width)]);
						dNorth = solve2dQuadratic(x, y, coef);

						//West
						x = selectX(distanceFuncPtr, height, width, iNew, jNew_minus);
						y = selectY(distanceFuncPtr, height, width, iNew, jNew_minus);
						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_minus, iNew, width)]);
						dWest = solve2dQuadratic(x, y, coef);

						//South
						x = selectX(distanceFuncPtr, height, width, iNew_plus, jNew);
						y = selectY(distanceFuncPtr, height, width, iNew_plus, jNew);
						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_plus, width)]);
						dSouth = solve2dQuadratic(x, y, coef);

					}
					else {

						iNew_minus = iNew - 1;
						iNew_plus = iNew + 1;
						jNew_plus = jNew + 1;
						jNew_minus = jNew - 1;

						//North
						x = selectX(distanceFuncPtr, height, width, iNew_minus, jNew);
						y = selectY(distanceFuncPtr, height, width, iNew_minus, jNew);
						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_minus, width)]);
						dNorth = solve2dQuadratic(x, y, coef);

						//West
						x = selectX(distanceFuncPtr, height, width, iNew, jNew_minus);
						y = selectY(distanceFuncPtr, height, width, iNew, jNew_minus);
						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_minus, iNew, width)]);
						dWest = solve2dQuadratic(x, y, coef);

						//East
						x = selectX(distanceFuncPtr, height, width, iNew, jNew_plus);
						y = selectY(distanceFuncPtr, height, width, iNew, jNew_plus);
						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_plus, iNew, width)]);
						dEast = solve2dQuadratic(x, y, coef);

						//South
						x = selectX(distanceFuncPtr, height, width, iNew_plus, jNew);
						y = selectY(distanceFuncPtr, height, width, iNew_plus, jNew);
						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_plus, width)]);
						dSouth = solve2dQuadratic(x, y, coef);
					}
				}
			}
		}

		//Test and update of neighbors
		//North
		if (labelArray[x_new(jNew, iNew_minus, width)] == 3) {
			distanceFuncPtr[x_new(jNew, iNew_minus, width)] = dNorth;
			i_inProcess.push_back(iNew_minus); j_inProcess.push_back(jNew);
			tempTimeFunc.push_back(dNorth);
			labelArray[x_new(jNew, iNew_minus, width)] = 2;
		}
		else {
			if (labelArray[x_new(jNew, iNew_minus, width)] == 2) {
				if (dNorth < distanceFuncPtr[x_new(jNew, iNew_minus, width)]) {
					distanceFuncPtr[x_new(jNew, iNew_minus, width)] = dNorth;
				}
			}
		}

		//West
		if (labelArray[x_new(jNew_minus, iNew, width)] == 3) {
			distanceFuncPtr[x_new(jNew_minus, iNew, width)] = dWest;
			i_inProcess.push_back(iNew); j_inProcess.push_back(jNew_minus);
			tempTimeFunc.push_back(dWest);
			labelArray[x_new(jNew_minus, iNew, width)] = 2;
		}
		else {
			if (labelArray[x_new(jNew_minus, iNew, width)] == 2) {
				if (dWest < distanceFuncPtr[x_new(jNew_minus, iNew, width)]) {
					distanceFuncPtr[x_new(jNew_minus, iNew, width)] = dWest;
				}
			}
		}

		//East
		if (labelArray[x_new(jNew_plus, iNew, width)] == 3) {
			distanceFuncPtr[x_new(jNew_plus, iNew, width)] = dEast;
			i_inProcess.push_back(iNew); j_inProcess.push_back(jNew_plus);
			tempTimeFunc.push_back(dEast);
			labelArray[x_new(jNew_plus, iNew, width)] = 2;
		}
		else {
			if (labelArray[x_new(jNew_plus, iNew, width)] == 2) {
				if (dEast < distanceFuncPtr[x_new(jNew_plus, iNew, width)]) {
					distanceFuncPtr[x_new(jNew_plus, iNew, width)] = dEast;
				}
			}
		}

		//South
		if (labelArray[x_new(jNew, iNew_plus, width)] == 3) {
			distanceFuncPtr[x_new(jNew, iNew_plus, width)] = dSouth;
			i_inProcess.push_back(iNew_plus); j_inProcess.push_back(jNew);
			tempTimeFunc.push_back(dSouth);
			labelArray[x_new(jNew, iNew_plus, width)] = 2;
		}
		else {
			if (labelArray[x_new(jNew, iNew_plus, width)] == 2) {
				if (dSouth < distanceFuncPtr[x_new(jNew, iNew_plus, width)]) {
					distanceFuncPtr[x_new(jNew, iNew_plus, width)] = dSouth;
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
		cout << "\nAll the pixels have been visited" << endl;
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
	dataType tau = 0.1, dist_min = INFINITY, tol = 1.0;
	size_t i_init = seedPoints[0].y, j_init = seedPoints[0].x, i_end = seedPoints[1].y, j_end = seedPoints[1].x;
	size_t i_current, j_current;
	dataType iNew = 0.0;
	dataType jNew = 0.0;

	dataType * gradientVectorX = new dataType[dim2d];
	dataType * gradientVectorY = new dataType[dim2d];

	computeImageGradient(distanceFuncPtr, gradientVectorX, gradientVectorY, height, width, 1.0);

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
	while(dist_min > tol && cpt < 1000000);

	cout << "\nDistance to the end point : " << dist_min << endl;
	cout << "\nNumber of iterations : " << cpt << endl;

	delete[] gradientVectorX;
	delete[] gradientVectorY;

	return true;
}

bool findPath2d(dataType* distanceFuncPtr, dataType* resultedPath, const size_t height, const size_t width, Point2D* seedPoints) {

	if (distanceFuncPtr == NULL || resultedPath == NULL || seedPoints == NULL)
		return false;

	size_t k, hn, coordWest, coordEast, coordNorth, coordSouth, coordCurrent, cpt = 0;
	size_t i_init = seedPoints[0].y, j_init = seedPoints[0].x, i_end = seedPoints[1].y, j_end = seedPoints[1].x;
	dataType min_value, pathLenght = sqrt((i_end - i_init) * (i_end - i_init) + (j_end - j_init) * (j_end - j_init)), tol = 1.0;
	vector<dataType> pointAction;
	vector<size_t> Indice_i, Indice_j;

	short * labelArray = new short[height * width];
	for (k = 0; k < height * width; k++)
		labelArray[k] = 0;

	coordCurrent = x_new(j_init, i_init, width); 
	labelArray[coordCurrent] = 1;

	while (pathLenght > tol && cpt < 1000000) {

		coordNorth = x_new(j_init, i_init - 1, width);
		coordSouth = x_new(j_init, i_init + 1, width);
		coordWest = x_new(j_init - 1, i_init, width);
		coordEast = x_new(j_init + 1, i_init - 1, width);

		if (labelArray[coordNorth] != 1) {
			pointAction.push_back(distanceFuncPtr[coordNorth]); Indice_j.push_back(j_init); Indice_i.push_back(i_init - 1);
		}
		if (labelArray[coordSouth] != 1) {
			pointAction.push_back(distanceFuncPtr[coordSouth]); Indice_j.push_back(j_init); Indice_i.push_back(i_init + 1);
		}
		if (labelArray[coordWest] != 1) {
			pointAction.push_back(distanceFuncPtr[coordWest]); Indice_j.push_back(j_init - 1); Indice_i.push_back(i_init);
		}
		if (labelArray[coordEast] != 1) {
			pointAction.push_back(distanceFuncPtr[coordEast]); Indice_j.push_back(j_init + 1); Indice_i.push_back(i_init);
		}
		
		hn = pointAction.size();
		min_value = INFINITY;
		for (k = 0; k < hn; k++) {
			if (pointAction[k] < min_value && pointAction[k]) {
				min_value = pointAction[k];
				j_init = Indice_j[k]; i_init = Indice_i[k];
			}
		}
		coordCurrent = x_new(j_init, i_init, width);
		resultedPath[coordCurrent] = 1;
		labelArray[coordCurrent] = 1;
		pathLenght = sqrt((i_end - i_init) * (i_end - i_init) + (j_end - j_init) * (j_end - j_init));
		cpt++;

		pointAction.clear(); Indice_i.clear(); Indice_j.clear();
	}
	return true;
}

//bool bruteForce2d(dataType* imageDataPtr, dataType* distanceFuncPtr, const size_t height, const size_t width, dataType backGround) {
//
//	if (imageDataPtr == NULL || distanceFuncPtr == NULL) {
//		return false;
//	}
//
//	dataType dist = 0.0, minDist = 0.0;
//	size_t i1 = 0, j1 = 0, x1 = 0, i2 = 0, j2 = 0, x2 = 0;
//
//	for (i1 = 0; i1 < height; i1++) {
//		for (j1 = 0; j1 < width; j2++) {
//
//			x1 = x_new(i1, j1, height);
//			minDist = INFINITY;
//
//			if (imageDataPtr[x1] != backGround) {
//
//				for (i2 = 0; i2 < height; i2++) {
//					for (j2 = 0; j2 < width; j2++) {
//
//						x2 = x_new(i2, j2, height);
//
//						if (imageDataPtr[x2] != backGround) {
//							dist = sqrt(pow(i1 - i2, 2) + pow(j1 - j2, 2));
//							if (dist != 0 && minDist > dist) {
//								minDist = dist;
//								distanceFuncPtr[x1] = minDist;
//							}
//						}
//						else {
//							distanceFuncPtr[x1] = 0;
//						}
//					}
//				}
//			}
//			else {
//				distanceFuncPtr[x1] = 0;
//			}
//		}
//	}
//	return true;
//}

//Functions for 3D images
// aU^2 - 2U(X+Y+Z) + (X^2 + Y^2 + Z^2 - W) = 0
dataType solve3dQuadratic(dataType X, dataType Y, dataType Z, dataType W) {

	dataType sol = 0.0, a, b, c, delta;

	a = 3;
	if (X == BIG_VALUE) {
		X = 0; a--;
	}
	if (Y == BIG_VALUE) {
		Y = 0; a--;
	}
	if (Z == BIG_VALUE) {
		Z = 0; a--;
	}

	b = -2 * (X + Y + Z); c = X * X + Y * Y + Z * Z - W;
	delta = b * b  - 4 * a * c;

	if (delta >= 0) {
		sol = (-b + sqrt(delta)) / (2 * a);
	}
	else {
		sol = min(X, min(Y, Z)) + W;
	}

	if (sol < 0) {
		cout << "The solution is negative " << endl; //If everything is OK, it never happen
		return 0;
	}
	else {
		return sol;
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

//bool compute3dPotential(dataType** imageDataPtr, dataType** potentialFuncPtr, const size_t length, const size_t width, const size_t height, Point3d* seedPoints) {
//
//	if (imageDataPtr == NULL || potentialFuncPtr == NULL || seedPoints == NULL)
//		return false;
//
//	size_t i, j, k, in = seedPoints->y, jn = seedPoints->x, kn =  seedPoints->z;
//
//	dataType seedVal = imageDataPtr[k][x_new(jn, in, width)];
//	size_t currentIndx = 0;
//	dataType epsylon = 0.01;
//
//	//Computation of potential function
//	for (k = 0; k < height; k++) {
//		for (i = 0; i < length; j++) {
//			for (j = 0; j < width; j++) {
//				currentIndx = x_new(j, i, width);
//				potentialFuncPtr[k][currentIndx] = abs(seedVal - imageDataPtr[k][currentIndx]);
//			}
//		}
//	}
//
//	//find the max potential
//	dataType max_potential = -1;
//	for (k = 0; k < height; k++) {
//		for (i = 0; i < length; i++) {
//			for (j = 0; j < width; j++) {
//				currentIndx = x_new(j, i, width);
//				if (potentialFuncPtr[k][currentIndx] > max_potential) {
//					max_potential = potentialFuncPtr[k][currentIndx];
//				}
//			}
//		}
//	}
//
//	//Normalization
//	for (k = 0; k < height; k++) {
//		for (i = 0; i < length; i++) {
//			for (j = 0; j < width; j++) {
//				currentIndx = x_new(j, i, width);
//				potentialFuncPtr[k][currentIndx] = (dataType)(epsylon + potentialFuncPtr[k][currentIndx] / max_potential);
//			}
//		}
//	}
//
//	return true;
//}

//bool fastMarching3d_N(dataType** imageDataPtr, dataType** distanceFuncPtr, dataType** potentialFuncPtr, const size_t length, const size_t height, const size_t width, Point3d* seedPoints)
//{
//	
//	if (imageDataPtr == NULL || distanceFuncPtr == NULL || potentialFuncPtr == NULL || seedPoints == NULL) {
//		return false;
//	}
//
//	size_t i = 0, j = 0, k = 0, xd, dim2D = length * width, cpt = 0;
//	vector<size_t> i_Processed, i_inProcess;
//	vector<size_t> j_Processed, j_inProcess;
//	vector<size_t> k_Processed, k_inProcess;
//	dataType x = 0.0, y = 0.0, minSolution = 0.0, coef = 0.0, dist = 0.0;
//	size_t nbNeighborsFound = 0;
//	size_t h_n = 0, iSol = 0, jSol = 0, kSol = 0;
//	size_t iNew = 0, jNew = 0, iNew_minus = 0, iNew_plus = 0, jNew_minus = 0, jNew_plus = 0;
//	vector<dataType> tempTimeFunc;
//	dataType dNorth = 0.0, dSouth = 0.0, dEast = 0.0, dWest = 0.0, dTop = 0.0, dBottom = 0.0;
//
//	short** labelArray = (short**)malloc(height * sizeof(short*));
//	for (k = 0; k < height; k++) {
//		labelArray[k] = (short*)malloc(dim2D * sizeof(short));
//	}
//
//	//Compute the potential function
//	compute3dPotential(imageDataPtr, potentialFuncPtr, length, width, height, seedPoints);
//
//	//STEP 1
//	//In labelAray we have : 1 ---> already processed, 2 ---> in process and 3 ---> not processed
//	for (k = 0; k < height; k++) {
//		for (i = 0; i < length; i++) {
//			for (j = 0; j < width; j++) {
//				xd = x_new(j, i, width);
//				distanceFuncPtr[k][xd] = BIG_VALUE;
//				labelArray[k][xd] = 3;
//			}
//		}
//	}
//	//--------------------End of STEP 1 -----------------------------------
//
//	//STEP 2
//	//Add the starting point to the list
//	i = seedPoints[0].y; 
//	j = seedPoints[0].x; 
//	k = seedPoints[0].z;
//	distanceFuncPtr[k][x_new(j, i, width)] = 0;
//	i_Processed.push_back(i); j_Processed.push_back(j); k_Processed.push_back(k);
//	tempTimeFunc.push_back(0);
//
//	////Compute the solution for neighbors in the stack
//	//h_n = i_inProcess.size(); // h_n = j_inProcess.size();
//	//dataType coefSpeed = 1.0;
//	//for (k = 0; k < h_n; k++) {
//	//	x = selectX(distanceFuncPtr, height, width, i_inProcess[k], j_inProcess[k]);
//	//	y = selectY(distanceFuncPtr, height, width, i_inProcess[k], j_inProcess[k]);
//	//	coef = (dataType)(1.0 / potentialFuncPtr[x_new(j_inProcess[k], i_inProcess[k], width)]);
//	//	dist = solve2dQuadratic(x, y, coef);
//	//	tempTimeFunc.push_back(dist);
//	//}
//
//	//---------------------End of STEP 2 -------------------------------------
//
//	//STEP 3
//	while (tempTimeFunc.size() != 0) {
//
//		//h_n = i_inProcess.size();
//		if (k == 0) {
//
//			if (i == 0) {
//
//				if (j == 0) {
//
//					//East
//					if (labelArray[k][x_new(j, i + 1, width)] == 3) {
//						k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i + 1);
//						labelArray[k][x_new(j, i + 1, width)] = 2;
//					}
//
//					//South
//					if (labelArray[k][x_new(j + 1, i, width)] == 3) {
//						k_inProcess.push_back(k); j_inProcess.push_back(j + 1); i_inProcess.push_back(i);
//						labelArray[k][x_new(j + 1, i, width)] = 2;
//					}
//
//					//Bottom
//					if (labelArray[k + 1][x_new(j, i, width)] == 3) {
//						k_inProcess.push_back(k + 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//						labelArray[k + 1][x_new(j, i, width)] = 2;
//					}
//
//				}
//				else {
//					if (j == width - 1) {
//
//						//West
//						if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//							k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//							labelArray[k][x_new(j, i - 1, width)] = 2;
//						}
//
//						//South
//						if (labelArray[k][x_new(j + 1, i, width)] == 3) {
//							k_inProcess.push_back(k); j_inProcess.push_back(j + 1); i_inProcess.push_back(i);
//							labelArray[k][x_new(j + 1, i, width)] = 2;
//						}
//
//						//Bottom
//						if (labelArray[k + 1][x_new(j, i, width)] == 3) {
//							k_inProcess.push_back(k + 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//							labelArray[k + 1][x_new(j, i, width)] = 2;
//						}
//
//					}
//					else {
//
//						//East
//						if (labelArray[k][x_new(j, i + 1, width)] == 3) {
//							k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i + 1);
//							labelArray[k][x_new(j, i + 1, width)] = 2;
//						}
//
//						//West
//						if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//							k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//							labelArray[k][x_new(j, i - 1, width)] = 2;
//						}
//
//						//South
//						if (labelArray[k][x_new(j + 1, i, width)] == 3) {
//							k_inProcess.push_back(k); j_inProcess.push_back(j + 1); i_inProcess.push_back(i);
//							labelArray[k][x_new(j + 1, i, width)] = 2;
//						}
//
//						//Bottom
//						if (labelArray[k + 1][x_new(j, i, width)] == 3) {
//							k_inProcess.push_back(k + 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//							labelArray[k + 1][x_new(j, i, width)] = 2;
//						}
//					}
//				}
//			}
//			else {
//				if (i == length - 1) {
//
//					if (j == 0) {
//
//						//East
//						if (labelArray[k][x_new(j, i + 1, width)] == 3) {
//							k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i + 1);
//							labelArray[k][x_new(j, i + 1, width)] = 2;
//						}
//
//						//South
//						if (labelArray[k][x_new(j + 1, i, width)] == 3) {
//							k_inProcess.push_back(k); j_inProcess.push_back(j + 1); i_inProcess.push_back(i);
//							labelArray[k][x_new(j + 1, i, width)] = 2;
//						}
//
//						//Bottom
//						if (labelArray[k + 1][x_new(j, i, width)] == 3) {
//							k_inProcess.push_back(k + 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//							labelArray[k + 1][x_new(j, i, width)] = 2;
//						}
//
//					}
//					else {
//						if (j == width - 1) {
//
//							//West
//							if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//								labelArray[k][x_new(j, i - 1, width)] = 2;
//							}
//
//							//North
//							if (labelArray[k][x_new(j - 1, i, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j - 1); i_inProcess.push_back(i);
//								labelArray[k][x_new(j + 1, i, width)] = 2;
//							}
//
//							//Bottom
//							if (labelArray[k + 1][x_new(j, i, width)] == 3) {
//								k_inProcess.push_back(k + 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//								labelArray[k + 1][x_new(j, i, width)] = 2;
//							}
//
//						}
//						else {
//							//West
//							if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//								labelArray[k][x_new(j, i - 1, width)] = 2;
//							}
//
//							//North
//							if (labelArray[k][x_new(j - 1, i, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j - 1); i_inProcess.push_back(i);
//								labelArray[k][x_new(j - 1, i, width)] = 2;
//							}
//
//							//South
//							if (labelArray[k][x_new(j + 1, i, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j + 1); i_inProcess.push_back(i);
//								labelArray[k][x_new(j + 1, i, width)] = 2;
//							}
//
//							//Bottom
//							if (labelArray[k + 1][x_new(j, i, width)] == 3) {
//								k_inProcess.push_back(k + 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//								labelArray[k + 1][x_new(j, i, width)] = 2;
//							}
//						}
//					}
//				}
//				else {
//					if (j == 0) {
//
//						//West
//						if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//							k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//							labelArray[k][x_new(j, i - 1, width)] = 2;
//						}
//
//						//East
//						if (labelArray[k][x_new(j, i + 1, width)] == 3) {
//							k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i + 1);
//							labelArray[k][x_new(j, i + 1, width)] = 2;
//						}
//
//						//South
//						if (labelArray[k][x_new(j + 1, i, width)] == 3) {
//							k_inProcess.push_back(k); j_inProcess.push_back(j + 1); i_inProcess.push_back(i);
//							labelArray[k][x_new(j + 1, i, width)] = 2;
//						}
//
//						//Bottom
//						if (labelArray[k + 1][x_new(j, i, width)] == 3) {
//							k_inProcess.push_back(k + 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//							labelArray[k + 1][x_new(j, i, width)] = 2;
//						}
//
//					}
//					else {
//						if (j == width - 1) {
//
//							//West
//							if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//								labelArray[k][x_new(j, i - 1, width)] = 2;
//							}
//
//							//East
//							if (labelArray[k][x_new(j, i + 1, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i + 1);
//								labelArray[k][x_new(j, i + 1, width)] = 2;
//							}
//
//							//North
//							if (labelArray[k][x_new(j - 1, i, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j - 1); i_inProcess.push_back(i);
//								labelArray[k][x_new(j - 1, i, width)] = 2;
//							}
//
//							//Bottom
//							if (labelArray[k + 1][x_new(j, i, width)] == 3) {
//								k_inProcess.push_back(k + 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//								labelArray[k + 1][x_new(j, i, width)] = 2;
//							}
//
//						}
//						else {
//
//							//West
//							if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//								labelArray[k][x_new(j, i - 1, width)] = 2;
//							}
//
//							//East
//							if (labelArray[k][x_new(j, i + 1, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i + 1);
//								labelArray[k][x_new(j, i + 1, width)] = 2;
//							}
//
//							//North
//							if (labelArray[k][x_new(j - 1, i, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j - 1); i_inProcess.push_back(i);
//								labelArray[k][x_new(j - 1, i, width)] = 2;
//							}
//
//							//South
//							if (labelArray[k][x_new(j + 1, i, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j + 1); i_inProcess.push_back(i);
//								labelArray[k][x_new(j + 1, i, width)] = 2;
//							}
//
//							//Bottom
//							if (labelArray[k + 1][x_new(j, i, width)] == 3) {
//								k_inProcess.push_back(k + 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//								labelArray[k + 1][x_new(j, i, width)] = 2;
//							}
//
//						}
//					}
//				}
//			}
//		}
//		else {
//			if (k == height - 1) {
//
//				if (i == 0) {
//
//					if (j == 0) {
//
//						//East
//						if (labelArray[k][x_new(j, i + 1, width)] == 3) {
//							k_inProcess.push_back(k + 1); j_inProcess.push_back(j); i_inProcess.push_back(i + 1);
//							labelArray[k][x_new(j, i + 1, width)] = 2;
//						}
//
//						//South
//						if (labelArray[k][x_new(j + 1, i, width)] == 3) {
//							k_inProcess.push_back(k); j_inProcess.push_back(j + 1); i_inProcess.push_back(i);
//							labelArray[k][x_new(j + 1, i, width)] = 2;
//						}
//
//						//Top
//						if (labelArray[k - 1][x_new(j, i, width)] == 3) {
//							k_inProcess.push_back(k - 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//							labelArray[k - 1][x_new(j, i, width)] = 2;
//						}
//
//					}
//					else {
//						if (j == width - 1) {
//
//							//West
//							if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//								labelArray[k][x_new(j, i + 1, width)] = 2;
//							}
//
//							//South
//							if (labelArray[k][x_new(j + 1, i, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j + 1); i_inProcess.push_back(i);
//								labelArray[k][x_new(j + 1, i, width)] = 2;
//							}
//
//							//Top
//							if (labelArray[k - 1][x_new(j, i, width)] == 3) {
//								k_inProcess.push_back(k - 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//								labelArray[k - 1][x_new(j, i, width)] = 2;
//							}
//
//						}
//						else {
//
//							//East
//							if (labelArray[k][x_new(j, i + 1, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i + 1);
//								labelArray[k][x_new(j, i + 1, width)] = 2;
//							}
//
//							//West
//							if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//								labelArray[k][x_new(j, i + 1, width)] = 2;
//							}
//
//							//South
//							if (labelArray[k][x_new(j + 1, i, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j + 1); i_inProcess.push_back(i);
//								labelArray[k][x_new(j + 1, i, width)] = 2;
//							}
//
//							//Top
//							if (labelArray[k - 1][x_new(j, i, width)] == 3) {
//								k_inProcess.push_back(k - 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//								labelArray[k - 1][x_new(j, i, width)] = 2;
//							}
//						}
//					}
//				}
//				else {
//					if (i == length - 1) {
//
//						if (j == 0) {
//
//							//West
//							if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//								labelArray[k][x_new(j, i - 1, width)] = 2;
//							}
//
//							//South
//							if (labelArray[k][x_new(j + 1, i, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j + 1); i_inProcess.push_back(i);
//								labelArray[k][x_new(j + 1, i, width)] = 2;
//							}
//
//							//Top
//							if (labelArray[k - 1][x_new(j, i, width)] == 3) {
//								k_inProcess.push_back(k - 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//								labelArray[k - 1][x_new(j, i, width)] = 2;
//							}
//
//						}
//						else {
//							if (j == width - 1) {
//
//								//West
//								if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//									labelArray[k][x_new(j, i - 1, width)] = 2;
//								}
//
//								//North
//								if (labelArray[k][x_new(j - 1, i, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j - 1); i_inProcess.push_back(i);
//									labelArray[k][x_new(j + 1, i, width)] = 2;
//								}
//
//								//Top
//								if (labelArray[k - 1][x_new(j, i, width)] == 3) {
//									k_inProcess.push_back(k - 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//									labelArray[k - 1][x_new(j, i, width)] = 2;
//								}
//
//							}
//							else {
//								//West
//								if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//									labelArray[k][x_new(j, i - 1, width)] = 2;
//								}
//
//								//North
//								if (labelArray[k][x_new(j - 1, i, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j - 1); i_inProcess.push_back(i);
//									labelArray[k][x_new(j - 1, i, width)] = 2;
//								}
//
//								//South
//								if (labelArray[k][x_new(j + 1, i, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j + 1); i_inProcess.push_back(i);
//									labelArray[k][x_new(j + 1, i, width)] = 2;
//								}
//
//								//Top
//								if (labelArray[k - 1][x_new(j, i, width)] == 3) {
//									k_inProcess.push_back(k - 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//									labelArray[k - 1][x_new(j, i, width)] = 2;
//								}
//							}
//						}
//					}
//					else {
//						if (j == 0) {
//
//							//West
//							if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//								labelArray[k][x_new(j, i - 1, width)] = 2;
//							}
//
//							//East
//							if (labelArray[k][x_new(j, i + 1, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i + 1);
//								labelArray[k][x_new(j, i + 1, width)] = 2;
//							}
//
//							//South
//							if (labelArray[k][x_new(j + 1, i, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j + 1); i_inProcess.push_back(i);
//								labelArray[k][x_new(j + 1, i, width)] = 2;
//							}
//
//							//Top
//							if (labelArray[k - 1][x_new(j, i, width)] == 3) {
//								k_inProcess.push_back(k - 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//								labelArray[k - 1][x_new(j, i, width)] = 2;
//							}
//
//						}
//						else {
//							if (j == width - 1) {
//
//								//West
//								if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//									labelArray[k][x_new(j, i - 1, width)] = 2;
//								}
//
//								//East
//								if (labelArray[k][x_new(j, i + 1, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i + 1);
//									labelArray[k][x_new(j, i + 1, width)] = 2;
//								}
//
//								//North
//								if (labelArray[k][x_new(j - 1, i, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j - 1); i_inProcess.push_back(i);
//									labelArray[k][x_new(j - 1, i, width)] = 2;
//								}
//
//								//Top
//								if (labelArray[k - 1][x_new(j, i, width)] == 3) {
//									k_inProcess.push_back(k - 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//									labelArray[k - 1][x_new(j, i, width)] = 2;
//								}
//
//							}
//							else {
//
//								//West
//								if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//									labelArray[k][x_new(j, i - 1, width)] = 2;
//								}
//
//								//East
//								if (labelArray[k][x_new(j, i + 1, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i + 1);
//									labelArray[k][x_new(j, i + 1, width)] = 2;
//								}
//
//								//North
//								if (labelArray[k][x_new(j - 1, i, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j - 1); i_inProcess.push_back(i);
//									labelArray[k][x_new(j - 1, i, width)] = 2;
//								}
//
//								//South
//								if (labelArray[k][x_new(j + 1, i, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j + 1); i_inProcess.push_back(i);
//									labelArray[k][x_new(j + 1, i, width)] = 2;
//								}
//
//								//Top
//								if (labelArray[k - 1][x_new(j, i, width)] == 3) {
//									k_inProcess.push_back(k - 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//									labelArray[k - 1][x_new(j, i, width)] = 2;
//								}
//
//							}
//						}
//					}
//				}
//			}
//			else {
//				if (i == 0) {
//
//					if (j == 0) {
//
//						//East
//						if (labelArray[k][x_new(j, i + 1, width)] == 3) {
//							k_inProcess.push_back(k + 1); j_inProcess.push_back(j); i_inProcess.push_back(i + 1);
//							labelArray[k][x_new(j, i + 1, width)] = 2;
//						}
//
//						//South
//						if (labelArray[k][x_new(j + 1, i, width)] == 3) {
//							k_inProcess.push_back(k); j_inProcess.push_back(j + 1); i_inProcess.push_back(i);
//							labelArray[k][x_new(j + 1, i, width)] = 2;
//						}
//
//						//Top
//						if (labelArray[k - 1][x_new(j, i, width)] == 3) {
//							k_inProcess.push_back(k - 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//							labelArray[k - 1][x_new(j, i, width)] = 2;
//						}
//
//						//Bottom
//						if (labelArray[k + 1][x_new(j, i, width)] == 3) {
//							k_inProcess.push_back(k + 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//							labelArray[k + 1][x_new(j, i, width)] = 2;
//						}
//
//
//					}
//					else {
//						if (j == width - 1) {
//
//							//West
//							if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//								labelArray[k][x_new(j, i + 1, width)] = 2;
//							}
//
//							//South
//							if (labelArray[k][x_new(j + 1, i, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j + 1); i_inProcess.push_back(i);
//								labelArray[k][x_new(j + 1, i, width)] = 2;
//							}
//
//							//Top
//							if (labelArray[k - 1][x_new(j, i, width)] == 3) {
//								k_inProcess.push_back(k - 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//								labelArray[k - 1][x_new(j, i, width)] = 2;
//							}
//
//							//Bottom
//							if (labelArray[k + 1][x_new(j, i, width)] == 3) {
//								k_inProcess.push_back(k + 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//								labelArray[k + 1][x_new(j, i, width)] = 2;
//							}
//
//						}
//						else {
//
//							//East
//							if (labelArray[k][x_new(j, i + 1, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i + 1);
//								labelArray[k][x_new(j, i + 1, width)] = 2;
//							}
//
//							//West
//							if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//								labelArray[k][x_new(j, i + 1, width)] = 2;
//							}
//
//							//South
//							if (labelArray[k][x_new(j + 1, i, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j + 1); i_inProcess.push_back(i);
//								labelArray[k][x_new(j + 1, i, width)] = 2;
//							}
//
//							//Top
//							if (labelArray[k - 1][x_new(j, i, width)] == 3) {
//								k_inProcess.push_back(k - 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//								labelArray[k - 1][x_new(j, i, width)] = 2;
//							}
//
//							//Bottom
//							if (labelArray[k + 1][x_new(j, i, width)] == 3) {
//								k_inProcess.push_back(k + 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//								labelArray[k + 1][x_new(j, i, width)] = 2;
//							}
//						}
//					}
//				}
//				else {
//					if (i == length - 1) {
//
//						if (j == 0) {
//
//							//West
//							if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//								labelArray[k][x_new(j, i - 1, width)] = 2;
//							}
//
//							//South
//							if (labelArray[k][x_new(j + 1, i, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j + 1); i_inProcess.push_back(i);
//								labelArray[k][x_new(j + 1, i, width)] = 2;
//							}
//
//							//Top
//							if (labelArray[k - 1][x_new(j, i, width)] == 3) {
//								k_inProcess.push_back(k - 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//								labelArray[k - 1][x_new(j, i, width)] = 2;
//							}
//
//							//Bottom
//							if (labelArray[k + 1][x_new(j, i, width)] == 3) {
//								k_inProcess.push_back(k + 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//								labelArray[k + 1][x_new(j, i, width)] = 2;
//							}
//
//						}
//						else {
//							if (j == width - 1) {
//
//								//West
//								if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//									labelArray[k][x_new(j, i - 1, width)] = 2;
//								}
//
//								//North
//								if (labelArray[k][x_new(j - 1, i, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j - 1); i_inProcess.push_back(i);
//									labelArray[k][x_new(j + 1, i, width)] = 2;
//								}
//
//								//Top
//								if (labelArray[k - 1][x_new(j, i, width)] == 3) {
//									k_inProcess.push_back(k - 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//									labelArray[k - 1][x_new(j, i, width)] = 2;
//								}
//
//								//Bottom
//								if (labelArray[k + 1][x_new(j, i, width)] == 3) {
//									k_inProcess.push_back(k + 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//									labelArray[k + 1][x_new(j, i, width)] = 2;
//								}
//
//							}
//							else {
//
//								//West
//								if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//									labelArray[k][x_new(j, i - 1, width)] = 2;
//								}
//
//								//North
//								if (labelArray[k][x_new(j - 1, i, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j - 1); i_inProcess.push_back(i);
//									labelArray[k][x_new(j - 1, i, width)] = 2;
//								}
//
//								//South
//								if (labelArray[k][x_new(j + 1, i, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j + 1); i_inProcess.push_back(i);
//									labelArray[k][x_new(j + 1, i, width)] = 2;
//								}
//
//								//Top
//								if (labelArray[k - 1][x_new(j, i, width)] == 3) {
//									k_inProcess.push_back(k - 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//									labelArray[k - 1][x_new(j, i, width)] = 2;
//								}
//
//								//Bottom
//								if (labelArray[k + 1][x_new(j, i, width)] == 3) {
//									k_inProcess.push_back(k + 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//									labelArray[k + 1][x_new(j, i, width)] = 2;
//								}
//							}
//						}
//					}
//					else {
//						if (j == 0) {
//
//							//West
//							if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//								labelArray[k][x_new(j, i - 1, width)] = 2;
//							}
//
//							//East
//							if (labelArray[k][x_new(j, i + 1, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i + 1);
//								labelArray[k][x_new(j, i + 1, width)] = 2;
//							}
//
//							//South
//							if (labelArray[k][x_new(j + 1, i, width)] == 3) {
//								k_inProcess.push_back(k); j_inProcess.push_back(j + 1); i_inProcess.push_back(i);
//								labelArray[k][x_new(j + 1, i, width)] = 2;
//							}
//
//							//Top
//							if (labelArray[k - 1][x_new(j, i, width)] == 3) {
//								k_inProcess.push_back(k - 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//								labelArray[k - 1][x_new(j, i, width)] = 2;
//							}
//
//							//Bottom
//							if (labelArray[k + 1][x_new(j, i, width)] == 3) {
//								k_inProcess.push_back(k + 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//								labelArray[k + 1][x_new(j, i, width)] = 2;
//							}
//
//						}
//						else {
//							if (j == width - 1) {
//
//								//West
//								if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//									labelArray[k][x_new(j, i - 1, width)] = 2;
//								}
//
//								//East
//								if (labelArray[k][x_new(j, i + 1, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i + 1);
//									labelArray[k][x_new(j, i + 1, width)] = 2;
//								}
//
//								//North
//								if (labelArray[k][x_new(j - 1, i, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j - 1); i_inProcess.push_back(i);
//									labelArray[k][x_new(j - 1, i, width)] = 2;
//								}
//
//								//Top
//								if (labelArray[k - 1][x_new(j, i, width)] == 3) {
//									k_inProcess.push_back(k - 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//									labelArray[k - 1][x_new(j, i, width)] = 2;
//								}
//
//								//Bottom
//								if (labelArray[k + 1][x_new(j, i, width)] == 3) {
//									k_inProcess.push_back(k + 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//									labelArray[k + 1][x_new(j, i, width)] = 2;
//								}
//
//							}
//							else {
//
//								//West
//								if (labelArray[k][x_new(j, i - 1, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i - 1);
//									labelArray[k][x_new(j, i - 1, width)] = 2;
//								}
//
//								//East
//								if (labelArray[k][x_new(j, i + 1, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j); i_inProcess.push_back(i + 1);
//									labelArray[k][x_new(j, i + 1, width)] = 2;
//								}
//
//								//North
//								if (labelArray[k][x_new(j - 1, i, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j - 1); i_inProcess.push_back(i);
//									labelArray[k][x_new(j - 1, i, width)] = 2;
//								}
//
//								//South
//								if (labelArray[k][x_new(j + 1, i, width)] == 3) {
//									k_inProcess.push_back(k); j_inProcess.push_back(j + 1); i_inProcess.push_back(i);
//									labelArray[k][x_new(j + 1, i, width)] = 2;
//								}
//
//								//Top
//								if (labelArray[k - 1][x_new(j, i, width)] == 3) {
//									k_inProcess.push_back(k - 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//									labelArray[k - 1][x_new(j, i, width)] = 2;
//								}
//
//								//Bottom
//								if (labelArray[k + 1][x_new(j, i, width)] == 3) {
//									k_inProcess.push_back(k + 1); j_inProcess.push_back(j); i_inProcess.push_back(i);
//									labelArray[k + 1][x_new(j, i, width)] = 2;
//								}
//
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//
//		//Find the minimal solution
//		minSolution = INFINITY;
//		for (k = 0; k < h_n; k++) {
//			if (minSolution > tempTimeFunc[k]) {
//				minSolution = tempTimeFunc[k];
//				iSol = k; iNew = i_inProcess[k];
//				jSol = k; jNew = j_inProcess[k];
//			}
//		}
//
//		//cout << "\ndistance  : " << minSolution << endl;
//
//		//Set the distance to the processed pixel
//		distanceFuncPtr[x_new(jNew, iNew, width)] = minSolution;
//		labelArray[x_new(jNew, iNew, width)] = 1;
//		i_Processed.push_back(iNew); j_Processed.push_back(jNew);
//
//		//Remove the processed pixel to the stack
//		if (iSol == 0 || jSol == 0) {
//			i_inProcess.erase(i_inProcess.begin());
//			j_inProcess.erase(j_inProcess.begin());
//			tempTimeFunc.erase(tempTimeFunc.begin());
//		}
//		else {
//			i_inProcess.erase(i_inProcess.begin() + iSol);
//			j_inProcess.erase(j_inProcess.begin() + jSol);
//			tempTimeFunc.erase(tempTimeFunc.begin() + jSol);
//		}
//
//		//Compute solution for the neigbors of the selected point
//
//		//STEP 4
//		//Find the neighbors of the processed pixel and compute they time function
//		if (iNew == 0) {
//			if (jNew == 0) {
//
//				jNew_plus = jNew + 1;
//				iNew_plus = iNew + 1;
//
//				//East
//				x = selectX(distanceFuncPtr, height, width, iNew, jNew_plus);
//				y = selectY(distanceFuncPtr, height, width, iNew, jNew_plus);
//				coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_plus, iNew, width)]);
//				dEast = solve2dQuadratic(x, y, coef);
//
//				//South
//				x = selectX(distanceFuncPtr, height, width, iNew_plus, jNew);
//				y = selectY(distanceFuncPtr, height, width, iNew_plus, jNew);
//				coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_plus, width)]);
//				dSouth = solve2dQuadratic(x, y, coef);
//
//			}
//			else {
//				if (jNew == (width - 1)) {
//
//					iNew_plus = iNew + 1;
//					jNew_minus = jNew - 1;
//
//					//West
//					x = selectX(distanceFuncPtr, height, width, iNew, jNew_minus);
//					y = selectY(distanceFuncPtr, height, width, iNew, jNew_minus);
//					coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_minus, iNew, width)]);
//					dWest = solve2dQuadratic(x, y, coef);
//
//					//South
//					x = selectX(distanceFuncPtr, height, width, iNew_plus, jNew);
//					y = selectY(distanceFuncPtr, height, width, iNew_plus, jNew);
//					coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_plus, width)]);
//					dSouth = solve2dQuadratic(x, y, coef);
//
//				}
//				else {
//
//					iNew_plus = iNew + 1;
//					jNew_plus = jNew + 1;
//					jNew_minus = jNew - 1;
//
//					//West
//					x = selectX(distanceFuncPtr, height, width, iNew, jNew_minus);
//					y = selectY(distanceFuncPtr, height, width, iNew, jNew_minus);
//					coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_minus, iNew, width)]);
//					dWest = solve2dQuadratic(x, y, coef);
//
//					//East
//					x = selectX(distanceFuncPtr, height, width, iNew, jNew_plus);
//					y = selectY(distanceFuncPtr, height, width, iNew, jNew_plus);
//					coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_plus, iNew, width)]);
//					dEast = solve2dQuadratic(x, y, coef);
//
//					//South
//					x = selectX(distanceFuncPtr, height, width, iNew_plus, jNew);
//					y = selectY(distanceFuncPtr, height, width, iNew_plus, jNew);
//					coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_plus, width)]);
//					dSouth = solve2dQuadratic(x, y, coef);
//
//				}
//			}
//		}
//		else {
//			if (iNew == (height - 1)) {
//				if (jNew == 0) {
//
//					iNew_minus = iNew - 1;
//					jNew_plus = jNew + 1;
//
//					//North
//					x = selectX(distanceFuncPtr, height, width, iNew_minus, jNew);
//					y = selectY(distanceFuncPtr, height, width, iNew_minus, jNew);
//					coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_minus, width)]);
//					dNorth = solve2dQuadratic(x, y, coef);
//
//					//East
//					x = selectX(distanceFuncPtr, height, width, iNew, jNew_plus);
//					y = selectY(distanceFuncPtr, height, width, iNew, jNew_plus);
//					coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_plus, iNew, width)]);
//					dEast = solve2dQuadratic(x, y, coef);
//
//				}
//				else {
//					if (jNew == (width - 1)) {
//
//						iNew_minus = iNew - 1;
//						jNew_minus = jNew - 1;
//
//						//North
//						x = selectX(distanceFuncPtr, height, width, iNew_minus, jNew);
//						y = selectY(distanceFuncPtr, height, width, iNew_minus, jNew);
//						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_minus, width)]);
//						dNorth = solve2dQuadratic(x, y, coef);
//
//						//West
//						x = selectX(distanceFuncPtr, height, width, iNew, jNew_minus);
//						y = selectY(distanceFuncPtr, height, width, iNew, jNew_minus);
//						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_minus, iNew, width)]);
//						dWest = solve2dQuadratic(x, y, coef);
//
//					}
//					else {
//
//						iNew_minus = iNew - 1;
//						jNew_plus = jNew + 1;
//						jNew_minus = jNew - 1;
//
//						//West
//						x = selectX(distanceFuncPtr, height, width, iNew, jNew_minus);
//						y = selectY(distanceFuncPtr, height, width, iNew, jNew_minus);
//						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_minus, iNew, width)]);
//						dWest = solve2dQuadratic(x, y, coef);
//
//						//North
//						x = selectX(distanceFuncPtr, height, width, iNew_minus, jNew);
//						y = selectY(distanceFuncPtr, height, width, iNew_minus, jNew);
//						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_minus, width)]);
//						dNorth = solve2dQuadratic(x, y, coef);
//
//						//East
//						x = selectX(distanceFuncPtr, height, width, iNew, jNew_plus);
//						y = selectY(distanceFuncPtr, height, width, iNew, jNew_plus);
//						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_plus, iNew, width)]);
//						dEast = solve2dQuadratic(x, y, coef);
//
//					}
//				}
//			}
//			else {
//				if (jNew == 0) {
//
//					iNew_minus = iNew - 1;
//					iNew_plus = iNew + 1;
//					jNew_plus = jNew + 1;
//
//					//North
//					x = selectX(distanceFuncPtr, height, width, iNew_minus, jNew);
//					y = selectY(distanceFuncPtr, height, width, iNew_minus, jNew);
//					coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_minus, width)]);
//					dNorth = solve2dQuadratic(x, y, coef);
//
//					//East
//					x = selectX(distanceFuncPtr, height, width, iNew, jNew_plus);
//					y = selectY(distanceFuncPtr, height, width, iNew, jNew_plus);
//					coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_plus, iNew, width)]);
//					dEast = solve2dQuadratic(x, y, coef);
//
//					//South
//					x = selectX(distanceFuncPtr, height, width, iNew_plus, jNew);
//					y = selectY(distanceFuncPtr, height, width, iNew_plus, jNew);
//					coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_plus, width)]);
//					dSouth = solve2dQuadratic(x, y, coef);
//
//				}
//				else {
//					if (jNew == (width - 1)) {
//
//						iNew_minus = iNew - 1;
//						iNew_plus = iNew + 1;
//						jNew_minus = jNew - 1;
//
//						//North
//						x = selectX(distanceFuncPtr, height, width, iNew_minus, jNew);
//						y = selectY(distanceFuncPtr, height, width, iNew_minus, jNew);
//						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_minus, width)]);
//						dNorth = solve2dQuadratic(x, y, coef);
//
//						//West
//						x = selectX(distanceFuncPtr, height, width, iNew, jNew_minus);
//						y = selectY(distanceFuncPtr, height, width, iNew, jNew_minus);
//						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_minus, iNew, width)]);
//						dWest = solve2dQuadratic(x, y, coef);
//
//						//South
//						x = selectX(distanceFuncPtr, height, width, iNew_plus, jNew);
//						y = selectY(distanceFuncPtr, height, width, iNew_plus, jNew);
//						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_plus, width)]);
//						dSouth = solve2dQuadratic(x, y, coef);
//
//					}
//					else {
//
//						iNew_minus = iNew - 1;
//						iNew_plus = iNew + 1;
//						jNew_plus = jNew + 1;
//						jNew_minus = jNew - 1;
//
//						//North
//						x = selectX(distanceFuncPtr, height, width, iNew_minus, jNew);
//						y = selectY(distanceFuncPtr, height, width, iNew_minus, jNew);
//						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_minus, width)]);
//						dNorth = solve2dQuadratic(x, y, coef);
//
//						//West
//						x = selectX(distanceFuncPtr, height, width, iNew, jNew_minus);
//						y = selectY(distanceFuncPtr, height, width, iNew, jNew_minus);
//						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_minus, iNew, width)]);
//						dWest = solve2dQuadratic(x, y, coef);
//
//						//East
//						x = selectX(distanceFuncPtr, height, width, iNew, jNew_plus);
//						y = selectY(distanceFuncPtr, height, width, iNew, jNew_plus);
//						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew_plus, iNew, width)]);
//						dEast = solve2dQuadratic(x, y, coef);
//
//						//South
//						x = selectX(distanceFuncPtr, height, width, iNew_plus, jNew);
//						y = selectY(distanceFuncPtr, height, width, iNew_plus, jNew);
//						coef = (dataType)(1.0 / potentialFuncPtr[x_new(jNew, iNew_plus, width)]);
//						dSouth = solve2dQuadratic(x, y, coef);
//					}
//				}
//			}
//		}
//
//		//Test and update of neighbors
//		//North
//		if (labelArray[x_new(jNew, iNew_minus, width)] == 3) {
//			distanceFuncPtr[x_new(jNew, iNew_minus, width)] = dNorth;
//			i_inProcess.push_back(iNew_minus); j_inProcess.push_back(jNew);
//			tempTimeFunc.push_back(dNorth);
//			labelArray[x_new(jNew, iNew_minus, width)] = 2;
//		}
//		else {
//			if (labelArray[x_new(jNew, iNew_minus, width)] == 2) {
//				if (dNorth < distanceFuncPtr[x_new(jNew, iNew_minus, width)]) {
//					distanceFuncPtr[x_new(jNew, iNew_minus, width)] = dNorth;
//				}
//			}
//		}
//
//		//West
//		if (labelArray[x_new(jNew_minus, iNew, width)] == 3) {
//			distanceFuncPtr[x_new(jNew_minus, iNew, width)] = dWest;
//			i_inProcess.push_back(iNew); j_inProcess.push_back(jNew_minus);
//			tempTimeFunc.push_back(dWest);
//			labelArray[x_new(jNew_minus, iNew, width)] = 2;
//		}
//		else {
//			if (labelArray[x_new(jNew_minus, iNew, width)] == 2) {
//				if (dWest < distanceFuncPtr[x_new(jNew_minus, iNew, width)]) {
//					distanceFuncPtr[x_new(jNew_minus, iNew, width)] = dWest;
//				}
//			}
//		}
//
//		//East
//		if (labelArray[x_new(jNew_plus, iNew, width)] == 3) {
//			distanceFuncPtr[x_new(jNew_plus, iNew, width)] = dEast;
//			i_inProcess.push_back(iNew); j_inProcess.push_back(jNew_plus);
//			tempTimeFunc.push_back(dEast);
//			labelArray[x_new(jNew_plus, iNew, width)] = 2;
//		}
//		else {
//			if (labelArray[x_new(jNew_plus, iNew, width)] == 2) {
//				if (dEast < distanceFuncPtr[x_new(jNew_plus, iNew, width)]) {
//					distanceFuncPtr[x_new(jNew_plus, iNew, width)] = dEast;
//				}
//			}
//		}
//
//		//South
//		if (labelArray[x_new(jNew, iNew_plus, width)] == 3) {
//			distanceFuncPtr[x_new(jNew, iNew_plus, width)] = dSouth;
//			i_inProcess.push_back(iNew_plus); j_inProcess.push_back(jNew);
//			tempTimeFunc.push_back(dSouth);
//			labelArray[x_new(jNew, iNew_plus, width)] = 2;
//		}
//		else {
//			if (labelArray[x_new(jNew, iNew_plus, width)] == 2) {
//				if (dSouth < distanceFuncPtr[x_new(jNew, iNew_plus, width)]) {
//					distanceFuncPtr[x_new(jNew, iNew_plus, width)] = dSouth;
//				}
//			}
//		}
//
//	}
//	for (k = 0; k < dim2D; k++) {
//		if (labelArray[k] == 3) {
//			cpt++;
//		}
//	}
//	if (cpt == 0) {
//		cout << "\nAll the pixels have been visited" << endl;
//	}
//	else {
//		cout << "\n" << cpt << " pixels haven't been visited" << endl;
//	}
//
//	//cout << "\n" << cpt << " pixels with distance < 0 : " << endl;
//	//cout << "\nNumber of processed Point :" << i_Processed.size() << endl;
//
//	//free(labelArray);
//	delete[] labelArray;
//
//	return true;
//}
