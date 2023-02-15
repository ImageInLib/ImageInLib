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

// aU^2 -2U(X+Y) + (X^2 + Y^2 - W) = 0
dataType solve2dQuadratic(dataType X, dataType Y, dataType W) {

	dataType sol, a, b, c, delta;

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

dataType selectX(dataType * distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t I, const size_t J) {

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

dataType selectY(dataType * distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t I, const size_t J) {

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

bool computeImageGradient(dataType * imageDataPtr, dataType * imageGradient, const size_t height, const size_t width, dataType h) {
	
	if (imageDataPtr == NULL || imageGradient == NULL) {
		return false;
	}

	size_t i, j, i_ext, j_ext, xd, xd_ext, dim2d = height * width;
	dataType ux = 0.0, uy = 0.0, diff = 0.0, racine = 0.0;

	////Manual reflexion
	//size_t height_ext = height + 2, width_ext = width + 2, dim2d_ext = height_ext * width_ext;
	//dataType * extendedImage = new dataType[dim2d_ext];
	////Set corners
	//extendedImage[0] = imageDataPtr[0];
	//extendedImage[x_new(0, height_ext - 1, width_ext)] = imageDataPtr[x_new(0, height - 1, width)];
	//extendedImage[x_new(width_ext - 1, 0, width_ext)] = imageDataPtr[x_new(width - 1, 0, width)];
	//extendedImage[x_new(width_ext - 1, height_ext - 1, width_ext)] = imageDataPtr[x_new(width - 1, height - 1, width)];
	////Set the middle
	//for (i = 0, i_ext = 1; i < height; i++, i_ext++) {
	//	for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
	//		xd = x_new(j, i, width);
	//		xd_ext = x_new(j_ext, i_ext, width_ext);
	//		extendedImage[xd_ext] = imageDataPtr[xd];
	//	}
	//}
	////Set top parts of the extended Image
	//for (i = 0; i < height; i++) {
	//	extendedImage[x_new(0, i + 1, width_ext)] = imageDataPtr[x_new(0, i, width)];
	//	extendedImage[x_new(width_ext - 1, i + 1, width_ext)] = imageDataPtr[x_new(width - 1, i, width)];
	//}
	//for (j = 0; j < width; j++) {
	//	extendedImage[x_new(j + 1, 0, width_ext)] = imageDataPtr[x_new(j, 0, width)];
	//	extendedImage[x_new(j + 1, height_ext - 1, width_ext)] = imageDataPtr[x_new(j, height - 1, width)];
	//}
	////Using central difference to compute the gradient
	//for (i = 0, i_ext = 1; i < height; i++, i_ext++) {
	//	for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
	//		xd = x_new(j, i, width);
	//		ux = (extendedImage[x_new(j + 1, i, width_ext)] - extendedImage[x_new(j - 1, i, width_ext)]) / (2 * h);
	//		uy = (extendedImage[x_new(j, i + 1, width_ext)] - extendedImage[x_new(j, i - 1, width_ext)]) / (2 * h);
	//		diff = (dataType)(ux * ux + uy * uy);
	//		racine = sqrt(diff);
	//		//imageGradient[xd] = (dataType)sqrt(pow(ux, 2) + pow(uy, 2));
	//		imageGradient[xd] = (dataType)sqrt(diff);
	//	}
	//}

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
			imageGradient[xd] = sqrt(ux * ux + uy * uy);
		}
	}

	//delete[] extendedImage;

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

bool computePotential(dataType * imageDataPtr, dataType * potentialFuncPtr, const size_t height, const size_t width, Point2D * seedPoints) {

	if (imageDataPtr == NULL || potentialFuncPtr == NULL || seedPoints == NULL)
		return false;

	//dataType* imageGradient = new dataType[height * width];
	//computeImageGradient(imageDataPtr, imageGradient, height, width, 1.0);
	//dataType norm_gradient = computeImageNorm2d(imageDataPtr, height, width);

	size_t i, j, in = seedPoints->y, jn = seedPoints->x;

	dataType seedVal = imageDataPtr[x_new(jn, in, width)];
	size_t currentIndx = 0;

	//Computation of potential function
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			currentIndx = x_new(j, i, width);
			potentialFuncPtr[currentIndx] = (dataType)(0.01 + abs(seedVal - imageDataPtr[currentIndx]));
			//potentialFuncPtr[currentIndx] = (dataType)(1 / (1 + 1 * norm_gradient * norm_gradient) );
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
	//		potentialFuncPtr[currentIndx] = (dataType)(0.01 + potentialFuncPtr[currentIndx] / max_potential);
	//	}
	//}

	//delete[] imageGradient;

	return true;
}

bool fastMarching2d(dataType * imageDataPtr, dataType * distanceFuncPtr, dataType * potentialFuncPtr, const size_t height, const size_t width, Point2D * seedPoints)
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
	i = seedPoints->y; j = seedPoints->x;
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

	dataType coefSpeed = 1.0;

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

	//for (k = 0; k < dim2D; k++) {
	//	if (distanceFuncPtr[k] < 0) {
	//		distanceFuncPtr[k] = 0;
	//		cpt++;
	//	}
	//}
	//cout << "\n" << cpt << " pixels with distance < 0 : " << endl;
	//cout << "\nNumber of processed Point :" << i_Processed.size() << endl;

	//free(labelArray);
	delete[] labelArray;

	return true;
}

bool shortestPath2d(dataType * distanceFuncPtr, dataType * resultedPath, const size_t height, const size_t width, dataType h, Point2D* seedPoints) {

	if (distanceFuncPtr == NULL || resultedPath == NULL || seedPoints == NULL)
		return false;

	size_t i, j, xd, dim2d = height * width;
	dataType norme_of_gradient = 0.0, sum_gradient = 0.0, tau = 0.8, dist_min = INFINITY, tol = 0.01;
	size_t i_init = seedPoints[0].y, j_init = seedPoints[0].x, i_end = seedPoints[1].y, j_end = seedPoints[1].x;
	size_t i_current, j_current, iSol, jSol;

	dataType * imageGradient = new dataType[dim2d];

	computeImageGradient(distanceFuncPtr, imageGradient, height, width, h);

	//Compute norm of gradient
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			xd = x_new(j, i, width);
			sum_gradient = sum_gradient + imageGradient[xd] * imageGradient[xd];
		}
	}
	norme_of_gradient = sqrt(sum_gradient);

	//Compute the normalization
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			xd = x_new(j, i, width);
			if (norme_of_gradient != 0) {
				imageGradient[xd] = (dataType(imageGradient[xd] / norme_of_gradient));
			}
		}
	}

	//Find the closest point till the last point
	size_t cpt = 1;
	i_current = i_end; j_current = j_end;
	resultedPath[x_new(j_current, i_current, width)] = 1;
	do{
		iSol = (size_t)(i_current - tau * imageGradient[x_new(j_current, i_current, width)]);
		jSol = (size_t)(j_current - tau * imageGradient[x_new(j_current, i_current, width)]);
		dist_min = sqrt((iSol - i_init) * (iSol - i_init) + (jSol - j_init) * (jSol - j_init));
		resultedPath[x_new(jSol, iSol, width)] = 1;
		i_current = iSol; j_current = jSol;
		cpt++;
	}
	while(dist_min > tol && cpt < 20000);

	delete[] imageGradient;

	return true;
}

bool findPath2d(dataType* distanceFuncPtr, dataType* resultedPath, const size_t height, const size_t width, Point2D* seedPoints) {

	if (distanceFuncPtr == NULL || resultedPath == NULL || seedPoints == NULL)
		return false;

	size_t k, hn, coordWest, coordEast, coordNorth, coordSouth, coordCurrent, cpt = 0;
	size_t i_init = seedPoints[0].y, j_init = seedPoints[0].x, i_end = seedPoints[1].y, j_end = seedPoints[1].x;
	dataType min_value, pathLenght = sqrt((i_end - i_init) * (i_end - i_init) + (j_end - j_init) * (j_end - j_init)), tol = 0.01;
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
