#include <iostream>
#include <climits>
#include <crtdbg.h>
#include <corecrt_malloc.h>
#include<cmath>
#include <omp.h>
#include<vector>
#include "distanceForPathFinding.h"

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

	size_t i, j, i_ext, j_ext, xd, dim2d = height * width, height_ext = height + 2, width_ext = width + 2, dim2d_ext =  height_ext * width_ext;
	dataType ux = 0.0, uy = 0.0;

	dataType * extendedArray = new dataType[dim2d_ext];
	if (extendedArray == NULL)
		return false;

	//Initialization
	for (i = 0; i < dim2d_ext; i++) {
		extendedArray[i] = 0;
	}

	copyDataTo2dExtendedArea(imageDataPtr, extendedArray, height, width);
	reflection2D(extendedArray, height_ext, width_ext);
	 
	//Using central difference to compute the gradient
	for (i = 0, i_ext = 1; i < height; i++, i_ext++) {
		for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
			xd = x_new(i, j, height);
			ux = (extendedArray[x_new(i_ext, j_ext + 1, height_ext)] - extendedArray[x_new(i_ext, j_ext - 1, height_ext)]) / (2 * h);
			uy = (extendedArray[x_new(i_ext + 1, j_ext, height_ext)] - extendedArray[x_new(i_ext - 1, j_ext, height_ext)]) / (2 * h);
			imageGradient[xd] = sqrt(pow(ux, 2) + pow(uy, 2));
		}
	}

	delete[] extendedArray;

	return true;
}

bool computePotential(dataType * imageDataPtr, dataType * potentialFuncPtr, const size_t height, const size_t width, Point2D * seedPoints) {

	if (imageDataPtr == NULL || potentialFuncPtr == NULL || seedPoints == NULL)
		return false;

	size_t i, j, in = seedPoints->y, jn = seedPoints->x;

	dataType seedVal = imageDataPtr[x_new(jn, in, width)];
	size_t currentIndx = 0;

	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			currentIndx = x_new(j, i, width);
			potentialFuncPtr[currentIndx] = (dataType)(20.0 + abs(seedVal - imageDataPtr[currentIndx]));
			//potentialFuncPtr[currentIndx] = 0.01  + imageDataPtr[currentIndx];
			/*if(i > 200)
				potentialFuncPtr[currentIndx] = (dataType)(0.001 + abs(seedVal - imageDataPtr[currentIndx]));
			else
				potentialFuncPtr[currentIndx] = 0.1;*/
		}
	}

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
