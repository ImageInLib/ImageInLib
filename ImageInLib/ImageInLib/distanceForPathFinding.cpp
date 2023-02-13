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
		j_minus = distanceFuncPtr[x_new(I, J - 1, dimI)];
	}

	if (J == dimJ - 1) {
		j_plus = BIG_VALUE;
	}
	else {
		j_plus = distanceFuncPtr[x_new(I, J + 1, dimI)];
	}

	return min(j_minus, j_plus);
}

dataType selectY(dataType * distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t I, const size_t J) {

	dataType i_minus, i_plus;

	if (I == 0) {
		i_minus = BIG_VALUE;
	}
	else {
		i_minus = distanceFuncPtr[x_new(I - 1, J, dimI)];
	}

	if (I == dimI - 1) {
		i_plus = BIG_VALUE;
	}
	else {
		i_plus = distanceFuncPtr[x_new(I + 1, J, dimI)];
	}

	return min(i_minus, i_plus);
}

bool fastMarching2d(dataType* imageDataPtr, dataType * distanceFuncPtr, const size_t height, const size_t width, Point2D * seedPoints)
{
	//short* labelArray = (short*)malloc(height * width * sizeof(short));
	short * labelArray = new short[height * width];

	if (imageDataPtr == NULL || distanceFuncPtr == NULL || seedPoints == NULL || labelArray == NULL) {
		return false;
	}

	size_t i = 0, j = 0, k = 0, dim2D = height * width, cpt = 0;
	vector<size_t> i_Processed, i_inProcess;
	vector<size_t> j_Processed, j_inProcess;
	dataType x = 0.0, y = 0.0, speed = 1.0, space = 1.0, minSolution = 0.0, coef = 0.0, dist = 0.0;
	size_t nbNeighborsFound = 0;
	size_t h_n = 0, iSol = 0, jSol = 0, iNew = 0, jNew = 0, iNew_minus = 0, iNew_plus = 0, jNew_minus = 0, jNew_plus = 0;
	vector<dataType> tempTimeFunc;
	dataType dNorth = 0.0, dSouth = 0.0, dEast = 0.0, dWest = 0.0;

	//cout << "\nNumber of pixels : " << height * width  << endl;

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
	distanceFuncPtr[x_new(i, j, height)] = 0;
	i_Processed.push_back(i); j_Processed.push_back(j);

	size_t iminus = i - 1, iplus = i + 1, jminus = j - 1, jplus = j + 1;
	if (i == 0) {
		if (j == 0) {

			//East
			if (labelArray[x_new(i, jplus, height)] == 3) {
				i_inProcess.push_back(i); j_inProcess.push_back(jplus);
				labelArray[x_new(i, jplus, height)] = 2;
				nbNeighborsFound++;
			}
			//South
			if (labelArray[x_new(iplus, j, height)] == 3) {
				i_inProcess.push_back(iplus); j_inProcess.push_back(j);
				labelArray[x_new(iplus, j, height)] = 2;
				nbNeighborsFound++;
			}
		}
		else {
			if ( j == (width - 1) ) {
				//West
				if (labelArray[x_new(i, jminus, height)] == 3) {
					i_inProcess.push_back(i); j_inProcess.push_back(jminus);
					labelArray[x_new(i, jminus, height)] = 2;
					nbNeighborsFound++;
				}
				//South
				if (labelArray[x_new(iplus, j, height)] == 3) {
					i_inProcess.push_back(iplus); j_inProcess.push_back(j);
					labelArray[x_new(iplus, j, height)] = 2;
					nbNeighborsFound++;
				}
			}
			else {

				//East
				if (labelArray[x_new(i, jplus, height)] == 3) {
					i_inProcess.push_back(i); j_inProcess.push_back(jplus);
					labelArray[x_new(i, jplus, height)] = 2;
					nbNeighborsFound++;
				}
				//West
				if (labelArray[x_new(i, jminus, height)] == 3) {
					i_inProcess.push_back(i); j_inProcess.push_back(jminus);
					labelArray[x_new(i, jminus, height)] = 2;
					nbNeighborsFound++;
				}
				//South
				if (labelArray[x_new(iplus, j, height)] == 3) {
					i_inProcess.push_back(iplus); j_inProcess.push_back(j);
					labelArray[x_new(iplus, j, height)] = 2;
					nbNeighborsFound++;
				}

			}
		}
	}
	else {
		if (i == (height - 1)) {
			if (j == 0) {

				//East
				if (labelArray[x_new(i, jplus, height)] == 3) {
					i_inProcess.push_back(i); j_inProcess.push_back(jplus);
					labelArray[x_new(i, jplus, height)] = 2;
					nbNeighborsFound++;
				}
				//North
				if (labelArray[x_new(iminus, j, height)] == 3) {
					i_inProcess.push_back(iminus); j_inProcess.push_back(j);
					labelArray[x_new(iminus, j, height)] = 2;
					nbNeighborsFound++;
				}
				
			}
			else {
				if (j == (width - 1)) {

					//North
					if (labelArray[x_new(iminus, j, height)] == 3) {
						i_inProcess.push_back(iminus); j_inProcess.push_back(j);
						labelArray[x_new(iminus, j, height)] = 2;
						nbNeighborsFound++;
					}
					//West
					if (labelArray[x_new(i, jminus, height)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(jminus);
						labelArray[x_new(i, jminus, height)] = 2;
						nbNeighborsFound++;
					}

				}
				else {

					//East
					if (labelArray[x_new(i, jplus, height)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(jplus);
						labelArray[x_new(i, jplus, height)] = 2;
						nbNeighborsFound++;
					}
					//North
					if (labelArray[x_new(iminus, j, height)] == 3) {
						i_inProcess.push_back(iminus); j_inProcess.push_back(j);
						labelArray[x_new(iminus, j, height)] = 2;
						nbNeighborsFound++;
					}
					//West
					if (labelArray[x_new(i, jminus, height)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(jminus);
						labelArray[x_new(i, jminus, height)] = 2;
						nbNeighborsFound++;
					}

				}
			}
		}
		else {
			if (j == 0) {

				//East
				if (labelArray[x_new(i, jplus, height)] == 3) {
					i_inProcess.push_back(i); j_inProcess.push_back(jplus);
					labelArray[x_new(i, jplus, height)] = 2;
					nbNeighborsFound++;
				}
				//North
				if (labelArray[x_new(iminus, j, height)] == 3) {
					i_inProcess.push_back(iminus); j_inProcess.push_back(j);
					labelArray[x_new(iminus, j, height)] = 2;
					nbNeighborsFound++;
				}
				//South
				if (labelArray[x_new(iplus, j, height)] == 3) {
					i_inProcess.push_back(iplus); j_inProcess.push_back(j);
					labelArray[x_new(iplus, j, height)] = 2;
					nbNeighborsFound++;
				}

			}
			else {
				if (j == (width - 1) ) {

					//North
					if (labelArray[x_new(iminus, j, height)] == 3) {
						i_inProcess.push_back(iminus); j_inProcess.push_back(j);
						labelArray[x_new(iminus, j, height)] = 2;
						nbNeighborsFound++;
					}
					//West
					if (labelArray[x_new(i, jminus, height)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(jminus);
						labelArray[x_new(i, jminus, height)] = 2;
						nbNeighborsFound++;
					}
					//South
					if (labelArray[x_new(iplus, j, height)] == 3) {
						i_inProcess.push_back(iplus); j_inProcess.push_back(j);
						labelArray[x_new(iplus, j, height)] = 2;
						nbNeighborsFound++;
					}
				}
				else {

					//East
					if (labelArray[x_new(i, jplus, height)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(jplus);
						labelArray[x_new(i, jplus, height)] = 2;
						nbNeighborsFound++;
					}
					//North
					if (labelArray[x_new(iminus, j, height)] == 3) {
						i_inProcess.push_back(iminus); j_inProcess.push_back(j);
						labelArray[x_new(iminus, j, height)] = 2;
						nbNeighborsFound++;
					}
					//West
					if (labelArray[x_new(i, jminus, height)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(jminus);
						labelArray[x_new(i, jminus, height)] = 2;
						nbNeighborsFound++;
					}
					//South
					if (labelArray[x_new(iplus, j, height)] == 3) {
						i_inProcess.push_back(iplus); j_inProcess.push_back(j);
						labelArray[x_new(iplus, j, height)] = 2;
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

		coef = pow((space / speed), 2);

		dist = solve2dQuadratic(x, y, coef);
		tempTimeFunc.push_back(dist);
	}

	//Update the label of seed point as processed
	labelArray[x_new(i, j, height)] = 1;

	//cout << "Number of Neigbors found : " << nbNeighborsFound << endl;
	//---------------------End of STEP 2 -------------------------------------

	//STEP 3
	while (i_inProcess.size() != 0 || j_inProcess.size() != 0 ) {

		h_n = i_inProcess.size();
		// i_inProcess and j_inProcess have the same size

		//Find the minimal solution
		minSolution = INFINITY;
		for (k = 0; k < h_n; k++) {
			if (minSolution >= tempTimeFunc[k]) {
				minSolution = tempTimeFunc[k];
				iSol = k; iNew = i_inProcess[k];
				jSol = k; jNew = j_inProcess[k];
			}
		}

		//cout << "\ndistance  : " << minSolution << endl;
		
		//Set the distance to the processed pixel
		distanceFuncPtr[x_new(iNew, jNew, height)] = minSolution;
		labelArray[x_new(iNew, jNew, height)] = 1;
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
		 
		//iNew_minus = iNew - 1; 
		//iNew_plus = iNew + 1;
		//jNew_minus = jNew - 1; 
		//jNew_plus = jNew + 1;

		//STEP 4
		//Find the neighbors of the processed pixel and compute they time function
		if (iNew == 0) {
			if (jNew == 0) {

				jNew_plus = jNew + 1; 
				iNew_plus = iNew + 1;

				//East
				x = selectX(distanceFuncPtr, height, width, iNew, jNew_plus);
				y = selectY(distanceFuncPtr, height, width, iNew, jNew_plus);
				coef = pow((space / speed), 2);
				dEast = solve2dQuadratic(x, y, coef);

				//South
				x = selectX(distanceFuncPtr, height, width, iNew_plus, jNew);
				y = selectY(distanceFuncPtr, height, width, iNew_plus, jNew);
				coef = pow((space / speed), 2);
				dSouth = solve2dQuadratic(x, y, coef);

			}
			else {
				if (jNew == (width - 1) ) {

					iNew_plus = iNew + 1;
					jNew_minus = jNew - 1;

					//West
					x = selectX(distanceFuncPtr, height, width, iNew, jNew_minus);
					y = selectY(distanceFuncPtr, height, width, iNew, jNew_minus);
					coef = pow((space / speed), 2);
					dWest = solve2dQuadratic(x, y, coef);

					//South
					x = selectX(distanceFuncPtr, height, width, iNew_plus, jNew);
					y = selectY(distanceFuncPtr, height, width, iNew_plus, jNew);
					coef = pow((space / speed), 2);
					dSouth = solve2dQuadratic(x, y, coef);

				}
				else {

					iNew_plus = iNew + 1;
					jNew_plus = jNew + 1;
					jNew_minus = jNew - 1;

					//West
					x = selectX(distanceFuncPtr, height, width, iNew, jNew_minus);
					y = selectY(distanceFuncPtr, height, width, iNew, jNew_minus);
					coef = pow((space / speed), 2);
					dWest = solve2dQuadratic(x, y, coef);

					//East
					x = selectX(distanceFuncPtr, height, width, iNew, jNew_plus);
					y = selectY(distanceFuncPtr, height, width, iNew, jNew_plus);
					coef = pow((space / speed), 2);
					dEast = solve2dQuadratic(x, y, coef);

					//South
					x = selectX(distanceFuncPtr, height, width, iNew_plus, jNew);
					y = selectY(distanceFuncPtr, height, width, iNew_plus, jNew);
					coef = pow((space / speed), 2);
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
					coef = pow((space / speed), 2);
					dNorth = solve2dQuadratic(x, y, coef);

					//East
					x = selectX(distanceFuncPtr, height, width, iNew, jNew_plus);
					y = selectY(distanceFuncPtr, height, width, iNew, jNew_plus);
					coef = pow((space / speed), 2);
					dEast = solve2dQuadratic(x, y, coef);

				}
				else {
					if (jNew == (width - 1) ) {

						iNew_minus = iNew - 1;
						jNew_minus = jNew - 1;

						//North
						x = selectX(distanceFuncPtr, height, width, iNew_minus, jNew);
						y = selectY(distanceFuncPtr, height, width, iNew_minus, jNew);
						coef = pow((space / speed), 2);
						dNorth = solve2dQuadratic(x, y, coef);

						//West
						x = selectX(distanceFuncPtr, height, width, iNew, jNew_minus);
						y = selectY(distanceFuncPtr, height, width, iNew, jNew_minus);
						coef = pow((space / speed), 2);
						dWest = solve2dQuadratic(x, y, coef);

					}
					else {

						iNew_minus = iNew - 1;
						jNew_plus = jNew + 1;
						jNew_minus = jNew - 1;

						//West
						x = selectX(distanceFuncPtr, height, width, iNew, jNew_minus);
						y = selectY(distanceFuncPtr, height, width, iNew, jNew_minus);
						coef = pow((space / speed), 2);
						dWest = solve2dQuadratic(x, y, coef);

						//North
						x = selectX(distanceFuncPtr, height, width, iNew_minus, jNew);
						y = selectY(distanceFuncPtr, height, width, iNew_minus, jNew);
						coef = pow((space / speed), 2);
						dNorth = solve2dQuadratic(x, y, coef);

						//East
						x = selectX(distanceFuncPtr, height, width, iNew, jNew_plus);
						y = selectY(distanceFuncPtr, height, width, iNew, jNew_plus);
						coef = pow((space / speed), 2);
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
					coef = pow((space / speed), 2);
					dNorth = solve2dQuadratic(x, y, coef);

					//East
					x = selectX(distanceFuncPtr, height, width, iNew, jNew_plus);
					y = selectY(distanceFuncPtr, height, width, iNew, jNew_plus);
					coef = pow((space / speed), 2);
					dEast = solve2dQuadratic(x, y, coef);

					//South
					x = selectX(distanceFuncPtr, height, width, iNew_plus, jNew);
					y = selectY(distanceFuncPtr, height, width, iNew_plus, jNew);
					coef = pow((space / speed), 2);
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
						coef = pow((space / speed), 2);
						dNorth = solve2dQuadratic(x, y, coef);

						//West
						x = selectX(distanceFuncPtr, height, width, iNew, jNew_minus);
						y = selectY(distanceFuncPtr, height, width, iNew, jNew_minus);
						coef = pow((space / speed), 2);
						dWest = solve2dQuadratic(x, y, coef);

						//South
						x = selectX(distanceFuncPtr, height, width, iNew_plus, jNew);
						y = selectY(distanceFuncPtr, height, width, iNew_plus, jNew);
						coef = pow((space / speed), 2);
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
						coef = pow((space / speed), 2);
						dNorth = solve2dQuadratic(x, y, coef);

						//West
						x = selectX(distanceFuncPtr, height, width, iNew, jNew_minus);
						y = selectY(distanceFuncPtr, height, width, iNew, jNew_minus);
						coef = pow((space / speed), 2);
						dWest = solve2dQuadratic(x, y, coef);

						//East
						x = selectX(distanceFuncPtr, height, width, iNew, jNew_plus);
						y = selectY(distanceFuncPtr, height, width, iNew, jNew_plus);
						coef = pow((space / speed), 2);
						dEast = solve2dQuadratic(x, y, coef);

						//South
						x = selectX(distanceFuncPtr, height, width, iNew_plus, jNew);
						y = selectY(distanceFuncPtr, height, width, iNew_plus, jNew);
						coef = pow((space / speed), 2);
						dSouth = solve2dQuadratic(x, y, coef);
					}
				}
			}
		}

		//Test and update of neighbors
		//North
		if (labelArray[x_new(iNew_minus, jNew, height)] == 3) {
			distanceFuncPtr[x_new(iNew_minus, jNew, height)] = dNorth;
			i_inProcess.push_back(iNew_minus); j_inProcess.push_back(jNew);
			tempTimeFunc.push_back(dNorth);
			labelArray[x_new(iNew_minus, jNew, height)] = 2;
		}
		else {
			if (labelArray[x_new(iNew_minus, jNew, height)] == 2) {
				if (dNorth < distanceFuncPtr[x_new(iNew_minus, jNew, height)]) {
					distanceFuncPtr[x_new(iNew_minus, jNew, height)] = dNorth;
				}
			}
		}

		//West
		if (labelArray[x_new(iNew, jNew_minus, height)] == 3) {
			distanceFuncPtr[x_new(iNew, jNew_minus, height)] = dWest;
			i_inProcess.push_back(iNew); j_inProcess.push_back(jNew_minus);
			tempTimeFunc.push_back(dWest);
			labelArray[x_new(iNew, jNew_minus, height)] = 2;
		}
		else {
			if (labelArray[x_new(iNew, jNew_minus, height)] == 2) {
				if (dWest < distanceFuncPtr[x_new(iNew, jNew_minus, height)]) {
					distanceFuncPtr[x_new(iNew, jNew_minus, height)] = dWest;
				}
			}
		}

		//East
		if (labelArray[x_new(iNew, jNew_plus, height)] == 3) {
			distanceFuncPtr[x_new(iNew, jNew_plus, height)] = dEast;
			i_inProcess.push_back(iNew); j_inProcess.push_back(jNew_plus);
			tempTimeFunc.push_back(dEast);
			labelArray[x_new(iNew, jNew_plus, height)] = 2;
		}
		else {
			if (labelArray[x_new(iNew, jNew_plus, height)] == 2) {
				if (dEast < distanceFuncPtr[x_new(iNew, jNew_plus, height)]) {
					distanceFuncPtr[x_new(iNew, jNew_plus, height)] = dEast;
				}
			}
		}

		//South
		if (labelArray[x_new(iNew_plus, jNew, height)] == 3) {
			distanceFuncPtr[x_new(iNew_plus, jNew, height)] = dSouth;
			i_inProcess.push_back(iNew_plus); j_inProcess.push_back(jNew);
			tempTimeFunc.push_back(dSouth);
			labelArray[x_new(iNew_plus, jNew, height)] = 2;
		}
		else {
			if (labelArray[x_new(iNew_plus, jNew, height)] == 2) {
				if (dSouth < distanceFuncPtr[x_new(iNew_plus, jNew, height)]) {
					distanceFuncPtr[x_new(iNew_plus, jNew, height)] = dSouth;
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
