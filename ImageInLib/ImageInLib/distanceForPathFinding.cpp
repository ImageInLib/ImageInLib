#include <iostream>
#include <climits>
#include <crtdbg.h>
#include <corecrt_malloc.h>
#include<cmath>
#include <omp.h>
#include<vector>
#include "distanceForPathFinding.h"

using namespace std;

bool fastMarching2d(dataType* imageDataPtr, dataType* distanceFuncPtr, const size_t height, const size_t width, Point2D* seedPoints)
{

	short* labelArray = (short*)malloc(height * width * sizeof(short));

	if (imageDataPtr == NULL || distanceFuncPtr == NULL || seedPoints == NULL || labelArray == NULL) {
		return false;
	}

	size_t i = 0, j = 0, k = 0, dim2D = width * height, cpt = 0;
	vector<size_t> i_Processed, i_inProcess;
	vector<size_t> j_Processed, j_inProcess;
	dataType x, y, speed = 1.0, space = 1.0, minSolution = 0.0, coef = 0.0, dist = 0.0;
	size_t nbNeighborsFound = 0;
	size_t h_n = 0, w_n = 0, iSol = 0, jSol = 0, iNew = 0, jNew = 0;
	vector<dataType> tempTimeFunc;
	dataType a, b, c, delta;

	//STEP 1
	//Initialization : distance is set to infinity
	//In labelAray we have : 1 ---> already processed, 2 ---> in process and 3 ---> not processed
	for (k = 0; k < dim2D; k++) {
		distanceFuncPtr[k] = INFINITY;
		labelArray[k] = 3;
	}
	//--------------------End of STEP 1 -----------------------------------

	//STEP 2
	//Add the neighbors of the seed point in the list of pixels to be processed
	i = seedPoints->x; j = seedPoints->y;
	distanceFuncPtr[x_new(i, j, height)] = 0;
	i_Processed.push_back(i); j_Processed.push_back(j);

	if (i == 0) {
		if (j == 0) {
			//East
			if (labelArray[x_new(i, j + 1, height)] == 3) {
				i_inProcess.push_back(i); j_inProcess.push_back(j + 1);
				labelArray[x_new(i, j + 1, height)] = 2;
				nbNeighborsFound++;
			}
			//South
			if (labelArray[x_new(i + 1, j, height)] == 3) {
				i_inProcess.push_back(i + 1); j_inProcess.push_back(j);
				labelArray[x_new(i + 1, j, height)] = 2;
				nbNeighborsFound++;
			}
		}
		else {
			if (j == width - 1) {
				//West
				if (labelArray[x_new(i, j - 1, height)] == 3) {
					i_inProcess.push_back(i); j_inProcess.push_back(j - 1);
					labelArray[x_new(i, j - 1, height)] = 2;
					nbNeighborsFound++;
				}
				//South
				if (labelArray[x_new(i + 1, j, height)] == 3) {
					i_inProcess.push_back(i + 1); j_inProcess.push_back(j);
					labelArray[x_new(i + 1, j, height)] = 2;
					nbNeighborsFound++;
				}
			}
			else {
				//West
				if (labelArray[x_new(i, j - 1, height)] == 3) {
					i_inProcess.push_back(i); j_inProcess.push_back(j - 1);
					labelArray[x_new(i, j - 1, height)] = 2;
					nbNeighborsFound++;
				}
				//East
				if (labelArray[x_new(i, j + 1, height)] == 3) {
					i_inProcess.push_back(i); j_inProcess.push_back(j + 1);
					labelArray[x_new(i, j + 1, height)] = 2;
					nbNeighborsFound++;
				}
				//South
				if (labelArray[x_new(i + 1, j, height)] == 3) {
					i_inProcess.push_back(i + 1); j_inProcess.push_back(j);
					labelArray[x_new(i + 1, j, height)] = 2;
					nbNeighborsFound++;
				}
			}
		}
	}
	else {
		if (i == height - 1) {
			if (j == 0) {
				//North
				if (labelArray[x_new(i - 1, j, height)] == 3) {
					i_inProcess.push_back(i - 1); j_inProcess.push_back(j);
					labelArray[x_new(i - 1, j, height)] = 2;
					nbNeighborsFound++;
				}
				//East
				if (labelArray[x_new(i, j + 1, height)] == 3) {
					i_inProcess.push_back(i); j_inProcess.push_back(j + 1);
					labelArray[x_new(i, j + 1, height)] = 2;
					nbNeighborsFound++;
				}
			}
			else {
				if (j == width - 1) {
					//North
					if (labelArray[x_new(i - 1, j, height)] == 3) {
						i_inProcess.push_back(i - 1); j_inProcess.push_back(j);
						labelArray[x_new(i - 1, j, height)] = 2;
						nbNeighborsFound++;
					}
					//West
					if (labelArray[x_new(i, j - 1, height)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(j - 1);
						labelArray[x_new(i, j - 1, height)] = 2;
						nbNeighborsFound++;
					}
				}
				else {
					//West
					if (labelArray[x_new(i, j - 1, height)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(j - 1);
						labelArray[x_new(i, j - 1, height)] = 2;
						nbNeighborsFound++;
					}
					//North
					if (labelArray[x_new(i - 1, j, height)] == 3) {
						i_inProcess.push_back(i - 1); j_inProcess.push_back(j);
						labelArray[x_new(i - 1, j, height)] = 2;
						nbNeighborsFound++;
					}
					//East
					if (labelArray[x_new(i, j + 1, height)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(j + 1);
						labelArray[x_new(i, j + 1, height)] = 2;
						nbNeighborsFound++;
					}
				}
			}
		}
		else {
			if (j == 0) {
				//North
				if (labelArray[x_new(i - 1, j, height)] == 3) {
					i_inProcess.push_back(i - 1); j_inProcess.push_back(j);
					labelArray[x_new(i - 1, j, height)] = 2;
					nbNeighborsFound++;
				}
				//East
				if (labelArray[x_new(i, j + 1, height)] == 3) {
					i_inProcess.push_back(i); j_inProcess.push_back(j + 1);
					labelArray[x_new(i, j + 1, height)] = 2;
					nbNeighborsFound++;
				}
				//South
				if (labelArray[x_new(i + 1, j, height)] == 3) {
					i_inProcess.push_back(i + 1); j_inProcess.push_back(j);
					labelArray[x_new(i + 1, j, height)] = 2;
					nbNeighborsFound++;
				}
			}
			else {
				if (j == width - 1) {
					//North
					if (labelArray[x_new(i - 1, j, height)] == 3) {
						i_inProcess.push_back(i - 1); j_inProcess.push_back(j);
						labelArray[x_new(i - 1, j, height)] = 2;
						nbNeighborsFound++;
					}
					//West
					if (labelArray[x_new(i, j - 1, height)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(j - 1);
						labelArray[x_new(i, j - 1, height)] = 2;
						nbNeighborsFound++;
					}
					//South
					if (labelArray[x_new(i + 1, j, height)] == 3) {
						i_inProcess.push_back(i + 1); j_inProcess.push_back(j);
						labelArray[x_new(i + 1, j, height)] = 2;
						nbNeighborsFound++;
					}
				}
				else {
					//North
					if (labelArray[x_new(i - 1, j, height)] == 3) {
						i_inProcess.push_back(i - 1); j_inProcess.push_back(j);
						labelArray[x_new(i - 1, j, height)] = 2;
						nbNeighborsFound++;
					}
					//West
					if (labelArray[x_new(i, j - 1, height)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(j - 1);
						labelArray[x_new(i, j - 1, height)] = 2;
						nbNeighborsFound++;
					}
					//East
					if (labelArray[x_new(i, j + 1, height)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(j + 1);
						labelArray[x_new(i, j + 1, height)] = 2;
						nbNeighborsFound++;
					}
					//South
					if (labelArray[x_new(i + 1, j, height)] == 3) {
						i_inProcess.push_back(i + 1); j_inProcess.push_back(j);
						labelArray[x_new(i + 1, j, height)] = 2;
						nbNeighborsFound++;
					}
				}
			}
		}
	}

	//Update the label of seed point as processed
	labelArray[x_new(i, j, height)] = 1;

	cout << "Number of Neigbors found : " << nbNeighborsFound << endl;
	//---------------------End of STEP 2 -------------------------------------

	//STEP 3
	while (i_inProcess.size() != 0 || j_inProcess.size() != 0 ) {

		h_n = i_inProcess.size();
		//w_n = j_inProcess.size(); // they have the same value

		//Solve the Eikonale equation for all the neighbors in the stack
		for (k = 0; k < h_n; k++) {

			if (i_inProcess[k] == 0 || i_inProcess[k] == (height - 1) ) {
				x = 0;
			}
			else {
				x = min(distanceFuncPtr[x_new(i_inProcess[k] + 1, j_inProcess[k], height)], distanceFuncPtr[x_new(i_inProcess[k] - 1, j_inProcess[k], height)]);
			}
			if (j_inProcess[k] == 0 || j_inProcess[k] == (width - 1) ) {
				y = 0;
			}
			else {
				y = min(distanceFuncPtr[x_new(i_inProcess[k], j_inProcess[k] + 1, height)], distanceFuncPtr[x_new(i_inProcess[k], j_inProcess[k] - 1, height)]);
			}

			coef = space / speed; 

			a = 2;
			if (x == INFINITY) {
				a--; x = 0;
			}
			if (y == INFINITY) {
				a--; y = 0;
			}
			b = -2 * (x + y);
			c = pow(x, 2) + pow(y, 2) - pow(coef, 2);
			delta = pow(b, 2) - 4 * a * c;
			
			//We need to select the largest quadratic solution
			//J.A Sethian, A Fast Marching Level Set method for Monotonically advancing fronts, 1995, page 8.
			//link to article ---> http://ugweb.cs.ualberta.ca/~vis/courses/CompVis/readings/modelrec/sethian95fastlev.pdf 
			if (delta >= 0) {
				dist = (dataType)( ( (- b + sqrt(delta)) / (2 * a) ) );
				//cout << "\nT1  : " << (dataType)(( (- b + sqrt(delta)) / (2 * a)) ) << endl;
				//cout << "\nT2  : " << (dataType)(( (- b - sqrt(delta)) / (2 * a)) ) << endl;
			}
			else {
				dist = (dataType)(min(x, y) + pow(coef, 2));
			}

			/*if (sqrt(2) * coef > abs(x - y)) {
				dist = 0.5 * (x + y) + 0.5 * sqrt(pow(x + y, 2) - 2 * (pow(x, 2) + pow(x, 2) - pow(coef, 2)));
			}
			else {
				dist = min(x, y) + pow(coef, 2);
			}*/

			tempTimeFunc.push_back(dist);
			
			//cout << "\ndistance  : " << dist << endl;

		}

		//Find the minimal solution  ---> same article page 12
		minSolution = INFINITY;
		for (k = 0; k < h_n; k++) {
			if (minSolution >= tempTimeFunc[k]) {
				minSolution = tempTimeFunc[k];
				iSol = k; iNew = i_inProcess[iSol];
				jSol = k; jNew = j_inProcess[jSol];
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
		}
		else {
			if (iSol == i_inProcess.size() - 1 || jSol == j_inProcess.size() - 1) {
				i_inProcess.erase(i_inProcess.end());
				j_inProcess.erase(j_inProcess.end());
			}
			else {
				i_inProcess.erase(i_inProcess.begin() + iSol);
				j_inProcess.erase(j_inProcess.begin() + jSol);
			}
		}
		tempTimeFunc.clear();

		//STEP 4
		//Find the neighbors of the processed pixel
		if (iNew == 0) {
			if (jNew == 0) {
				//East
				if (labelArray[x_new(iNew, jNew + 1, height)] == 3) {
					i_inProcess.push_back(iNew); j_inProcess.push_back(jNew + 1);
					labelArray[x_new(iNew, jNew + 1, height)] = 2;
				}
				//South
				if (labelArray[x_new(iNew + 1, jNew, height)] == 3) {
					i_inProcess.push_back(iNew + 1); j_inProcess.push_back(jNew);
					labelArray[x_new(iNew + 1, jNew, height)] = 2;
				}
			}
			else {
				if (jNew == width - 1) {
					//West
					if (labelArray[x_new(iNew, jNew - 1, height)] == 3) {
						i_inProcess.push_back(iNew); j_inProcess.push_back(jNew - 1);
						labelArray[x_new(iNew, jNew - 1, height)] = 2;
					}
					//South
					if (labelArray[x_new(iNew + 1, jNew, height)] == 3) {
						i_inProcess.push_back(iNew + 1); j_inProcess.push_back(jNew);
						labelArray[x_new(iNew + 1, jNew, height)] = 2;
					}
				}
				else {
					//West
					if (labelArray[x_new(iNew, jNew - 1, height)] == 3) {
						i_inProcess.push_back(iNew); j_inProcess.push_back(jNew - 1);
						labelArray[x_new(iNew, jNew - 1, height)] = 2;
					}
					//East
					if (labelArray[x_new(iNew, jNew + 1, height)] == 3) {
						i_inProcess.push_back(iNew); j_inProcess.push_back(jNew + 1);
						labelArray[x_new(iNew, jNew + 1, height)] = 2;
					}
					//South
					if (labelArray[x_new(iNew + 1, jNew, height)] == 3) {
						i_inProcess.push_back(iNew + 1); j_inProcess.push_back(jNew);
						labelArray[x_new(iNew + 1, jNew, height)] = 2;
					}
				}
			}
		}
		else {
			if (iNew == height - 1) {
				if (jNew == 0) {
					//North
					if (labelArray[x_new(iNew - 1, jNew, height)] == 3) {
						i_inProcess.push_back(iNew - 1); j_inProcess.push_back(jNew);
						labelArray[x_new(iNew - 1, jNew, height)] = 2;
					}
					//East
					if (labelArray[x_new(iNew, jNew + 1, height)] == 3) {
						i_inProcess.push_back(iNew); j_inProcess.push_back(jNew + 1);
						labelArray[x_new(iNew, jNew + 1, height)] = 2;
					}
				}
				else {
					if (jNew == width - 1) {
						//North
						if (labelArray[x_new(iNew - 1, jNew, height)] == 3) {
							i_inProcess.push_back(iNew - 1); j_inProcess.push_back(jNew);
							labelArray[x_new(iNew - 1, jNew, height)] = 2;
						}
						//West
						if (labelArray[x_new(iNew, jNew - 1, height)] == 3) {
							i_inProcess.push_back(iNew); j_inProcess.push_back(jNew - 1);
							labelArray[x_new(iNew, jNew - 1, height)] = 2;
						}
					}
					else {
						//West
						if (labelArray[x_new(iNew, jNew - 1, height)] == 3) {
							i_inProcess.push_back(iNew); j_inProcess.push_back(jNew - 1);
							labelArray[x_new(iNew, jNew - 1, height)] = 2;
						}
						//North
						if (labelArray[x_new(iNew - 1, jNew, height)] == 3) {
							i_inProcess.push_back(iNew - 1); j_inProcess.push_back(jNew);
							labelArray[x_new(iNew - 1, jNew, height)] = 2;
						}
						//East
						if (labelArray[x_new(iNew, jNew + 1, height)] == 3) {
							i_inProcess.push_back(iNew); j_inProcess.push_back(jNew + 1);
							labelArray[x_new(iNew, jNew + 1, height)] = 2;
						}
					}
				}
			}
			else {
				if (jNew == 0) {
					//North
					if (labelArray[x_new(iNew - 1, jNew, height)] == 3) {
						i_inProcess.push_back(iNew - 1); j_inProcess.push_back(jNew);
						labelArray[x_new(iNew - 1, jNew, height)] = 2;
					}
					//East
					if (labelArray[x_new(iNew, jNew + 1, height)] == 3) {
						i_inProcess.push_back(iNew); j_inProcess.push_back(jNew + 1);
						labelArray[x_new(iNew, jNew + 1, height)] = 2;
					}
					//South
					if (labelArray[x_new(iNew + 1, jNew, height)] == 3) {
						i_inProcess.push_back(iNew + 1); j_inProcess.push_back(jNew);
						labelArray[x_new(iNew + 1, jNew, height)] = 2;
					}
				}
				else {
					if (jNew == width - 1) {
						//North
						if (labelArray[x_new(iNew - 1, jNew, height)] == 3) {
							i_inProcess.push_back(iNew - 1); j_inProcess.push_back(jNew);
							labelArray[x_new(iNew - 1, jNew, height)] = 2;
						}
						//West
						if (labelArray[x_new(iNew, jNew - 1, height)] == 3) {
							i_inProcess.push_back(iNew); j_inProcess.push_back(jNew - 1);
							labelArray[x_new(iNew, jNew - 1, height)] = 2;
						}
						//South
						if (labelArray[x_new(iNew + 1, jNew, height)] == 3) {
							i_inProcess.push_back(iNew + 1); j_inProcess.push_back(jNew);
							labelArray[x_new(iNew + 1, jNew, height)] = 2;
						}
					}
					else {
						//North
						if (labelArray[x_new(iNew - 1, jNew, height)] == 3) {
							i_inProcess.push_back(iNew - 1); j_inProcess.push_back(jNew);
							labelArray[x_new(iNew - 1, jNew, height)] = 2;
						}
						/*else {
							if (labelArray[x_new(iNew - 1, jNew, imageHeight)] == 2) {
								if ((iNew - 1) == 0 || (iNew - 1) == (imageHeight - 1) ) {
									x = 0;
								}
								else {
									x = min(outputImage[x_new((iNew - 1) + 1, jNew, imageHeight)], outputImage[x_new((iNew - 1) - 1, jNew, imageHeight)]);
								}
								if (jNew == 0 || jNew == (imageWidth - 1)) {
									y = 0;
								}
								else {
									y = min(outputImage[x_new(iNew - 1, jNew + 1, imageHeight)], outputImage[x_new(iNew - 1, jNew - 1, imageHeight)]);
								}
								if (sqrt(2) * coef > abs(x - y)) {
									outputImage[x_new(iNew - 1, jNew, imageHeight)] = 0.5 * (x + y) + 0.5 * sqrt(pow(x + y, 2) - 2 * (pow(x, 2) + pow(x, 2) - pow(coef, 2)));
								}
								else {
									outputImage[x_new(iNew - 1, jNew, imageHeight)] = min(x, y) + pow(coef, 2);
								}
							}
						}*/
						//West
						if (labelArray[x_new(iNew, jNew - 1, height)] == 3) {
							i_inProcess.push_back(iNew); j_inProcess.push_back(jNew - 1);
							labelArray[x_new(iNew, jNew - 1, height)] = 2;
						}
						/*else {
							if (labelArray[x_new(iNew, jNew - 1, imageHeight)] == 2) {
								if (iNew == 0 || iNew == (imageHeight - 1)) {
									x = 0;
								}
								else {
									x = min(outputImage[x_new(iNew + 1, (jNew - 1), imageHeight)], outputImage[x_new(iNew - 1, (jNew - 1), imageHeight)]);
								}
								if ((jNew - 1) == 0 || (jNew - 1) == (imageWidth - 1)) {
									y = 0;
								}
								else {
									y = min(outputImage[x_new(iNew, (jNew - 1) + 1, imageHeight)], outputImage[x_new(iNew, (jNew - 1) - 1, imageHeight)]);
								}
								if (sqrt(2) * coef > abs(x - y)) {
									outputImage[x_new(iNew - 1, jNew, imageHeight)] = 0.5 * (x + y) + 0.5 * sqrt(pow(x + y, 2) - 2 * (pow(x, 2) + pow(x, 2) - pow(coef, 2)));
								}
								else {
									outputImage[x_new(iNew - 1, jNew, imageHeight)] = min(x, y) + pow(coef, 2);
								}
							}
						}*/
						//East
						if (labelArray[x_new(iNew, jNew + 1, height)] == 3) {
							i_inProcess.push_back(iNew); j_inProcess.push_back(jNew + 1);
							labelArray[x_new(iNew, jNew + 1, height)] = 2;
						}
						/*else {
							if (labelArray[x_new(iNew, jNew + 1, imageHeight)] == 2) {
								if (iNew == 0 || iNew == (imageHeight - 1)) {
									x = 0;
								}
								else {
									x = min(outputImage[x_new(iNew + 1, (jNew + 1), imageHeight)], outputImage[x_new(iNew - 1, (jNew + 1), imageHeight)]);
								}
								if ((jNew + 1) == 0 || (jNew + 1) == (imageWidth - 1) ) {
									y = 0;
								}
								else {
									y = min(outputImage[x_new(iNew, (jNew + 1) + 1, imageHeight)], outputImage[x_new(iNew, (jNew + 1) - 1, imageHeight)]);
								}
								if (sqrt(2) * coef > abs(x - y)) {
									outputImage[x_new(iNew - 1, jNew, imageHeight)] = 0.5 * (x + y) + 0.5 * sqrt(pow(x + y, 2) - 2 * (pow(x, 2) + pow(x, 2) - pow(coef, 2)));
								}
								else {
									outputImage[x_new(iNew - 1, jNew, imageHeight)] = min(x, y) + pow(coef, 2);
								}
							}
						}*/
						//South
						if (labelArray[x_new(iNew + 1, jNew, height)] == 3) {
							i_inProcess.push_back(iNew + 1); j_inProcess.push_back(jNew);
							labelArray[x_new(iNew + 1, jNew, height)] = 2;
						}
						/*else {
							if (labelArray[x_new(iNew + 1, jNew, imageHeight)] == 2) {
								if ((iNew + 1) == 0 || (iNew + 1) == (imageHeight - 1)) {
									x = 0;
								}
								else {
									x = min(outputImage[x_new((iNew + 1) + 1, jNew, imageHeight)], outputImage[x_new((iNew + 1) - 1, jNew, imageHeight)]);
								}
								if (jNew == 0 || jNew == (imageWidth - 1)) {
									y = 0;
								}
								else {
									y = min(outputImage[x_new((iNew + 1), jNew + 1, imageHeight)], outputImage[x_new((iNew + 1), jNew - 1, imageHeight)]);
								}
								if (sqrt(2) * coef > abs(x - y)) {
									outputImage[x_new(iNew - 1, jNew, imageHeight)] = 0.5 * (x + y) + 0.5 * sqrt(pow(x + y, 2) - 2 * (pow(x, 2) + pow(x, 2) - pow(coef, 2)));
								}
								else {
									outputImage[x_new(iNew - 1, jNew, imageHeight)] = min(x, y) + pow(coef, 2);
								}
							}
						}*/
					}
				}
			}
		}

		//Empty the solution list after finding the good value
		//tempTimeFunc.clear();
		/*while (tempTimeFunc.size() > 0) {
			tempTimeFunc.pop_back();
		}*/

	}

	for (k = 0; k < dim2D; k++) {
		if (distanceFuncPtr[k] == INFINITY) {
			distanceFuncPtr[k] = 0;
			cpt++;
		}
	}
	cout << "\n" << cpt << " pixels with distance = infinity : " << endl;

	free(labelArray);

	return true;
}
