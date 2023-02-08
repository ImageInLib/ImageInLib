#include <iostream>
#include <climits>
#include <crtdbg.h>
#include <corecrt_malloc.h>
#include<cmath>
#include <omp.h>
#include<vector>
#include "distanceForPathFinding.h"

using namespace std;

bool fastMarching2d(dataType* inputImagePtr, dataType* outputImage, const size_t imageHeight, const size_t imageWidth, Point2D* seedPoints)
{
	// The image should be threshoded before (or at least define the background pixels).

	short* labelArray = (short*)malloc(imageHeight * imageWidth * sizeof(short));

	if (inputImagePtr == NULL || outputImage == NULL || seedPoints == NULL || labelArray == NULL) {
		return false;
	}

	size_t i = 0, j = 0, dim2D = imageWidth * imageHeight;
	vector<size_t> i_Processed, i_inProcess;
	vector<size_t> j_Processed, j_inProcess;
	dataType x, y, timeFunction = 0, speed = 1.0, space = 1.0, minSolution = 0, coef = 0, dist = 0;
	size_t nbNeighborsFound = 0;
	size_t h_n = 0, w_n = 0, iSol = 0, jSol = 0, iNew = 0, jNew = 0;
	vector<dataType> tempTimeFunc;
	dataType a, b, c, delta;

	//STEP 1
	//Initialization : distance is set to infinity
	//In labelAray we have : 1 ---> already processed, 2 ---> in process and 3 ---> not processed
	for (i = 0; i < dim2D; i++) {
		outputImage[i] = INFINITY;
		labelArray[i] = 3;
	}

	//--------------------End of STEP 1 -----------------------------------

	//STEP 2
	//Add the neighbors of the seed point in the list of pixels to be processed
	i = seedPoints->x; j = seedPoints->y;
	outputImage[x_new(i, j, imageHeight)] = 0;
	i_Processed.push_back(i); j_Processed.push_back(j);

	if (i == 0) {
		if (j == 0) {
			//East
			if (labelArray[x_new(i, j + 1, imageHeight)] == 3) {
				i_inProcess.push_back(i); j_inProcess.push_back(j + 1);
				labelArray[x_new(i, j + 1, imageHeight)] = 2;
				nbNeighborsFound++;
			}
			//South
			if (labelArray[x_new(i + 1, j, imageHeight)] == 3) {
				i_inProcess.push_back(i + 1); j_inProcess.push_back(j);
				labelArray[x_new(i + 1, j, imageHeight)] = 2;
				nbNeighborsFound++;
			}
		}
		else {
			if (j == imageWidth - 1) {
				//West
				if (labelArray[x_new(i, j - 1, imageHeight)] == 3) {
					i_inProcess.push_back(i); j_inProcess.push_back(j - 1);
					labelArray[x_new(i, j - 1, imageHeight)] = 2;
					nbNeighborsFound++;
				}
				//South
				if (labelArray[x_new(i + 1, j, imageHeight)] == 3) {
					i_inProcess.push_back(i + 1); j_inProcess.push_back(j);
					labelArray[x_new(i + 1, j, imageHeight)] = 2;
					nbNeighborsFound++;
				}
			}
			else {
				//West
				if (labelArray[x_new(i, j - 1, imageHeight)] == 3) {
					i_inProcess.push_back(i); j_inProcess.push_back(j - 1);
					labelArray[x_new(i, j - 1, imageHeight)] = 2;
					nbNeighborsFound++;
				}
				//East
				if (labelArray[x_new(i, j + 1, imageHeight)] == 3) {
					i_inProcess.push_back(i); j_inProcess.push_back(j + 1);
					labelArray[x_new(i, j + 1, imageHeight)] = 2;
					nbNeighborsFound++;
				}
				//South
				if (labelArray[x_new(i + 1, j, imageHeight)] == 3) {
					i_inProcess.push_back(i + 1); j_inProcess.push_back(j);
					labelArray[x_new(i + 1, j, imageHeight)] = 2;
					nbNeighborsFound++;
				}
			}
		}
	}
	else {
		if (i == imageHeight - 1) {
			if (j == 0) {
				//North
				if (labelArray[x_new(i - 1, j, imageHeight)] == 3) {
					i_inProcess.push_back(i - 1); j_inProcess.push_back(j);
					labelArray[x_new(i - 1, j, imageHeight)] = 2;
					nbNeighborsFound++;
				}
				//East
				if (labelArray[x_new(i, j + 1, imageHeight)] == 3) {
					i_inProcess.push_back(i); j_inProcess.push_back(j + 1);
					labelArray[x_new(i, j + 1, imageHeight)] = 2;
					nbNeighborsFound++;
				}
			}
			else {
				if (j == imageWidth - 1) {
					//North
					if (labelArray[x_new(i - 1, j, imageHeight)] == 3) {
						i_inProcess.push_back(i - 1); j_inProcess.push_back(j);
						labelArray[x_new(i - 1, j, imageHeight)] = 2;
						nbNeighborsFound++;
					}
					//West
					if (labelArray[x_new(i, j - 1, imageHeight)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(j - 1);
						labelArray[x_new(i, j - 1, imageHeight)] = 2;
						nbNeighborsFound++;
					}
				}
				else {
					//West
					if (labelArray[x_new(i, j - 1, imageHeight)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(j - 1);
						labelArray[x_new(i, j - 1, imageHeight)] = 2;
						nbNeighborsFound++;
					}
					//North
					if (labelArray[x_new(i - 1, j, imageHeight)] == 3) {
						i_inProcess.push_back(i - 1); j_inProcess.push_back(j);
						labelArray[x_new(i - 1, j, imageHeight)] = 2;
						nbNeighborsFound++;
					}
					//East
					if (labelArray[x_new(i, j + 1, imageHeight)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(j + 1);
						labelArray[x_new(i, j + 1, imageHeight)] = 2;
						nbNeighborsFound++;
					}
				}
			}
		}
		else {
			if (j == 0) {
				//North
				if (labelArray[x_new(i - 1, j, imageHeight)] == 3) {
					i_inProcess.push_back(i - 1); j_inProcess.push_back(j);
					labelArray[x_new(i - 1, j, imageHeight)] = 2;
					nbNeighborsFound++;
				}
				//East
				if (labelArray[x_new(i, j + 1, imageHeight)] == 3) {
					i_inProcess.push_back(i); j_inProcess.push_back(j + 1);
					labelArray[x_new(i, j + 1, imageHeight)] = 2;
					nbNeighborsFound++;
				}
				//South
				if (labelArray[x_new(i + 1, j, imageHeight)] == 3) {
					i_inProcess.push_back(i + 1); j_inProcess.push_back(j);
					labelArray[x_new(i + 1, j, imageHeight)] = 2;
					nbNeighborsFound++;
				}
			}
			else {
				if (j == imageWidth - 1) {
					//North
					if (labelArray[x_new(i - 1, j, imageHeight)] == 3) {
						i_inProcess.push_back(i - 1); j_inProcess.push_back(j);
						labelArray[x_new(i - 1, j, imageHeight)] = 2;
						nbNeighborsFound++;
					}
					//West
					if (labelArray[x_new(i, j - 1, imageHeight)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(j - 1);
						labelArray[x_new(i, j - 1, imageHeight)] = 2;
						nbNeighborsFound++;
					}
					//South
					if (labelArray[x_new(i + 1, j, imageHeight)] == 3) {
						i_inProcess.push_back(i + 1); j_inProcess.push_back(j);
						labelArray[x_new(i + 1, j, imageHeight)] = 2;
						nbNeighborsFound++;
					}
				}
				else {
					//North
					if (labelArray[x_new(i - 1, j, imageHeight)] == 3) {
						i_inProcess.push_back(i - 1); j_inProcess.push_back(j);
						labelArray[x_new(i - 1, j, imageHeight)] = 2;
						nbNeighborsFound++;
					}
					//West
					if (labelArray[x_new(i, j - 1, imageHeight)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(j - 1);
						labelArray[x_new(i, j - 1, imageHeight)] = 2;
						nbNeighborsFound++;
					}
					//East
					if (labelArray[x_new(i, j + 1, imageHeight)] == 3) {
						i_inProcess.push_back(i); j_inProcess.push_back(j + 1);
						labelArray[x_new(i, j + 1, imageHeight)] = 2;
						nbNeighborsFound++;
					}
					//South
					if (labelArray[x_new(i + 1, j, imageHeight)] == 3) {
						i_inProcess.push_back(i + 1); j_inProcess.push_back(j);
						labelArray[x_new(i + 1, j, imageHeight)] = 2;
						nbNeighborsFound++;
					}
				}
			}
		}
	}

	//Update the label of seed point as processed
	labelArray[x_new(i, j, imageHeight)] = 1;

	cout << "\nNumber of Neigbors found : " << nbNeighborsFound << endl;
	//---------------------End of STEP 2 -------------------------------------

	//STEP 3
	while (i_inProcess.size() != 0 || j_inProcess.size() != 0 ) {

		h_n = i_inProcess.size();
		//w_n = j_inProcess.size(); // they have the same value

		//Solve the Eikonale equation for all the neighbors in the stack
		for (i = 0; i < h_n; i++) {
			j = i;

			if (i_inProcess[i] == 0 || i_inProcess[i] == imageHeight - 1) {
				x = 0;
			}
			else {
				x = min(outputImage[x_new(i_inProcess[i] + 1, j_inProcess[j], imageHeight)], outputImage[x_new(i_inProcess[i] - 1, j_inProcess[j], imageHeight)]);
			}

			if (j_inProcess[j] == 0 || j_inProcess[j] == imageWidth - 1) {
				y = 0;
			}
			else {
				y = min(outputImage[x_new(i_inProcess[i], j_inProcess[j] + 1, imageHeight)], outputImage[x_new(i_inProcess[i], j_inProcess[j] - 1, imageHeight)]);
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

			if (delta >= 0) {
				dist = (dataType)( ( (-b + sqrt(delta)) / (2 * a) ) );
			}
			else {
				dist = (dataType)(min(x, y) + pow(coef, 2));
			}
			tempTimeFunc.push_back(dist);

			//cout << "\ndistance  : " << dist << endl;

		}

		//Find the minimal solution
		minSolution = INFINITY;
		for (i = 0; i < h_n; i++) {
			if (minSolution > tempTimeFunc[i]) {
				minSolution = tempTimeFunc[i];
				iSol = i; iNew = i_inProcess[iSol];
				jSol = i; jNew = j_inProcess[jSol];
			}
		}

		//Set the distance to the processed pixel
		outputImage[x_new(iNew, jNew, imageHeight)] = minSolution;
		labelArray[x_new(iNew, jNew, imageHeight)] = 1;
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

		//STEP 4
		//Find the neighbors of the processed pixel
		if (iNew == 0) {
			if (jNew == 0) {
				//East
				if (labelArray[x_new(iNew, jNew + 1, imageHeight)] == 3) {
					i_inProcess.push_back(iNew); j_inProcess.push_back(jNew + 1);
					labelArray[x_new(iNew, jNew + 1, imageHeight)] = 2;
				}
				//South
				if (labelArray[x_new(iNew + 1, jNew, imageHeight)] == 3) {
					i_inProcess.push_back(iNew + 1); j_inProcess.push_back(jNew);
					labelArray[x_new(iNew + 1, jNew, imageHeight)] = 2;
				}
			}
			else {
				if (jNew == imageWidth - 1) {
					//West
					if (labelArray[x_new(iNew, jNew - 1, imageHeight)] == 3) {
						i_inProcess.push_back(iNew); j_inProcess.push_back(jNew - 1);
						labelArray[x_new(iNew, jNew - 1, imageHeight)] = 2;
					}
					//South
					if (labelArray[x_new(iNew + 1, jNew, imageHeight)] == 3) {
						i_inProcess.push_back(iNew + 1); j_inProcess.push_back(jNew);
						labelArray[x_new(iNew + 1, jNew, imageHeight)] = 2;
					}
				}
				else {
					//West
					if (labelArray[x_new(iNew, jNew - 1, imageHeight)] == 3) {
						i_inProcess.push_back(iNew); j_inProcess.push_back(jNew - 1);
						labelArray[x_new(iNew, jNew - 1, imageHeight)] = 2;
					}
					//East
					if (labelArray[x_new(iNew, jNew + 1, imageHeight)] == 3) {
						i_inProcess.push_back(iNew); j_inProcess.push_back(jNew + 1);
						labelArray[x_new(iNew, jNew + 1, imageHeight)] = 2;
					}
					//South
					if (labelArray[x_new(iNew + 1, jNew, imageHeight)] == 3) {
						i_inProcess.push_back(iNew + 1); j_inProcess.push_back(jNew);
						labelArray[x_new(iNew + 1, jNew, imageHeight)] = 2;
					}
				}
			}
		}
		else {
			if (iNew == imageHeight - 1) {
				if (jNew == 0) {
					//North
					if (labelArray[x_new(iNew - 1, jNew, imageHeight)] == 3) {
						i_inProcess.push_back(iNew - 1); j_inProcess.push_back(jNew);
						labelArray[x_new(iNew - 1, jNew, imageHeight)] = 2;
					}
					//East
					if (labelArray[x_new(iNew, jNew + 1, imageHeight)] == 3) {
						i_inProcess.push_back(iNew); j_inProcess.push_back(jNew + 1);
						labelArray[x_new(iNew, jNew + 1, imageHeight)] = 2;
					}
				}
				else {
					if (jNew == imageWidth - 1) {
						//North
						if (labelArray[x_new(iNew - 1, jNew, imageHeight)] == 3) {
							i_inProcess.push_back(iNew - 1); j_inProcess.push_back(jNew);
							labelArray[x_new(iNew - 1, jNew, imageHeight)] = 2;
						}
						//West
						if (labelArray[x_new(iNew, jNew - 1, imageHeight)] == 3) {
							i_inProcess.push_back(iNew); j_inProcess.push_back(jNew - 1);
							labelArray[x_new(iNew, jNew - 1, imageHeight)] = 2;
						}
					}
					else {
						//West
						if (labelArray[x_new(iNew, jNew - 1, imageHeight)] == 3) {
							i_inProcess.push_back(iNew); j_inProcess.push_back(jNew - 1);
							labelArray[x_new(iNew, jNew - 1, imageHeight)] = 2;
						}
						//North
						if (labelArray[x_new(iNew - 1, jNew, imageHeight)] == 3) {
							i_inProcess.push_back(iNew - 1); j_inProcess.push_back(jNew);
							labelArray[x_new(iNew - 1, jNew, imageHeight)] = 2;
						}
						//East
						if (labelArray[x_new(iNew, jNew + 1, imageHeight)] == 3) {
							i_inProcess.push_back(iNew); j_inProcess.push_back(jNew + 1);
							labelArray[x_new(iNew, jNew + 1, imageHeight)] = 2;
						}
					}
				}
			}
			else {
				if (jNew == 0) {
					//North
					if (labelArray[x_new(iNew - 1, jNew, imageHeight)] == 3) {
						i_inProcess.push_back(iNew - 1); j_inProcess.push_back(jNew);
						labelArray[x_new(iNew - 1, jNew, imageHeight)] = 2;
					}
					//East
					if (labelArray[x_new(iNew, jNew + 1, imageHeight)] == 3) {
						i_inProcess.push_back(iNew); j_inProcess.push_back(jNew + 1);
						labelArray[x_new(iNew, jNew + 1, imageHeight)] = 2;
					}
					//South
					if (labelArray[x_new(iNew + 1, jNew, imageHeight)] == 3) {
						i_inProcess.push_back(iNew + 1); j_inProcess.push_back(jNew);
						labelArray[x_new(iNew + 1, jNew, imageHeight)] = 2;
					}
				}
				else {
					if (jNew == imageWidth - 1) {
						//North
						if (labelArray[x_new(iNew - 1, jNew, imageHeight)] == 3) {
							i_inProcess.push_back(iNew - 1); j_inProcess.push_back(jNew);
							labelArray[x_new(iNew - 1, jNew, imageHeight)] = 2;
						}
						//West
						if (labelArray[x_new(iNew, jNew - 1, imageHeight)] == 3) {
							i_inProcess.push_back(iNew); j_inProcess.push_back(jNew - 1);
							labelArray[x_new(iNew, jNew - 1, imageHeight)] = 2;
						}
						//South
						if (labelArray[x_new(iNew + 1, jNew, imageHeight)] == 3) {
							i_inProcess.push_back(iNew + 1); j_inProcess.push_back(jNew);
							labelArray[x_new(iNew + 1, jNew, imageHeight)] = 2;
						}
					}
					else {
						//North
						if (labelArray[x_new(iNew - 1, jNew, imageHeight)] == 3) {
							i_inProcess.push_back(iNew - 1); j_inProcess.push_back(jNew);
							labelArray[x_new(iNew - 1, jNew, imageHeight)] = 2;
						}
						//West
						if (labelArray[x_new(iNew, jNew - 1, imageHeight)] == 3) {
							i_inProcess.push_back(iNew); j_inProcess.push_back(jNew - 1);
							labelArray[x_new(iNew, jNew - 1, imageHeight)] = 2;
						}
						//East
						if (labelArray[x_new(iNew, jNew + 1, imageHeight)] == 3) {
							i_inProcess.push_back(iNew); j_inProcess.push_back(jNew + 1);
							labelArray[x_new(iNew, jNew + 1, imageHeight)] = 2;
						}
						//South
						if (labelArray[x_new(iNew + 1, jNew, imageHeight)] == 3) {
							i_inProcess.push_back(iNew + 1); j_inProcess.push_back(jNew);
							labelArray[x_new(iNew + 1, jNew, imageHeight)] = 2;
						}
					}
				}
			}
		}
	}

	short cpt = 0;
	for (i = 0; i < dim2D; i++) {
		if (outputImage[i] == INFINITY) {
			outputImage[i] = 0;
			cpt++;
		}
	}
	cout << "\n" << cpt << " pixels with distance = infinity : " << endl;

	free(labelArray);
	return true;
}
