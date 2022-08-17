#include <iostream>
#include <climits>
#include <crtdbg.h>
#include <corecrt_malloc.h>
#include <omp.h>
#include<vector>
#include "common_functions.h"
#include "Labelling.h"


#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

using namespace std;

//for 2D image using 2D arrays

bool labelling2D(int** imageDataPtr, int** segmentedImage, bool** statusArray, const size_t xDim, const size_t yDim, int object) {

	if (imageDataPtr == NULL)
		return false;
	if (segmentedImage == NULL)
		return false;
	if (statusArray == NULL)
		return false;

	vector<int> iStack, jStack, iTmpStack, jTmpStack;
	int i = 0, j = 0, iNew = 0, jNew = 0, k = 0, xd = 0, label = 1;

	//statusArray has false everywhere at the beginning
	//false---> non-processed and true---> processed
	//segmentedImage has 0 everywhere at the begining
	for (i = 0; i < xDim; i++) {
		for (j = 0; j < yDim; j++) {

			if (statusArray[i][j] == false) {
				if (imageDataPtr[i][j] == object) {
					//top neighbor
					if (i > 0 && statusArray[i - 1][j] == false) {
						if (imageDataPtr[i - 1][j] == object) {
							//element is in region add its coordinates in stacks
							iStack.push_back(i - 1);
							jStack.push_back(j);
						}
						else {
							//it's not object, so, update its status
							statusArray[i - 1][j] = true;
						}
					}
					//right neighbor
					if (j < yDim - 1 && statusArray[i][j + 1] == false) {
						if (imageDataPtr[i][j + 1] == object) {
							//element is in region add its coordinates in stacks
							iStack.push_back(i);
							jStack.push_back(j + 1);
						}
						else {
							statusArray[i][j + 1] = true;
						}
					}
					//bottom neighbor
					if (i < xDim - 1 && statusArray[i + 1][j] == false) {
						if (imageDataPtr[i + 1][j] == object) {
							//if element is in region add its coodinates in stacks
							iStack.push_back(i + 1);
							jStack.push_back(j);
						}
						else {
							statusArray[i + 1][j] = true;
						}
					}
					//left neighbor
					if (j > 0 && statusArray[i][j - 1] == false) {
						if (imageDataPtr[i][j - 1] == object) {
							//element is in region add its coodinates in stacks
							iStack.push_back(i);
							jStack.push_back(j - 1);
						}
						else {
							statusArray[i][j - 1] = true;
						}
					}
					//after checking all neighbors of current element
					//we give its label, and update its status
					segmentedImage[i][j] = label;
					statusArray[i][j] = true;
					//Now we check neighbors of its neighbors saved in the stacks
					//we start by the last added element in stacks
					//If there is no neighbor iStack.size() = jStack.size() = 0, and the while loop will no be ran
					while (iStack.size() > 0 && jStack.size() > 0) { //One is enought because they have same size
						iNew = iStack.size() - 1;
						jNew = jStack.size() - 1;
						//top neighbor
						if (iStack[iNew] > 0 && statusArray[ iStack[iNew] - 1 ][ jStack[jNew] ] == false) {
							if (imageDataPtr[iStack[iNew] - 1][jStack[jNew]] == object) {
								//If the element is in the current region, then save its coordinates in temporary stacks
								iTmpStack.push_back(iStack[iNew] - 1);
								jTmpStack.push_back(jStack[jNew]);
							}
							else {
								//Not in region, update status and go to the next neighbor
								statusArray[iStack[iNew] - 1][jStack[jNew]] = true;
							}
						}
						//right neighbor
						if (jStack[jNew] < yDim - 1 && statusArray[ iStack[iNew] ][ jStack[jNew] + 1 ] == false) {
							if (imageDataPtr[ iStack[iNew] ][ jStack[jNew] + 1 ] == object) {
								//If the element is in the current region, then save its coordinates in temporary stacks
								iTmpStack.push_back(iStack[iNew]);
								jTmpStack.push_back(jStack[jNew] + 1);
							}
							else {
								statusArray[ iStack[iNew] ][ jStack[jNew] + 1 ] = true;
							}
						}
						//bottom neighbor
						if (iStack[iNew] < xDim - 1 && statusArray[ iStack[iNew] + 1 ][ jStack[jNew] ] == false) {
							if (imageDataPtr[ iStack[iNew] + 1 ][ jStack[jNew] ] == object) {
								//If the element is in the current region, then save its coordinates in temporary stacks
								iTmpStack.push_back(iStack[iNew] + 1);
								jTmpStack.push_back(jStack[jNew]);
							}
							else {
								statusArray[ iStack[iNew] + 1 ][ jStack[jNew] ] = true;
							}
						}
						//left neighbor
						if (jStack[jNew] > 0 && statusArray[ iStack[iNew] ][ jStack[jNew] - 1 ] == false) {
							if (imageDataPtr[ iStack[iNew] ][ jStack[jNew] - 1 ] == object) {
								//If the element is in the current region, then save its coordinates in temporary stacks
								iTmpStack.push_back(iStack[iNew]);
								jTmpStack.push_back(jStack[jNew] - 1);
							}
							else {
								statusArray[ iStack[iNew] ][ jStack[jNew] - 1 ] = true;
							}
						}
						//updating of processed element befor removal
						segmentedImage[ iStack[iNew] ][ jStack[jNew] ] = label;
						statusArray[ iStack[iNew] ][ jStack[jNew] ] = true;
						//Remove the processed element of the initial stacks
						iStack.pop_back();
						jStack.pop_back();
						//Add new found neighbors in initial stacks
						//we can do it once because they have the same size
						for (k = 0; k < iTmpStack.size(); k++) {
							//if any neighbors have been found iTmpStack.size() = 0
							//and nothing will happen
							iStack.push_back(iTmpStack[k]);
						}
						for (k = 0; k < jTmpStack.size(); k++) {
							//if any neighbors have been found jTmpStack.size() = 0
							//and nothing will happen
							jStack.push_back(jTmpStack[k]);
						}
						//empty the temporary stacks
						//we can do it once because they have the same size
						while (iTmpStack.size() > 0) {
							iTmpStack.pop_back();
						}
						while (jTmpStack.size() > 0) {
							jTmpStack.pop_back();
						}
					}
					label++;
				}
				else {
					statusArray[i][j] = true;
				}
			}
		}
	}
	return true;
}

//Erosion ---> Shrink objects and remove boundaries pixels
bool erosion2D(int** imageDataPtr, const size_t xDim, const size_t yDim, int background, int object) {
	
	if (imageDataPtr == NULL)
		return false;

	size_t i, j, k;

	for (i = 0; i < xDim; i++) {
		for (j = 0; j < yDim; j++) {
			k = 0;
			if (i == xDim - 1) {
				if (j == yDim - 1) {
					continue;
				}
				else {
					if (imageDataPtr[i][j + 1] == background) {
						imageDataPtr[i][j] = background;
					}
				}
			}
			else {
				if (j == yDim - 1) {
					if (imageDataPtr[i + 1][j] == background) {
						imageDataPtr[i][j] = background;
					}
				}
				else {
					if (imageDataPtr[i][j + 1] == object) {
						k++;
					}

					if (imageDataPtr[i + 1][j] == object) {
						k++;
					}

					if (imageDataPtr[i + 1][j + 1] == object) {
						k++;
					}

					if (k == 3) {
						imageDataPtr[i][j] = object;
					}
					else {
						imageDataPtr[i][j] = background;
					}
				}
			}
		}
	}
	return true;
}

//Statistics after labelling
bool labellingStats(int** segmentedImage, int* CountingArray, const size_t xDim, const size_t yDim, const char* pathPtr) {
	if (segmentedImage == NULL) return false;
	if (CountingArray == NULL) return false;

	int i, j;

	FILE* file;
	if (fopen_s(&file, pathPtr, "w") != 0) {
		printf("Enable to open");
		return false;
	}

	int numObjectCells = 0;
	for (i = 0; i < xDim; i++) {
		for (j = 0; j < yDim; j++) {
			if (segmentedImage[i][j] > 0) {
				numObjectCells++;
			}
		}
	}
	int numOfRegions = 0;
	for (i = 0; i < numObjectCells; i++) {
		if (CountingArray[i] > 0) {
			numOfRegions++;
		}
	}

	fprintf(file, "       2D image \n");
	fprintf(file, " Number of Object Cells = %d \n", numObjectCells);
	fprintf(file, " Number of Regions = %d \n", numOfRegions);
	fprintf(file, "#############################################################\n");
	fprintf(file, "\n\n");
	fprintf(file, "Region Label        Count");
	for (i = 1; i < numObjectCells; i++) {
		if (CountingArray[i] > 0) {
			fprintf(file, "\n");
			fprintf(file, "   %d                 %d", i, CountingArray[i]);
		}
	}
	fclose(file);

	return true;
}


//========================
bool initialization2dArray(int** imageDataPtr, const size_t xDim, const size_t yDim, int value)
{
	if (imageDataPtr == NULL) {
		return false;
	}

	for (int i = 0; i < xDim; i++) {
		for (int j = 0; j < yDim; j++) {
			imageDataPtr[i][j] = value;
		}
	}
	return true;
}

//pixel values are 0 and 255
bool minorRegionRemoval(int** imageDataPtr, int** segmentedImage, int* CountingArray, const size_t xDim, const size_t yDim, int size)
{
	int i, j;
	if (imageDataPtr == NULL) {
		return false;
	}
	if (segmentedImage == NULL) {
		return false;
	}
	if (CountingArray == NULL) {
		return false;
	}

	for (i = 0; i < xDim; i++) {
		for (j = 0; j < yDim; j++) {
			if (segmentedImage[i][j] > 0) {
				if (CountingArray[segmentedImage[i][j]] < size) {
					imageDataPtr[i][j] = 255;
				}
			}
		}
	}

	return true;
}

bool rescalingTo2D(int** imageDataPtr, int xDim, int yDim, int minNew, int maxNew)
{
	int i, j;
	double factor;
	int minOld = INT_MAX;
	int maxOld = INT_MIN;

	for (i = 0; i < xDim; i++) {
		for (j = 0; j < yDim; j++) {
			if (imageDataPtr[i][j] > maxOld) {
				maxOld = imageDataPtr[i][j];
			}
			if (imageDataPtr[i][j] < minOld) {
				minOld = imageDataPtr[i][j];
			}
		}
	}

	factor = (double)minNew / (maxOld - minOld);

	for (i = 0; i < xDim; i++) {
		for (j = 0; j < yDim; j++) {
			imageDataPtr[i][j] = (int)(factor * (imageDataPtr[i][j] - maxOld) + maxNew + 0.5);
		}
	}

	return true;
}

bool sortArray(int* valuePtr, int sizeArray) {
	int i, j, max = 0, ind = 0, ech = 0;

	if (valuePtr == NULL)
		return false;

	for (i = 0; i < sizeArray; i++) {
		max = valuePtr[i];
		for (j = i; j < sizeArray; j++) {
			if (valuePtr[j] > max) {
				max = valuePtr[j];
				ind = j;
			}
		}
		ech = valuePtr[i];
		valuePtr[i] = max;
		valuePtr[ind] = ech;
	}

	return true;
}
//=========================


//for 3D images
//New algorithm
bool labelling3D(dataType** imageDataPtr, int** segmentedImage, bool** statusArray, const size_t xDim, const size_t yDim, const size_t zDim, dataType object) {

	if (imageDataPtr == NULL)
		return false;
	if (segmentedImage == NULL)
		return false;
	if (statusArray == NULL)
		return false;

	vector<size_t> iStack, jStack, kStack, iTmpStack, jTmpStack, kTmpStack;
	size_t i = 0, j = 0, k = 0, iNew = 0, jNew = 0, kNew = 0, n = 0;
	int label = 1;

	//statusArray has false everywhere in the beginning
	//false---> non-processed and true---> processed
	//segmentedImage has 0 everywhere
	for (k = 0; k < zDim; k++) {
		for (i = 0; i < xDim; i++) {
			for (j = 0; j < yDim; j++) {

				if (statusArray[k][x_new(i, j, xDim)] == false) {
					if (imageDataPtr[k][x_new(i, j, xDim)] == object) {
						//top neighbor
						if (k > 0 && statusArray[k - 1][x_new(i, j, xDim)] == false) {
							if (imageDataPtr[k - 1][x_new(i, j, xDim)] == object) {
								//if element is in region, then add its coordinates in stacks
								iStack.push_back(i);
								jStack.push_back(j);
								kStack.push_back(k - 1);
							}
							else {
								//it's not object, so, update its status
								statusArray[k - 1][x_new(i, j, xDim)] = true;
							}
						}
						//down neighbor
						if (k < zDim - 1 && statusArray[k + 1][x_new(i, j, xDim)] == false) {
							if (imageDataPtr[k + 1][x_new(i, j, xDim)] == object) {
								//if element is in region add its coordinates in stacks
								iStack.push_back(i);
								jStack.push_back(j);
								kStack.push_back(k + 1);
							}
							else {
								//it's not object, so, update its status
								statusArray[k + 1][x_new(i, j, xDim)] = true;
							}
						}
						//left neighbor
						if (i > 0 && statusArray[k][x_new(i - 1, j, xDim)] == false) {
							if (imageDataPtr[k][x_new(i - 1, j, xDim)] == object) {
								//if element is in region add its coordinates in stacks
								iStack.push_back(i - 1);
								jStack.push_back(j);
								kStack.push_back(k);
							}
							else {
								statusArray[k][x_new(i - 1, j, xDim)] = true;
							}
						}
						//right neighbor
						if (i < xDim - 1 && statusArray[k][x_new(i + 1, j, xDim)] == false) {
							if (imageDataPtr[k][x_new(i + 1, j, xDim)] == object) {
								//if element is in region add its coordinates in stacks
								iStack.push_back(i + 1);
								jStack.push_back(j);
								kStack.push_back(k);
							}
							else {
								statusArray[k][x_new(i + 1, j, xDim)] = true;
							}
						}
						//front neighbor
						if (j > 0 && statusArray[k][x_new(i, j - 1, xDim)] == false) {
							if (imageDataPtr[k][x_new(i, j - 1, xDim)] == object) {
								//if element is in region add its coodinates in stacks
								iStack.push_back(i);
								jStack.push_back(j - 1);
								kStack.push_back(k);
							}
							else {
								statusArray[k][x_new(i, j - 1, xDim)] = true;
							}
						}
						//behind neighbor
						if (j < yDim - 1 && statusArray[k][x_new(i, j + 1, xDim)] == false) {
							if (imageDataPtr[k][x_new(i, j + 1, xDim)] == object) {
								//if element is in region add its coodinates in stacks
								iStack.push_back(i);
								jStack.push_back(j + 1);
								kStack.push_back(k);
							}
							else {
								statusArray[k][x_new(i, j + 1, xDim)] = true;
							}
						}
						//after checking all neighbors of current element
						//I give its label, and update its status
						segmentedImage[k][x_new(i, j, xDim)] = label;
						statusArray[k][x_new(i, j, xDim)] = true;
						//Now I check neighbors of its neighbors saved in the stacks
						//I start by the last added element in stacks
						//If there is no neighbor iStack.size() = jStack.size() = kStack.size() = 0 , and the while loop will no be ran
						while (iStack.size() > 0 && jStack.size() > 0 && kStack.size() > 0) {
							//One is enought because they have same size
							iNew = iStack.size() - 1;
							jNew = jStack.size() - 1;
							kNew = kStack.size() - 1;
							//top neighbor
							if (kStack[kNew] > 0 && statusArray[kStack[kNew] - 1][x_new(iStack[iNew], jStack[jNew], xDim)] == false) {
								if (imageDataPtr[kStack[kNew] - 1][x_new(iStack[iNew], jStack[jNew], xDim)] == object) {
									//If the element is in the current region, then save its coordinates in temporary stacks
									iTmpStack.push_back(iStack[iNew]);
									jTmpStack.push_back(jStack[jNew]);
									kTmpStack.push_back(kStack[kNew] - 1);
								}
								else {
									//Not in region, update status and go to the next neighbor
									statusArray[kStack[kNew] - 1][x_new(iStack[iNew], jStack[jNew], xDim)] = true;
								}
							}
							//down neighbor
							if (kStack[kNew] < zDim - 1 && statusArray[kStack[kNew] + 1][x_new(iStack[iNew], jStack[jNew], xDim)] == false) {
								if (imageDataPtr[kStack[kNew] + 1][x_new(iStack[iNew], jStack[jNew], xDim)] == object) {
									//If the element is in the current region, then save its coordinates in temporary stacks
									iTmpStack.push_back(iStack[iNew]);
									jTmpStack.push_back(jStack[jNew]);
									kTmpStack.push_back(kStack[kNew] + 1);
								}
								else {
									//Not in region, update status and go to the next neighbor
									statusArray[kStack[kNew] + 1][x_new(iStack[iNew], jStack[jNew], xDim)] = true;
								}
							}
							//right neighbor
							if (iStack[iNew] > 0 && statusArray[kStack[kNew]][x_new(iStack[iNew] - 1, jStack[jNew], xDim)] == false) {
								if (imageDataPtr[kStack[kNew]][x_new(iStack[iNew] - 1, jStack[jNew], xDim)] == object) {
									//If the element is in the current region, then save its coordinates in temporary stacks
									iTmpStack.push_back(iStack[iNew] - 1);
									jTmpStack.push_back(jStack[jNew]);
									kTmpStack.push_back(kStack[kNew]);
								}
								else {
									//Not in region, update status and go to the next neighbor
									statusArray[kStack[kNew]][x_new(iStack[iNew] - 1, jStack[jNew], xDim)] = true;
								}
							}
							//left neighbor
							if (iStack[iNew] < xDim - 1 && statusArray[kStack[kNew]][x_new(iStack[iNew] + 1, jStack[jNew], xDim)] == false) {
								if (imageDataPtr[kStack[kNew]][x_new(iStack[iNew] + 1, jStack[jNew], xDim)] == object) {
									//If the element is in the current region, then save its coordinates in temporary stacks
									iTmpStack.push_back(iStack[iNew] + 1);
									jTmpStack.push_back(jStack[jNew]);
									kTmpStack.push_back(kStack[kNew]);
								}
								else {
									//Not in region, update status and go to the next neighbor
									statusArray[kStack[kNew]][x_new(iStack[iNew] + 1, jStack[jNew], xDim)] = true;
								}
							}
							//front neighbor
							if (jStack[jNew] > 0 && statusArray[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] - 1, xDim)] == false) {
								if (imageDataPtr[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] - 1, xDim)] == object) {
									//If the element is in the current region, then save its coordinates in temporary stacks
									iTmpStack.push_back(iStack[iNew]);
									jTmpStack.push_back(jStack[jNew] - 1);
									kTmpStack.push_back(kStack[kNew]);
								}
								else {
									//Not in region, update status and go to the next neighbor
									statusArray[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] - 1, xDim)] = true;
								}
							}
							//behind neighbor
							if (jStack[jNew] < yDim - 1 && statusArray[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] + 1, xDim)] == false) {
								if (imageDataPtr[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] + 1, xDim)] == object) {
									//If the element is in the current region, then save its coordinates in temporary stacks
									iTmpStack.push_back(iStack[iNew]);
									jTmpStack.push_back(jStack[jNew] + 1);
									kTmpStack.push_back(kStack[kNew]);
								}
								else {
									//Not in region, update status and go to the next neighbor
									statusArray[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] + 1, xDim)] = true;
								}
							}
							//updating of processed element befor removal
							segmentedImage[kStack[kNew]][x_new(iStack[iNew], jStack[jNew], xDim)] = label;
							statusArray[kStack[kNew]][x_new(iStack[iNew], jStack[jNew], xDim)] = true;
							//Remove the processed element of the initial stacks
							iStack.pop_back();
							jStack.pop_back();
							kStack.pop_back();
							//Add new found neighbors in initial stacks
							//we can do it once because they have the same size
							for (n = 0; n < iTmpStack.size(); n++) {
								//if any neighbors have been found iTmpStack.size() = 0
								//and nothing will happen
								iStack.push_back(iTmpStack[n]);
							}
							for (n = 0; n < jTmpStack.size(); n++) {
								//if any neighbors have been found jTmpStack.size() = 0
								//and nothing will happen
								jStack.push_back(jTmpStack[n]);
							}
							for (n = 0; n < kTmpStack.size(); n++) {
								//if any neighbors have been found jTmpStack.size() = 0
								//and nothing will happen
								kStack.push_back(kTmpStack[n]);
							}
							//empty the temporary stacks
							//we can do it once because they have the same size
							while (iTmpStack.size() > 0) {
								iTmpStack.pop_back();
							}
							while (jTmpStack.size() > 0) {
								jTmpStack.pop_back();
							}
							while (kTmpStack.size() > 0) {
								kTmpStack.pop_back();
							}
						}
						label++;
					}
					else {
						statusArray[k][x_new(i, j, xDim)] = true;
					}
				}

			}
		}
	}
	return true;
}

//Erosion ---> Shrink object and remove boundaries pixels
bool erosion3D(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background) {
	if (imageDataPtr == NULL)
		return false;

	size_t i, j, k, cpt;

	for (k = 0; k < zDim; k++) {
		for (i = 0; i < xDim; i++) {
			for (j = 0; j < yDim; j++) {

				cpt = 0;

				if (k == zDim - 1) {
					if (i == xDim - 1) {
						if (j == yDim - 1) {
							continue;
						}
						else {
							if (imageDataPtr[k][x_new(i, j + 1, xDim)] == background) {
								imageDataPtr[k][x_new(i, j, xDim)] = background;
							}
						}
					}
					else {
						if (j == yDim - 1) {
							if (imageDataPtr[k][x_new(i + 1, j, xDim)] == background) {
								imageDataPtr[k][x_new(i, j, xDim)] = background;
							}
						}
						else {
							cpt = 0;

							if (imageDataPtr[k][x_new(i + 1, j, xDim)] == object) {
								cpt++;
							}
							if (imageDataPtr[k][x_new(i, j + 1, xDim)] == object) {
								cpt++;
							}
							if (imageDataPtr[k][x_new(i + 1, j + 1, xDim)] == object) {
								cpt++;
							}
							if (cpt == 3) {
								imageDataPtr[k][x_new(i, j, xDim)] = object;
							}
							else {
								imageDataPtr[k][x_new(i, j, xDim)] = background;
							}
						}
					}
				}
				else {
					if (i == xDim - 1) {
						if (j == yDim - 1) {
							if (imageDataPtr[k + 1][x_new(i, j, xDim)] == background) {
								imageDataPtr[k][x_new(i, j, xDim)] = background;
							}
						}
						else {
							cpt = 0;

							if (imageDataPtr[k + 1][x_new(i, j, xDim)] == object) {
								cpt++;
							}
							if (imageDataPtr[k][x_new(i, j + 1, xDim)] == object) {
								cpt++;
							}
							if (imageDataPtr[k + 1][x_new(i, j + 1, xDim)] == object) {
								cpt++;
							}
							if (cpt == 3) {
								imageDataPtr[k][x_new(i, j, xDim)] = object;
							}
							else {
								imageDataPtr[k][x_new(i, j, xDim)] = background;
							}
						}
					}
					else {
						if (j == yDim - 1) {
							cpt = 0;

							if (imageDataPtr[k][x_new(i + 1, j, xDim)] == object) {
								cpt++;
							}
							if (imageDataPtr[k + 1][x_new(i, j, xDim)] == object) {
								cpt++;
							}
							if (imageDataPtr[k + 1][x_new(i + 1, j, xDim)] == object) {
								cpt++;
							}
							if (cpt == 3) {
								imageDataPtr[k][x_new(i, j, xDim)] = object;
							}
							else {
								imageDataPtr[k][x_new(i, j, xDim)] = background;
							}
						}
						else {
							cpt = 0;

							if (imageDataPtr[k][x_new(i, j + 1, xDim)] == object) {
								cpt++;
							}
							if (imageDataPtr[k][x_new(i + 1, j, xDim)] == object) {
								cpt++;
							}
							if (imageDataPtr[k][x_new(i + 1, j + 1, xDim)] == object) {
								cpt++;
							}
							if (imageDataPtr[k + 1][x_new(i, j, xDim)] == object) {
								cpt++;
							}
							if (imageDataPtr[k + 1][x_new(i + 1, j, xDim)] == object) {
								cpt++;
							}
							if (imageDataPtr[k + 1][x_new(i, j + 1, xDim)] == object) {
								cpt++;
							}
							if (imageDataPtr[k + 1][x_new(i + 1, j + 1, xDim)] == object) {
								cpt++;
							}
							if (cpt == 7) {
								imageDataPtr[k][x_new(i, j, xDim)] = object;
							}
							else {
								imageDataPtr[k][x_new(i, j, xDim)] = background;
							}
						}
					}
				}
			}
		}
	}

	return true;
}