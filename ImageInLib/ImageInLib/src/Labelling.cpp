//#include <iostream>
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

//for 2D images
bool regionLabelling(int** imageDataPtr, int** segmentedImage, int xDim, int yDim, int background, int object)
{
	if (imageDataPtr == NULL) {
		return false;
	}
	if (segmentedImage == NULL) {
		return false;
	}

	int regionCounter = 1, min = 0, max = 0, i, j, n, m;

	for (i = 0; i < xDim; i++) {
		for (j = 0; j < yDim; j++) {

			if (i == 0) {
				if (j == 0) {
					if (imageDataPtr[i][j] == object) {
						segmentedImage[i][j] = regionCounter;
						regionCounter++;
					}
				}
				else {
					if (imageDataPtr[i][j] == object) {
						if (imageDataPtr[i][j - 1] == background) {
							segmentedImage[i][j] = regionCounter;
							regionCounter++;
						}
						else {
							segmentedImage[i][j] = segmentedImage[i][j - 1];
						}
					}
				}
			}
			else {
				if (j == 0) {
					if (imageDataPtr[i][j] == object) {
						if (imageDataPtr[i - 1][j] == background) {
							segmentedImage[i][j] = regionCounter;
							regionCounter++;
						}
						else {
							segmentedImage[i][j] = segmentedImage[i - 1][j];
						}
					}
				}
				else {

					if (imageDataPtr[i][j] == object) {
						if (imageDataPtr[i - 1][j] == background && imageDataPtr[i][j - 1] == background) {
							segmentedImage[i][j] = regionCounter;
							regionCounter++;
						}
						if (imageDataPtr[i - 1][j] == object && imageDataPtr[i][j - 1] == background) {
							segmentedImage[i][j] = segmentedImage[i - 1][j];
						}
						if (imageDataPtr[i - 1][j] == background && imageDataPtr[i][j - 1] == object) {
							segmentedImage[i][j] = segmentedImage[i][j - 1];
						}
						if (imageDataPtr[i - 1][j] == object && imageDataPtr[i][j - 1] == object) {

							segmentedImage[i][j] = segmentedImage[i - 1][j];

							if (segmentedImage[i][j - 1] != segmentedImage[i - 1][j]) {

								max = segmentedImage[i - 1][j];
								min = segmentedImage[i][j - 1];

								if (max < segmentedImage[i][j - 1]) {
									max = segmentedImage[i][j - 1];
									min = segmentedImage[i - 1][j];
								}

								for (n = i; n >= 0; n--) {
									for (m = j; m >= 0; m--) {
										if (segmentedImage[n][m] == max) {
											segmentedImage[n][m] = min;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	return true;
}

bool labelling(int* imageDataPtr, int* segmentedImage, bool* statusArray, int xDim, int yDim, int object) {

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

			xd = x_new(i, j, xDim);

			if (statusArray[xd] == false) {
				if (imageDataPtr[xd] == object) {
					//top neighbor
					if (i > 0 && statusArray[x_new(i - 1, j, xDim)] == false) {
						if (imageDataPtr[x_new(i - 1, j, xDim)] == object) {
							//if element is in region add its coordinates in stacks
							iStack.push_back(i - 1);
							jStack.push_back(j);
						}
						else {
							//it's not object, so, update its status
							statusArray[x_new(i - 1, j, xDim)] = true;
						}
					}
					//right neighbor
					if (j < yDim - 1 && statusArray[x_new(i, j + 1, xDim)] == false) {
						if (imageDataPtr[x_new(i, j + 1, xDim)] == object) {
							//if element is in region add its coordinates in stacks
							iStack.push_back(i);
							jStack.push_back(j + 1);
						}
						else {
							statusArray[x_new(i, j + 1, xDim)] = true;
						}
					}
					//bottom neighbor
					if (i < xDim - 1 && statusArray[x_new(i + 1, j, xDim)] == false) {
						if (imageDataPtr[x_new(i + 1, j, xDim)] == object) {
							//if element is in region add its coodinates in stacks
							iStack.push_back(i + 1);
							jStack.push_back(j);
						}
						else {
							statusArray[x_new(i + 1, j, xDim)] = true;
						}
					}
					//left neighbor
					if (j > 0 && statusArray[x_new(i, j - 1, xDim)] == false) {
						if (imageDataPtr[x_new(i, j - 1, xDim)] == object) {
							//if element is in region add its coodinates in stacks
							iStack.push_back(i);
							jStack.push_back(j - 1);
						}
						else {
							statusArray[x_new(i, j - 1, xDim)] = true;
						}
					}
					//after checking all neighbors of current element
					//I give its label, and update its status
					segmentedImage[xd] = label;
					statusArray[xd] = true;
					//Now I check neighbors of its neighbors saved in the stacks
					//I start by the last added element in stacks
					//If there is no neighbor iStack.size() = jStack.size() = 0, and the while loop will no be ran
					while (iStack.size() > 0 && jStack.size() > 0) { //One is enought because they have same size
						iNew = iStack.size() - 1;
						jNew = jStack.size() - 1;
						//top neighbor
						if (iStack[iNew] > 0 && statusArray[x_new(iStack[iNew] - 1, jStack[jNew], xDim)] == false) {
							if (imageDataPtr[x_new(iStack[iNew] - 1, jStack[jNew], xDim)] == object) {
								//If the element is in the current region, then save its coordinates in temporary stacks
								iTmpStack.push_back(iStack[iNew] - 1);
								jTmpStack.push_back(jStack[jNew]);
							}
							else {
								//Not in region, update status and go to the next neighbor
								statusArray[x_new(iStack[iNew] - 1, jStack[jNew], xDim)] = true;
							}
						}
						//right neighbor
						if (jStack[jNew] < yDim - 1 && statusArray[x_new(iStack[iNew], jStack[jNew] + 1, xDim)] == false) {
							if (imageDataPtr[x_new(iStack[iNew], jStack[jNew] + 1, xDim)] == object) {
								//If the element is in the current region, then save its coordinates in temporary stacks
								iTmpStack.push_back(iStack[iNew]);
								jTmpStack.push_back(jStack[jNew] + 1);
							}
							else {
								statusArray[x_new(iStack[iNew], jStack[jNew] + 1, xDim)] = true;
							}
						}
						//bottom neighbor
						if (iStack[iNew] < xDim + 1 && statusArray[x_new(iStack[iNew] + 1, jStack[jNew], xDim)] == false) {
							if (imageDataPtr[x_new(iStack[iNew] + 1, jStack[jNew], xDim)] == object) {
								//If the element is in the current region, then save its coordinates in temporary stacks
								iTmpStack.push_back(iStack[iNew] + 1);
								jTmpStack.push_back(jStack[jNew]);
							}
							else {
								statusArray[x_new(iStack[iNew] + 1, jStack[jNew], xDim)] = true;
							}
						}
						//left neighbor
						if (jStack[jNew] > 0 && statusArray[x_new(iStack[iNew], jStack[jNew] - 1, xDim)] == false) {
							if (imageDataPtr[x_new(iStack[iNew], jStack[jNew] - 1, xDim)] == object) {
								//If the element is in the current region, then save its coordinates in temporary stacks
								iTmpStack.push_back(iStack[iNew]);
								jTmpStack.push_back(jStack[jNew] - 1);
							}
							else {
								statusArray[x_new(iStack[iNew], jStack[jNew] - 1, xDim)] = true;
							}
						}
						//updating of processed element befor removal
						segmentedImage[x_new(iStack[iNew], jStack[jNew], xDim)] = label;
						statusArray[x_new(iStack[iNew], jStack[jNew], xDim)] = true;
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
					statusArray[xd] = true;
				}
			}
		}
	}
	return true;
}

bool initialization2dArray(int** imageDataPtr, int xDim, int yDim, int value)
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

bool minorRegionRemoval(int** imageDataPtr, int** segmentedImage, int* CountingArray, int xDim, int yDim, int size)
{
	if (imageDataPtr == NULL) {
		return false;
	}
	if (segmentedImage == NULL) {
		return false;
	}
	if (CountingArray == NULL) {
		return false;
	}

	for (int i = 0; i < xDim; i++) {
		for (int j = 0; j < yDim; j++) {
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

//=================
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

//=================


//for 3D images

/*
bool fixEquivalence(dataType** segmentedImage, size_t xDim, size_t x, size_t y, size_t z, dataType minV, dataType maxV, bool parallize, size_t nbtreads) {

	int l, m, n;

	if (parallize == true) {
		omp_set_dynamic(0);
		omp_set_num_threads(nbtreads);

		if (z == 0) {
#pragma omp parallel
			{
				//do it just for x and y
#pragma omp for private(m,n) schedule(static) nowait
				for (m = x; m >= 0; m--) {
					for (n = y; n >= 0; n--) {
						if (segmentedImage[z][x_new(m, n, xDim)] == maxV) {
							segmentedImage[z][x_new(m, n, xDim)] = minV;
						}
					}
				}
			}
		}
		else {
			if (x == 0) {
				//do it for k and y
#pragma omp parallel
				{
#pragma omp for private(l,n) schedule(static) nowait
					for (l = z; l >= 0; l--) {
						for (n = y; n >= 0; n--) {
							if (segmentedImage[l][x_new(x, n, xDim)] == maxV) {
								segmentedImage[l][x_new(x, n, xDim)] = minV;
							}
						}
					}
				}
			}
			else {
				if (y == 0) {
					//do it for k and x
#pragma omp parallel
					{
#pragma omp for private(l,m) schedule(static) nowait
						for (l = z; l >= 0; l--) {
							for (m = x; m >= 0; m--) {
								if (segmentedImage[l][x_new(m, y, xDim)] == maxV) {
									segmentedImage[l][x_new(m, y, xDim)] = minV;
								}
							}
						}
					}
				}
				else {
					// do it for k,x,y
#pragma omp parallel
					{
#pragma omp for private(l,m,n) schedule(static) nowait
						for (l = z; l >= 0; l--) {
							for (m = x; m >= 0; m--) {
								for (n = y; n >= 0; n--) {
									if (segmentedImage[l][x_new(m, n, xDim)] == maxV) {
										segmentedImage[l][x_new(m, n, xDim)] = minV;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	else {
		if (z == 0) {
			//do it just for x and y
			for (m = x; m >= 0; m--) {
				for (n = y; n >= 0; n--) {
					if (segmentedImage[z][x_new(m, n, xDim)] == maxV) {
						segmentedImage[z][x_new(m, n, xDim)] = minV;
					}
				}
			}
		}
		else {
			if (x == 0) {
				//do it for k and y
				for (l = z; l >= 0; l--) {
					for (n = y; n >= 0; n--) {
						if (segmentedImage[l][x_new(x, n, xDim)] == maxV) {
							segmentedImage[l][x_new(x, n, xDim)] = minV;
						}
					}
				}
			}
			else {
				if (y == 0) {
					//do it for k and x
					for (l = z; l >= 0; l--) {
						for (m = x; m >= 0; m--) {
							if (segmentedImage[l][x_new(m, y, xDim)] == maxV) {
								segmentedImage[l][x_new(m, y, xDim)] = minV;
							}
						}
					}
				}
				else {
					// do it for k,x,y
					for (l = z; l >= 0; l--) {
						for (m = x; m >= 0; m--) {
							for (n = y; n >= 0; n--) {
								if (segmentedImage[l][x_new(m, n, xDim)] == maxV) {
									segmentedImage[l][x_new(m, n, xDim)] = minV;
								}
							}
						}
					}
				}
			}
		}

	}

	return true;
}
*/


bool fixEquivalence(dataType** segmentedImage, const size_t xDim, size_t x, size_t y, size_t z, dataType minV, dataType maxV, bool parallize, size_t nbtreads) {

	int l, m, n;

	if (parallize == true) {
		omp_set_dynamic(0);
		omp_set_num_threads(nbtreads);
#pragma omp parallel
		{
#pragma omp parallel for schedule(static) private(l,m,n)
			for (l = z; l >= 0; l--) {
				for (m = x; m >= 0; m--) {
					for (n = y; n >= 0; n--) {
						if (z == 0) {
							if (segmentedImage[z][x_new(m, n, xDim)] == maxV) {
								segmentedImage[z][x_new(m, n, xDim)] = minV;
							}
						}
						else {
							if (x == 0) {
								if (segmentedImage[l][x_new(x, n, xDim)] == maxV) {
									segmentedImage[l][x_new(x, n, xDim)] = minV;
								}
							}
							else {
								if (y == 0) {
									if (segmentedImage[l][x_new(m, y, xDim)] == maxV) {
										segmentedImage[l][x_new(m, y, xDim)] = minV;
									}
								}
								else {
									if (segmentedImage[l][x_new(m, n, xDim)] == maxV) {
										segmentedImage[l][x_new(m, n, xDim)] = minV;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	else {
		for (l = z; l >= 0; l--) {
			for (m = x; m >= 0; m--) {
				for (n = y; n >= 0; n--) {
					if (z == 0) {
						if (segmentedImage[z][x_new(m, n, xDim)] == maxV) {
							segmentedImage[z][x_new(m, n, xDim)] = minV;
						}
					}
					else {
						if (x == 0) {
							if (segmentedImage[l][x_new(x, n, xDim)] == maxV) {
								segmentedImage[l][x_new(x, n, xDim)] = minV;
							}
						}
						else {
							if (y == 0) {
								if (segmentedImage[l][x_new(m, y, xDim)] == maxV) {
									segmentedImage[l][x_new(m, y, xDim)] = minV;
								}
							}
							else {
								if (segmentedImage[l][x_new(m, n, xDim)] == maxV) {
									segmentedImage[l][x_new(m, n, xDim)] = minV;
								}
							}
						}
					}
				}
			}
		}
	}
	return true;
}

bool regionLabelling3D(dataType** imageDataPtr, dataType** segmentedImage, const size_t xDim, const size_t yDim, const size_t zDim, dataType background, dataType object, bool parallize, size_t nbtreads)
{

	if (imageDataPtr == NULL) {
		return false;
	}
	if (segmentedImage == NULL) {
		return false;
	}

	size_t k, i, j, xd;
	int l, m, n, p;
	dataType regionCounting = 1, minV = 0, maxV = 0;

	for (k = 0; k < zDim; k++) {
		for (i = 0; i < xDim; i++) {
			for (j = 0; j < yDim; j++) {

				xd = x_new(i, j, xDim);

				if (k == 0) {
					if (i == 0) {
						if (j == 0) {
							if (imageDataPtr[k][xd] == object) {
								segmentedImage[k][xd] = regionCounting;
								regionCounting++;
							}
						}
						else {
							if (imageDataPtr[k][xd] == object) {
								if (imageDataPtr[k][x_new(i, j - 1, xDim)] == background) {
									segmentedImage[k][xd] = regionCounting;
									regionCounting++;
								}
								else {
									segmentedImage[k][xd] = segmentedImage[k][x_new(i, j - 1, xDim)];
								}
							}
						}
					}
					else {
						if (j == 0) {
							if (imageDataPtr[k][xd] == object) {
								if (imageDataPtr[k][x_new(i - 1, j, xDim)] == background) {
									segmentedImage[k][xd] = regionCounting;
									regionCounting++;
								}
								else {
									segmentedImage[k][xd] = segmentedImage[k][x_new(i - 1, j, xDim)];
								}
							}
						}
						else {
							if (imageDataPtr[k][xd] == object) {
								if (imageDataPtr[k][x_new(i - 1, j, xDim)] == background && imageDataPtr[k][x_new(i, j - 1, xDim) == background]) {
									segmentedImage[k][xd] = regionCounting;
									regionCounting++;
								}
								if (imageDataPtr[k][x_new(i - 1, j, xDim)] == object && imageDataPtr[k][x_new(i, j - 1, xDim)] == background) {
									segmentedImage[k][xd] = segmentedImage[k][x_new(i - 1, j, xDim)];
								}
								if (imageDataPtr[k][x_new(i - 1, j, xDim)] == background && imageDataPtr[k][x_new(i, j - 1, xDim)] == object) {
									segmentedImage[k][xd] = segmentedImage[k][x_new(i, j - 1, xDim)];
								}
								if (imageDataPtr[k][x_new(i - 1, j, xDim)] == object && imageDataPtr[k][x_new(i, j - 1, xDim)] == object) {

									segmentedImage[k][xd] = segmentedImage[k][x_new(i, j - 1, xDim)];

									if (segmentedImage[k][x_new(i, j - 1, xDim)] != segmentedImage[k][x_new(i - 1, j, xDim)]) {
										/*
										max = segmentedImage[k][x_new(i - 1, j, xDim)];
										min = segmentedImage[k][x_new(i, j - 1, xDim)];

										if (max < segmentedImage[k][x_new(i, j - 1, xDim)]) {
											max = segmentedImage[k][x_new(i, j - 1, xDim)];
											min = segmentedImage[k][x_new(i - 1, j, xDim)];
										}

										for (m = i; m >= 0; m--) {
											for (n = j; n >= 0; n--) {
												if (segmentedImage[k][x_new(m, n, xDim)] == max) {
													segmentedImage[k][x_new(m, n, xDim)] = min;
												}
											}
										}
										*/
										minV = min(segmentedImage[k][x_new(i - 1, j, xDim)], segmentedImage[k][x_new(i, j - 1, xDim)]);
										maxV = max(segmentedImage[k][x_new(i - 1, j, xDim)], segmentedImage[k][x_new(i, j - 1, xDim)]);
										fixEquivalence(segmentedImage, xDim, i, j, k, minV, maxV, parallize, nbtreads);
									}
								}
							}
						}
					}
				}
				else {
					if (i == 0) {
						if (j == 0)
						{
							if (imageDataPtr[k][xd] == object) {
								if (imageDataPtr[k - 1][xd] == background) {
									segmentedImage[k][xd] = regionCounting;
									regionCounting++;
								}
								else {
									segmentedImage[k][xd] = segmentedImage[k - 1][xd];
								}
							}
						}
						else {
							if (imageDataPtr[k][xd] == object) {
								if (imageDataPtr[k - 1][xd] == background && imageDataPtr[k][x_new(i, j - 1, xDim)] == background) {
									segmentedImage[k][xd] = regionCounting;
									regionCounting++;
								}
								if (imageDataPtr[k - 1][xd] == object && imageDataPtr[k][x_new(i, j - 1, xDim)] == background) {
									segmentedImage[k][xd] = segmentedImage[k - 1][xd];
								}
								if (imageDataPtr[k - 1][xd] == background && imageDataPtr[k][x_new(i, j - 1, xDim)] == object) {
									segmentedImage[k][xd] = segmentedImage[k][x_new(i, j - 1, xDim)];
								}
								if (imageDataPtr[k - 1][xd] == object && imageDataPtr[k][x_new(i, j - 1, xDim)] == object) {

									segmentedImage[k][xd] = segmentedImage[k - 1][xd];

									if (segmentedImage[k - 1][xd] != segmentedImage[k][x_new(i, j - 1, xDim)]) {
										/*
										max = segmentedImage[k - 1][xd];
										min = segmentedImage[k][x_new(i, j - 1, xDim)];

										if (max < segmentedImage[k][x_new(i, j - 1, xDim)]) {
											max = segmentedImage[k][x_new(i, j - 1, xDim)];
											min = segmentedImage[k - 1][xd];
										}

										for (l = k; l >= 0; l--) {
											for (n = j; n >= 0; n--) {
												if (segmentedImage[l][x_new(i, n, xDim)] == max) {
													segmentedImage[l][x_new(i, n, xDim)] = min;
												}
											}
										}
										*/
										minV = min(segmentedImage[k - 1][xd], segmentedImage[k][x_new(i, j - 1, xDim)]);
										maxV = max(segmentedImage[k - 1][xd], segmentedImage[k][x_new(i, j - 1, xDim)]);
										fixEquivalence(segmentedImage, xDim, i, j, k, minV, maxV, parallize, nbtreads);
									}
								}
							}
						}
					}
					else {
						if (j == 0) {
							if (imageDataPtr[k][xd] == object) {
								if (imageDataPtr[k - 1][xd] == background && imageDataPtr[k][x_new(i - 1, j, xDim)] == background) {
									segmentedImage[k][xd] = regionCounting;
									regionCounting++;
								}
								if (imageDataPtr[k - 1][xd] == object && imageDataPtr[k][x_new(i - 1, j, xDim)] == background) {
									segmentedImage[k][xd] = segmentedImage[k - 1][xd];
								}
								if (imageDataPtr[k - 1][xd] == background && imageDataPtr[k][x_new(i - 1, j, xDim)] == object) {
									segmentedImage[k][xd] = segmentedImage[k][x_new(i - 1, j, xDim)];
								}
								if (imageDataPtr[k - 1][xd] == object && imageDataPtr[k][x_new(i - 1, j, xDim)] == object) {

									segmentedImage[k][xd] = segmentedImage[k - 1][xd];

									if (segmentedImage[k - 1][xd] != segmentedImage[k][x_new(i - 1, j, xDim)]) {
										/*
										max = segmentedImage[k - 1][xd];
										min = segmentedImage[k][x_new(i - 1, j, xDim)];

										if (max < segmentedImage[k][x_new(i - 1, j, xDim)]) {
											max = segmentedImage[k][x_new(i - 1, j, xDim)];
											min = segmentedImage[k - 1][xd];
										}

										for (l = k; l >= 0; l--) {
											for (m = i; m >= 0; m--) {
												if (segmentedImage[l][x_new(m, j, xDim)] == max) {
													segmentedImage[l][x_new(m, j, xDim)] = min;
												}
											}
										}
										*/
										minV = min(segmentedImage[k - 1][xd], segmentedImage[k][x_new(i - 1, j, xDim)]);
										maxV = max(segmentedImage[k - 1][xd], segmentedImage[k][x_new(i - 1, j, xDim)]);
										fixEquivalence(segmentedImage, xDim, i, j, k, minV, maxV, parallize, nbtreads);
									}
								}
							}
						}
						else {
							if (imageDataPtr[k][xd] == object) {
								if (imageDataPtr[k - 1][xd] == background && imageDataPtr[k][x_new(i - 1, j, xDim)] == background && imageDataPtr[k][x_new(i, j - 1, xDim) == background]) {
									segmentedImage[k][xd] = regionCounting;
									regionCounting++;
								}
								if (imageDataPtr[k - 1][xd] == object && imageDataPtr[k][x_new(i - 1, j, xDim)] == background && imageDataPtr[k][x_new(i, j - 1, xDim)] == background) {
									segmentedImage[k][xd] = segmentedImage[k - 1][xd];
								}
								if (imageDataPtr[k - 1][xd] == object && imageDataPtr[k][x_new(i - 1, j, xDim)] == object && imageDataPtr[k][x_new(i, j - 1, xDim)] == background) {

									segmentedImage[k][xd] = segmentedImage[k - 1][xd];

									if (segmentedImage[k - 1][xd] != segmentedImage[k][x_new(i - 1, j, xDim)]) {
										/*
										max = segmentedImage[k - 1][xd];
										min = segmentedImage[k][x_new(i - 1, j, xDim)];

										if (max < segmentedImage[k][x_new(i - 1, j, xDim)]) {
											max = segmentedImage[k][x_new(i - 1, j, xDim)];
											min = segmentedImage[k - 1][xd];
										}

										for (l = k; l >= 0; l--) {
											for (m = i; m >= 0; m--) {
												for (n = j; n >= 0; n--) {
													if (segmentedImage[l][x_new(m, n, xDim)] == max) {
														segmentedImage[l][x_new(m, n, xDim)] = min;
													}
												}
											}
										}
										*/
										minV = min(segmentedImage[k - 1][xd], segmentedImage[k][x_new(i - 1, j, xDim)]);
										maxV = max(segmentedImage[k - 1][xd], segmentedImage[k][x_new(i - 1, j, xDim)]);
										fixEquivalence(segmentedImage, xDim, i, j, k, minV, maxV, parallize, nbtreads);
									}
								}
								if (imageDataPtr[k - 1][xd] == object && imageDataPtr[k][x_new(i - 1, j, xDim)] == background && imageDataPtr[k][x_new(i, j - 1, xDim)] == object) {

									segmentedImage[k][xd] = segmentedImage[k - 1][xd];

									if (segmentedImage[k - 1][xd] != segmentedImage[k][x_new(i, j - 1, xDim)]) {
										/*
										max = segmentedImage[k - 1][xd];
										min = segmentedImage[k][x_new(i, j - 1, xDim)];

										if (max < segmentedImage[k][x_new(i, j - 1, xDim)]) {
											max = segmentedImage[k][x_new(i, j - 1, xDim)];
											min = segmentedImage[k - 1][xd];
										}

										for (l = k; l >= 0; l--) {
											for (m = i; m >= 0; m--) {
												for (n = j; n >= 0; n--) {
													if (segmentedImage[l][x_new(m, n, xDim)] == max) {
														segmentedImage[l][x_new(m, n, xDim)] = min;
													}
												}
											}
										}
										*/
										minV = min(segmentedImage[k - 1][xd], segmentedImage[k][x_new(i, j - 1, xDim)]);
										maxV = max(segmentedImage[k - 1][xd], segmentedImage[k][x_new(i, j - 1, xDim)]);
										fixEquivalence(segmentedImage, xDim, i, j, k, minV, maxV, parallize, nbtreads);
									}
								}
								if (imageDataPtr[k - 1][xd] == background && imageDataPtr[k][x_new(i - 1, j, xDim)] == object && imageDataPtr[k][x_new(i, j - 1, xDim)] == background) {
									segmentedImage[k][xd] = segmentedImage[k][x_new(i - 1, j, xDim)];
								}
								if (imageDataPtr[k - 1][xd] == background && imageDataPtr[k][x_new(i - 1, j, xDim)] == background && imageDataPtr[k][x_new(i, j - 1, xDim)] == object) {
									segmentedImage[k][xd] = segmentedImage[k][x_new(i, j - 1, xDim)];
								}
								if (imageDataPtr[k - 1][xd] == background && imageDataPtr[k][x_new(i - 1, j, xDim)] == object && imageDataPtr[k][x_new(i, j - 1, xDim)] == object) {

									segmentedImage[k][xd] = segmentedImage[k][x_new(i - 1, j, xDim)];

									if (segmentedImage[k][x_new(i - 1, j, xDim)] != segmentedImage[k][x_new(i, j - 1, xDim)]) {
										/*
										max = segmentedImage[k][x_new(i - 1, j, xDim)];
										min = segmentedImage[k][x_new(i, j - 1, xDim)];

										if (max < segmentedImage[k][x_new(i, j - 1, xDim)]) {
											max = segmentedImage[k][x_new(i, j - 1, xDim)];
											min = segmentedImage[k][x_new(i - 1, j, xDim)];
										}

										for (l = k; l >= 0; l--) {
											for (m = i; m >= 0; m--) {
												for (n = j; n >= 0; n--) {
													if (segmentedImage[l][x_new(m, n, xDim)] == max) {
														segmentedImage[l][x_new(m, n, xDim)] = min;
													}
												}
											}
										}
										*/
										minV = min(segmentedImage[k][x_new(i - 1, j, xDim)], segmentedImage[k][x_new(i, j - 1, xDim)]);
										maxV = max(segmentedImage[k][x_new(i - 1, j, xDim)], segmentedImage[k][x_new(i, j - 1, xDim)]);
										fixEquivalence(segmentedImage, xDim, i, j, k, minV, maxV, parallize, nbtreads);
									}
								}
								if (imageDataPtr[k - 1][xd] == object && imageDataPtr[k][x_new(i - 1, j, xDim)] == object && imageDataPtr[k][x_new(i, j - 1, xDim)] == object) {

									segmentedImage[k][xd] = segmentedImage[k - 1][xd];

									if (segmentedImage[k - 1][xd] != segmentedImage[k][x_new(i - 1, j, xDim)] || segmentedImage[k - 1][xd] != segmentedImage[k][x_new(i, j - 1, xDim)]) {

										dataType tab[3] = { segmentedImage[k - 1][xd], segmentedImage[k][x_new(i - 1, j, xDim)], segmentedImage[k][x_new(i, j - 1, xDim)] };
										minV = tab[0], maxV = tab[0];
										for (p = 0; p < 3; p++) {
											if (tab[p] > maxV) {
												maxV = tab[p];
											}
											if (tab[p] < minV) {
												minV = tab[p];
											}
										}
										/*
										for (l = k; l >= 0; l--) {
											for (m = i; m >= 0; m--) {
												for (n = j; n >= 0; n--) {
													if (segmentedImage[l][x_new(m, n, xDim)] == max) {
														segmentedImage[l][x_new(m, n, xDim)] = min;
													}
												}
											}
										}
										*/
										fixEquivalence(segmentedImage, xDim, i, j, k, minV, maxV, parallize, nbtreads);
									}
								}
							}
						}
					}
				}
			}
		}
	}

	return true;
}

bool labelling3D(dataType* imageDataPtr, dataType* segmentedImage, bool* statusArray, int xDim, int yDim, int zDim, dataType object) {

	if (imageDataPtr == NULL)
		return false;
	if (segmentedImage == NULL)
		return false;
	if (statusArray == NULL)
		return false;

	vector<int> iStack, jStack, kStack, iTmpStack, jTmpStack, kTmpStack;
	int i = 0, j = 0, k = 0, iNew = 0, jNew = 0, kNew = 0, n = 0, xd = 0, label = 1;

	//statusArray has false everywhere in the beginning
	//false---> non-processed and true---> processed
	//segmentedImage has 0 everywhere
	for (k = 0; k < zDim; k++) {
		for (i = 0; i < xDim; i++) {
			for (j = 0; j < yDim; j++) {

				xd = x_flat(i, j, k, xDim, yDim);

				if (statusArray[xd] == false) {
					if (imageDataPtr[xd] == object) {
						//top neighbor
						if (k > 0 && statusArray[x_flat(i, j, k - 1, xDim, yDim)] == false) {
							if (imageDataPtr[x_flat(i, j, k - 1, xDim, yDim)] == object) {
								//if element is in region add its coordinates in stacks
								iStack.push_back(i);
								jStack.push_back(j);
								kStack.push_back(k - 1);
							}
							else {
								//it's not object, so, update its status
								statusArray[x_flat(i, j, k - 1, xDim, yDim)] = true;
							}
						}
						//down neighbor
						if (k < zDim - 1 && statusArray[x_flat(i, j, k + 1, xDim, yDim)] == false) {
							if (imageDataPtr[x_flat(i, j, k + 1, xDim, yDim)] == object) {
								//if element is in region add its coordinates in stacks
								iStack.push_back(i);
								jStack.push_back(j);
								kStack.push_back(k + 1);
							}
							else {
								//it's not object, so, update its status
								statusArray[x_flat(i, j, k + 1, xDim, yDim)] = true;
							}
						}
						//left neighbor
						if (i > 0 && statusArray[x_flat(i - 1, j, k, xDim, yDim)] == false) {
							if (imageDataPtr[x_flat(i - 1, j, k, xDim, yDim)] == object) {
								//if element is in region add its coordinates in stacks
								iStack.push_back(i - 1);
								jStack.push_back(j);
								kStack.push_back(k);
							}
							else {
								statusArray[x_flat(i - 1, j, k, xDim, yDim)] = true;
							}
						}
						//right neighbor
						if (i < xDim - 1 && statusArray[x_flat(i + 1, j, k, xDim, yDim)] == false) {
							if (imageDataPtr[x_flat(i + 1, j, k, xDim, yDim)] == object) {
								//if element is in region add its coordinates in stacks
								iStack.push_back(i + 1);
								jStack.push_back(j);
								kStack.push_back(k);
							}
							else {
								statusArray[x_flat(i + 1, j, k, xDim, yDim)] = true;
							}
						}
						//front neighbor
						if (j > 0 && statusArray[x_flat(i, j - 1, k, xDim, yDim)] == false) {
							if (imageDataPtr[x_flat(i, j - 1, k, xDim, yDim)] == object) {
								//if element is in region add its coodinates in stacks
								iStack.push_back(i);
								jStack.push_back(j - 1);
								kStack.push_back(k);
							}
							else {
								statusArray[x_flat(i, j - 1, k, xDim, yDim)] = true;
							}
						}
						//behind neighbor
						if (j < yDim - 1 && statusArray[x_flat(i, j + 1, k, xDim, yDim)] == false) {
							if (imageDataPtr[x_flat(i, j + 1, k, xDim, yDim)] == object) {
								//if element is in region add its coodinates in stacks
								iStack.push_back(i);
								jStack.push_back(j + 1);
								kStack.push_back(k);
							}
							else {
								statusArray[x_flat(i, j + 1, k, xDim, yDim)] = true;
							}
						}
						//after checking all neighbors of current element
						//I give its label, and update its status
						segmentedImage[xd] = label;
						statusArray[xd] = true;
						//Now I check neighbors of its neighbors saved in the stacks
						//I start by the last added element in stacks
						//If there is no neighbor iStack.size() = jStack.size() = kStack.size() = 0 , and the while loop will no be ran
						while (iStack.size() > 0 && jStack.size() > 0 && kStack.size() > 0) {
							//One is enought because they have same size
							iNew = iStack.size() - 1;
							jNew = jStack.size() - 1;
							kNew = kStack.size() - 1;
							//top neighbor
							if (kStack[kNew] > 0 && statusArray[x_flat(iStack[iNew], jStack[jNew], kStack[kNew] - 1, xDim, yDim)] == false) {
								if (imageDataPtr[x_flat(iStack[iNew], jStack[jNew], kStack[kNew] - 1, xDim, yDim)] == object) {
									//If the element is in the current region, then save its coordinates in temporary stacks
									iTmpStack.push_back(iStack[iNew]);
									jTmpStack.push_back(jStack[jNew]);
									kTmpStack.push_back(kStack[kNew] - 1);
								}
								else {
									//Not in region, update status and go to the next neighbor
									statusArray[x_flat(iStack[iNew], jStack[jNew], kStack[kNew] - 1, xDim, yDim)] = true;
								}
							}
							//down neighbor
							if (kStack[kNew] < zDim - 1 && statusArray[x_flat(iStack[iNew], jStack[jNew], kStack[kNew] + 1, xDim, yDim)] == false) {
								if (imageDataPtr[x_flat(iStack[iNew], jStack[jNew], kStack[kNew] + 1, xDim, yDim)] == object) {
									//If the element is in the current region, then save its coordinates in temporary stacks
									iTmpStack.push_back(iStack[iNew]);
									jTmpStack.push_back(jStack[jNew]);
									kTmpStack.push_back(kStack[kNew] + 1);
								}
								else {
									//Not in region, update status and go to the next neighbor
									statusArray[x_flat(iStack[iNew], jStack[jNew], kStack[kNew] + 1, xDim, yDim)] = true;
								}
							}
							//right neighbor
							if (iStack[iNew] > 0 && statusArray[x_flat(iStack[iNew] - 1, jStack[jNew], kStack[kNew], xDim, yDim)] == false) {
								if (imageDataPtr[x_flat(iStack[iNew] - 1, jStack[jNew], kStack[kNew], xDim, yDim)] == object) {
									//If the element is in the current region, then save its coordinates in temporary stacks
									iTmpStack.push_back(iStack[iNew] - 1);
									jTmpStack.push_back(jStack[jNew]);
									kTmpStack.push_back(kStack[kNew]);
								}
								else {
									//Not in region, update status and go to the next neighbor
									statusArray[x_flat(iStack[iNew] - 1, jStack[jNew], kStack[kNew], xDim, yDim)] = true;
								}
							}
							//left neighbor
							if (iStack[iNew] < xDim - 1 && statusArray[x_flat(iStack[iNew] + 1, jStack[jNew], kStack[kNew], xDim, yDim)] == false) {
								if (imageDataPtr[x_flat(iStack[iNew] + 1, jStack[jNew], kStack[kNew], xDim, yDim)] == object) {
									//If the element is in the current region, then save its coordinates in temporary stacks
									iTmpStack.push_back(iStack[iNew] + 1);
									jTmpStack.push_back(jStack[jNew]);
									kTmpStack.push_back(kStack[kNew]);
								}
								else {
									//Not in region, update status and go to the next neighbor
									statusArray[x_flat(iStack[iNew] + 1, jStack[jNew], kStack[kNew], xDim, yDim)] = true;
								}
							}
							//front neighbor
							if (jStack[jNew] > 0 && statusArray[x_flat(iStack[iNew], jStack[jNew] - 1, kStack[kNew], xDim, yDim)] == false) {
								if (imageDataPtr[x_flat(iStack[iNew], jStack[jNew] - 1, kStack[kNew], xDim, yDim)] == object) {
									//If the element is in the current region, then save its coordinates in temporary stacks
									iTmpStack.push_back(iStack[iNew]);
									jTmpStack.push_back(jStack[jNew] - 1);
									kTmpStack.push_back(kStack[kNew]);
								}
								else {
									//Not in region, update status and go to the next neighbor
									statusArray[x_flat(iStack[iNew], jStack[jNew] - 1, kStack[kNew], xDim, yDim)] = true;
								}
							}
							//behind neighbor
							if (jStack[jNew] < yDim - 1 && statusArray[x_flat(iStack[iNew], jStack[jNew] + 1, kStack[kNew], xDim, yDim)] == false) {
								if (imageDataPtr[x_flat(iStack[iNew], jStack[jNew] + 1, kStack[kNew], xDim, yDim)] == object) {
									//If the element is in the current region, then save its coordinates in temporary stacks
									iTmpStack.push_back(iStack[iNew]);
									jTmpStack.push_back(jStack[jNew] + 1);
									kTmpStack.push_back(kStack[kNew]);
								}
								else {
									//Not in region, update status and go to the next neighbor
									statusArray[x_flat(iStack[iNew], jStack[jNew] + 1, kStack[kNew], xDim, yDim)] = true;
								}
							}
							//updating of processed element befor removal
							segmentedImage[x_flat(iStack[iNew], jStack[jNew], kStack[kNew], xDim, yDim)] = label;
							statusArray[x_flat(iStack[iNew], jStack[jNew], kStack[kNew], xDim, yDim)] = true;
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
						statusArray[xd] = true;
					}
				}

			}
		}
	}
	return true;
}

