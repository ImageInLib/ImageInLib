#include <iostream>
#include <climits>
#include <crtdbg.h>
#include <corecrt_malloc.h>
#include <omp.h>
#include<vector>
#include <common_vtk.h>
#include "common_functions.h"
#include "Labelling.h"
#include "morphological_change.h"

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

//for 3D images
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
								//the element is in current region, then add its coordinates in stacks
								kStack.push_back(k - 1);
								iStack.push_back(i);
								jStack.push_back(j);
							}
							else {
								//it's not object, so, update its status
								statusArray[k - 1][x_new(i, j, xDim)] = true;
							}
						}
						//down neighbor
						if (k < zDim - 1 && statusArray[k + 1][x_new(i, j, xDim)] == false) {
							if (imageDataPtr[k + 1][x_new(i, j, xDim)] == object) {
								//the element is in current region, then add its coordinates in stacks
								kStack.push_back(k + 1);
								iStack.push_back(i);
								jStack.push_back(j);
							}
							else {
								//it's not object, so, update its status
								statusArray[k + 1][x_new(i, j, xDim)] = true;
							}
						}
						//left neighbor
						if (i > 0 && statusArray[k][x_new(i - 1, j, xDim)] == false) {
							if (imageDataPtr[k][x_new(i - 1, j, xDim)] == object) {
								//the element is in region, then add its coordinates in stacks
								kStack.push_back(k);
								iStack.push_back(i - 1);
								jStack.push_back(j);
							}
							else {
								statusArray[k][x_new(i - 1, j, xDim)] = true;
							}
						}
						//right neighbor
						if (i < xDim - 1 && statusArray[k][x_new(i + 1, j, xDim)] == false) {
							if (imageDataPtr[k][x_new(i + 1, j, xDim)] == object) {
								//if element is in region add its coordinates in stacks
								kStack.push_back(k);
								iStack.push_back(i + 1);
								jStack.push_back(j);
							}
							else {
								statusArray[k][x_new(i + 1, j, xDim)] = true;
							}
						}
						//front neighbor
						if (j > 0 && statusArray[k][x_new(i, j - 1, xDim)] == false) {
							if (imageDataPtr[k][x_new(i, j - 1, xDim)] == object) {
								//if element is in region add its coodinates in stacks
								kStack.push_back(k);
								iStack.push_back(i);
								jStack.push_back(j - 1);
							}
							else {
								statusArray[k][x_new(i, j - 1, xDim)] = true;
							}
						}
						//behind neighbor
						if (j < yDim - 1 && statusArray[k][x_new(i, j + 1, xDim)] == false) {
							if (imageDataPtr[k][x_new(i, j + 1, xDim)] == object) {
								//the element is in region, then add its coodinates in stacks
								kStack.push_back(k);
								iStack.push_back(i);
								jStack.push_back(j + 1);
							}
							else {
								statusArray[k][x_new(i, j + 1, xDim)] = true;
							}
						}
						
						//after checking all neighbors of current element
						//we give its label, and update its status
						segmentedImage[k][x_new(i, j, xDim)] = label;
						statusArray[k][x_new(i, j, xDim)] = true;

						//Now we check neighbors of its neighbors saved in the stacks
						//I start by the last added element in stacks
						//If there is no neighbor iStack.size() = jStack.size() = kStack.size() = 0 , and the while loop will no be ran
						while (iStack.size() > 0 && jStack.size() > 0 && kStack.size() > 0) {
							//One is enought because they have same size

							//We work with the last element in the initial stacks
							kNew = kStack.size() - 1;
							iNew = iStack.size() - 1;
							jNew = jStack.size() - 1;

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
									//the element is in the current region, then save its coordinates in temporary stacks
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
									//the element is in the current region, then save its coordinates in temporary stacks
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
									//the element is in the current region, then save its coordinates in temporary stacks
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
							kStack.pop_back();
							iStack.pop_back();
							jStack.pop_back();

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

							//End of big while loop
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

bool regionGrowing(dataType** imageDataPtr, dataType** segmentedImage, bool** statusArray, const size_t xDim, const size_t yDim, const size_t zDim, dataType thres_min, dataType thres_max, Point3D * seedPoint) {

	if (imageDataPtr == NULL || segmentedImage == NULL || statusArray == NULL)
		return false;

	vector<size_t> iStack, jStack, kStack, iTmpStack, jTmpStack, kTmpStack;
	size_t i, j, k, iNew = 0, jNew = 0, kNew = 0, n = 0, dim2D = xDim * yDim;

	i = (size_t)seedPoint->x; j = (size_t)seedPoint->y; k = (size_t)seedPoint->z;

	//Processed the seed point
	//top neighbor
	if (k > 0 && statusArray[k - 1][x_new(i, j, xDim)] == false) {
		if (imageDataPtr[k - 1][x_new(i, j, xDim)] >= thres_min && imageDataPtr[k - 1][x_new(i, j, xDim)] <= thres_max) {
			//the element is in region, then add its coordinates in stacks
			kStack.push_back(k - 1);
			iStack.push_back(i);
			jStack.push_back(j);
		}
		else {
			statusArray[k - 1][x_new(i, j, xDim)] = true;
			segmentedImage[k - 1][x_new(i, j, xDim)] = 0;
		}
	}
	//bottom neighbor
	if (k < zDim - 1 && statusArray[k + 1][x_new(i, j, xDim)] == false) {
		if (imageDataPtr[k + 1][x_new(i, j, xDim)] >= thres_min && imageDataPtr[k + 1][x_new(i, j, xDim)] <= thres_max) {
			//the element is in region, then add its coordinates in stacks
			kStack.push_back(k + 1);
			iStack.push_back(i);
			jStack.push_back(j);
		}
		else {
			statusArray[k + 1][x_new(i, j, xDim)] = true;
			segmentedImage[k + 1][x_new(i, j, xDim)] = 0;
		}
	}
	//left neighbor
	if (i > 0 && statusArray[k][x_new(i - 1, j, xDim)] == false) {
		if (imageDataPtr[k][x_new(i - 1, j, xDim)] >= thres_min && imageDataPtr[k][x_new(i - 1, j, xDim)] <= thres_max) {
			//the element is in region, then add its coordinates in stacks
			kStack.push_back(k);
			iStack.push_back(i - 1);
			jStack.push_back(j);
		}
		else {
			statusArray[k][x_new(i - 1, j, xDim)] = true;
			segmentedImage[k][x_new(i - 1, j, xDim)] = 0;
		}
	}
	//right neighbor
	if (i < xDim - 1 && statusArray[k][x_new(i + 1, j, xDim)] == false) {
		if (imageDataPtr[k][x_new(i + 1, j, xDim)] >= thres_min && imageDataPtr[k][x_new(i + 1, j, xDim)] <= thres_max) {
			//the element is in region, then add its coordinates in stacks
			kStack.push_back(k);
			iStack.push_back(i + 1);
			jStack.push_back(j);
		}
		else {
			statusArray[k][x_new(i + 1, j, xDim)] = true;
			segmentedImage[k][x_new(i + 1, j, xDim)] = 0;
		}
	}
	//front neighbor
	if (j > 0 && statusArray[k][x_new(i, j - 1, xDim)] == false) {
		if (imageDataPtr[k][x_new(i, j - 1, xDim)] >= thres_min && imageDataPtr[k][x_new(i, j - 1, xDim)] <= thres_max) {
			//the element is in region, then add its coordinates in stacks
			kStack.push_back(k);
			iStack.push_back(i);
			jStack.push_back(j - 1);
		}
		else {
			statusArray[k][x_new(i, j - 1, xDim)] = true;
			segmentedImage[k][x_new(i, j - 1, xDim)] = 0;
		}
	}
	//behind neighbor
	if (j < yDim - 1 && statusArray[k][x_new(i, j + 1, xDim)] == false) {
		if (imageDataPtr[k][x_new(i, j + 1, xDim)] >= thres_min && imageDataPtr[k][x_new(i, j + 1, xDim)] <= thres_max) {
			//the element is in region, then add its coordinates in stacks
			kStack.push_back(k);
			iStack.push_back(i);
			jStack.push_back(j + 1);
		}
		else {
			statusArray[k][x_new(i, j + 1, xDim)] = true;
			segmentedImage[k][x_new(i, j + 1, xDim)] = 0;
		}
	}
	//Update the seed point status
	statusArray[k][x_new(i, j, xDim)] = true;
	//set  the seed point's label
	segmentedImage[k][x_new(i, j, xDim)] = imageDataPtr[k][x_new(i, j, xDim)];

	//Processed the other points
	while (iStack.size() > 0 && jStack.size() > 0 && kStack.size() > 0) {
		//One is enought because they have same size

		//We work with the last element in the initial stacks
		kNew = kStack.size() - 1;
		iNew = iStack.size() - 1;
		jNew = jStack.size() - 1;

		//top neighbor
		if (kStack[kNew] > 0 && statusArray[kStack[kNew] - 1][x_new(iStack[iNew], jStack[jNew], xDim)] == false) {
			if (imageDataPtr[kStack[kNew] - 1][x_new(iStack[iNew], jStack[jNew], xDim)] >= thres_min && imageDataPtr[kStack[kNew] - 1][x_new(iStack[iNew], jStack[jNew], xDim)] <= thres_max) {
				//If the element is in the current region, then save its coordinates in temporary stacks
				iTmpStack.push_back(iStack[iNew]);
				jTmpStack.push_back(jStack[jNew]);
				kTmpStack.push_back(kStack[kNew] - 1);
			}
			else {
				//Not in region, update status and go to the next neighbor
				statusArray[kStack[kNew] - 1][x_new(iStack[iNew], jStack[jNew], xDim)] = true;
				segmentedImage[kStack[kNew] - 1][x_new(iStack[iNew], jStack[jNew], xDim)] = 0;
			}
		}
		//down neighbor
		if (kStack[kNew] < zDim - 1 && statusArray[kStack[kNew] + 1][x_new(iStack[iNew], jStack[jNew], xDim)] == false) {
			if (imageDataPtr[kStack[kNew] + 1][x_new(iStack[iNew], jStack[jNew], xDim)] >= thres_min && imageDataPtr[kStack[kNew] + 1][x_new(iStack[iNew], jStack[jNew], xDim)] <= thres_max) {
				//the element is in the current region, then save its coordinates in temporary stacks
				iTmpStack.push_back(iStack[iNew]);
				jTmpStack.push_back(jStack[jNew]);
				kTmpStack.push_back(kStack[kNew] + 1);
			}
			else {
				//Not in region, update status and go to the next neighbor
				statusArray[kStack[kNew] + 1][x_new(iStack[iNew], jStack[jNew], xDim)] = true;
				segmentedImage[kStack[kNew] + 1][x_new(iStack[iNew], jStack[jNew], xDim)] = 0;
			}
		}
		//right neighbor
		if (iStack[iNew] > 0 && statusArray[kStack[kNew]][x_new(iStack[iNew] - 1, jStack[jNew], xDim)] == false) {
			if (imageDataPtr[kStack[kNew]][x_new(iStack[iNew] - 1, jStack[jNew], xDim)] >= thres_min && imageDataPtr[kStack[kNew]][x_new(iStack[iNew] - 1, jStack[jNew], xDim)] <= thres_max) {
				//the element is in the current region, then save its coordinates in temporary stacks
				iTmpStack.push_back(iStack[iNew] - 1);
				jTmpStack.push_back(jStack[jNew]);
				kTmpStack.push_back(kStack[kNew]);
			}
			else {
				//Not in region, update status and go to the next neighbor
				statusArray[kStack[kNew]][x_new(iStack[iNew] - 1, jStack[jNew], xDim)] = true;
				segmentedImage[kStack[kNew]][x_new(iStack[iNew] - 1, jStack[jNew], xDim)] = 0;
			}
		}
		//left neighbor
		if (iStack[iNew] < xDim - 1 && statusArray[kStack[kNew]][x_new(iStack[iNew] + 1, jStack[jNew], xDim)] == false) {
			if (imageDataPtr[kStack[kNew]][x_new(iStack[iNew] + 1, jStack[jNew], xDim)] >= thres_min && imageDataPtr[kStack[kNew]][x_new(iStack[iNew] + 1, jStack[jNew], xDim)] <= thres_max) {
				//If the element is in the current region, then save its coordinates in temporary stacks
				iTmpStack.push_back(iStack[iNew] + 1);
				jTmpStack.push_back(jStack[jNew]);
				kTmpStack.push_back(kStack[kNew]);
			}
			else {
				//Not in region, update status and go to the next neighbor
				statusArray[kStack[kNew]][x_new(iStack[iNew] + 1, jStack[jNew], xDim)] = true;
				segmentedImage[kStack[kNew]][x_new(iStack[iNew] + 1, jStack[jNew], xDim)] = 0;
			}
		}
		//front neighbor
		if (jStack[jNew] > 0 && statusArray[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] - 1, xDim)] == false) {
			if (imageDataPtr[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] - 1, xDim)] >= thres_min && imageDataPtr[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] - 1, xDim)] <= thres_max) {
				//If the element is in the current region, then save its coordinates in temporary stacks
				iTmpStack.push_back(iStack[iNew]);
				jTmpStack.push_back(jStack[jNew] - 1);
				kTmpStack.push_back(kStack[kNew]);
			}
			else {
				//Not in region, update status and give it label
				statusArray[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] - 1, xDim)] = true;
				segmentedImage[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] - 1, xDim)] = 0;
			}
		}
		//behind neighbor
		if (jStack[jNew] < yDim - 1 && statusArray[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] + 1, xDim)] == false) {
			if (imageDataPtr[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] + 1, xDim)] >= thres_min && imageDataPtr[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] + 1, xDim)] <= thres_max) {
				//the element is in the current region, then save its coordinates in temporary stacks
				iTmpStack.push_back(iStack[iNew]);
				jTmpStack.push_back(jStack[jNew] + 1);
				kTmpStack.push_back(kStack[kNew]);
			}
			else {
				//Not in region, update status and go to the next neighbor
				statusArray[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] + 1, xDim)] = true;
				segmentedImage[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] + 1, xDim)] = 0;
			}
		}

		//updating of processed element befor removal
		statusArray[kStack[kNew]][x_new(iStack[iNew], jStack[jNew], xDim)] = true;
		segmentedImage[kStack[kNew]][x_new(iStack[iNew], jStack[jNew], xDim)] = imageDataPtr[kStack[kNew]][x_new(iStack[iNew], jStack[jNew], xDim)];

		//Remove the processed element of the initial stacks
		kStack.pop_back();
		iStack.pop_back();
		jStack.pop_back();

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

		//End of big while loop
	}

	return true;
}
