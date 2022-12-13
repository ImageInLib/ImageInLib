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


bool regionGrowing(dataType** imageDataPtr, int** segmentedImage, bool** statusArray, const size_t xDim, const size_t yDim, const size_t zDim, dataType thres_min, dataType thres_max, Filter_Parameters smoothParameters) {

	if (imageDataPtr == NULL || segmentedImage == NULL || statusArray == NULL)
		return false;

	vector<size_t> iStack, jStack, kStack, iTmpStack, jTmpStack, kTmpStack;
	size_t i = 0, j = 0, k = 0, iNew = 0, jNew = 0, kNew = 0, n = 0, dim2D = xDim * yDim;
	int label = 1;

	dataType ux, uy, uz, quotient = 4 * smoothParameters.h;
	dataType u, uN, uS, uE, uW, uNW, uNE, uSE, uSW, Tu, TuN, TuS, TuE, TuW, TuNW, TuNE, TuSE, TuSW, Bu, BuN, BuS, BuE, BuW, BuNW, BuNE, BuSE, BuSW;

	dataType** n_Ptr = (dataType**)malloc(zDim * sizeof(dataType*));
	dataType** s_Ptr = (dataType**)malloc(zDim * sizeof(dataType*));
	dataType** w_Ptr = (dataType**)malloc(zDim * sizeof(dataType*));
	dataType** e_Ptr = (dataType**)malloc(zDim * sizeof(dataType*));
	dataType** b_Ptr = (dataType**)malloc(zDim * sizeof(dataType*));
	dataType** t_Ptr = (dataType**)malloc(zDim * sizeof(dataType*));

	for (k = 0; k < zDim; k++) {
		n_Ptr[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		s_Ptr[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		w_Ptr[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		e_Ptr[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		b_Ptr[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		t_Ptr[k] = (dataType*)malloc(dim2D * sizeof(dataType));
	}
	if (n_Ptr == NULL || s_Ptr == NULL || w_Ptr == NULL || w_Ptr == NULL || b_Ptr == NULL || t_Ptr == NULL)
		return false;

	const size_t xDim_ext = xDim + 2, yDim_ext = yDim + 2, zDim_ext = zDim + 2;

	dataType** extendedImagePtr = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	for (k = 0; k < zDim_ext; k++) {
		extendedImagePtr[k] = (dataType*)malloc(xDim_ext * yDim_ext * sizeof(dataType));
	}
	if (extendedImagePtr == NULL)
		return false;

	Image_Data presmoothData; presmoothData.height = zDim_ext; presmoothData.length = xDim_ext; presmoothData.width = yDim_ext;

	copyDataToExtendedArea(imageDataPtr, extendedImagePtr, zDim, xDim, yDim);
	reflection3D(extendedImagePtr, zDim_ext, xDim_ext, yDim_ext);
	rescaleNewRange(extendedImagePtr, xDim_ext, yDim_ext, zDim_ext, 0, 1);
	presmoothData.imageDataPtr = extendedImagePtr;
	imageGradient(extendedImagePtr, "C:/Users/Konan Allaly/Documents/Tests/output/", xDim_ext, yDim_ext, zDim_ext, 1.0);
	heatImplicitScheme(presmoothData, smoothParameters);

	size_t k_ext, i_ext, j_ext;
	//Edges detector coef
	for (k = 0, k_ext = 1; k < zDim; k++, k_ext++) {
		for (i = 0, i_ext = 1; i < xDim; i++, i_ext++) {
			for (j = 0, j_ext = 1; j < yDim; j++, j_ext++) {

				u = presmoothData.imageDataPtr[k_ext][x_new(i_ext, j_ext, xDim_ext)];
				uN = presmoothData.imageDataPtr[k_ext][x_new(i_ext, j_ext - 1, xDim_ext)];
				uS = presmoothData.imageDataPtr[k_ext][x_new(i_ext, j_ext + 1, xDim_ext)];
				uE = presmoothData.imageDataPtr[k_ext][x_new(i_ext + 1, j_ext, xDim_ext)];
				uW = presmoothData.imageDataPtr[k_ext][x_new(i_ext - 1, j_ext, xDim_ext)];
				uNW = presmoothData.imageDataPtr[k_ext][x_new(i_ext - 1, j_ext - 1, xDim_ext)];
				uNE = presmoothData.imageDataPtr[k_ext][x_new(i_ext + 1, j_ext - 1, xDim_ext)];
				uSE = presmoothData.imageDataPtr[k_ext][x_new(i_ext + 1, j_ext + 1, xDim_ext)];
				uSW = presmoothData.imageDataPtr[k_ext][x_new(i_ext - 1, j_ext + 1, xDim_ext)];

				Tu = presmoothData.imageDataPtr[k_ext - 1][x_new(i_ext, j_ext, xDim_ext)];
				TuN = presmoothData.imageDataPtr[k_ext - 1][x_new(i_ext, j_ext - 1, xDim_ext)];
				TuS = presmoothData.imageDataPtr[k_ext - 1][x_new(i_ext, j_ext + 1, xDim_ext)];
				TuE = presmoothData.imageDataPtr[k_ext - 1][x_new(i_ext + 1, j_ext, xDim_ext)];
				TuW = presmoothData.imageDataPtr[k_ext - 1][x_new(i_ext - 1, j_ext, xDim_ext)];
				TuNW = presmoothData.imageDataPtr[k_ext - 1][x_new(i_ext - 1, j_ext - 1, xDim_ext)];
				TuNE = presmoothData.imageDataPtr[k_ext - 1][x_new(i_ext + 1, j_ext - 1, xDim_ext)];
				TuSE = presmoothData.imageDataPtr[k_ext - 1][x_new(i_ext + 1, j_ext + 1, xDim_ext)];
				TuSW = presmoothData.imageDataPtr[k_ext - 1][x_new(i_ext - 1, j_ext + 1, xDim_ext)];

				Bu = presmoothData.imageDataPtr[k_ext + 1][x_new(i_ext, j_ext, xDim_ext)];
				BuN = presmoothData.imageDataPtr[k_ext + 1][x_new(i_ext, j_ext - 1, xDim_ext)];
				BuS = presmoothData.imageDataPtr[k_ext + 1][x_new(i_ext, j_ext + 1, xDim_ext)];
				BuE = presmoothData.imageDataPtr[k_ext + 1][x_new(i_ext + 1, j_ext, xDim_ext)];
				BuW = presmoothData.imageDataPtr[k_ext + 1][x_new(i_ext - 1, j_ext, xDim_ext)];
				BuNW = presmoothData.imageDataPtr[k_ext + 1][x_new(i_ext - 1, j_ext - 1, xDim_ext)];
				BuNE = presmoothData.imageDataPtr[k_ext + 1][x_new(i_ext + 1, j_ext - 1, xDim_ext)];
				BuSE = presmoothData.imageDataPtr[k_ext + 1][x_new(i_ext + 1, j_ext + 1, xDim_ext)];
				BuSW = presmoothData.imageDataPtr[k_ext + 1][x_new(i_ext - 1, j_ext + 1, xDim_ext)];

				//calculation of coefficients on the image data

				// Calculation of coefficients in east direction
				ux = (uE - u) / smoothParameters.h;
				uy = ((uN + uNE) - (uS + uSE)) / quotient;
				uz = ((Tu + TuE) - (Bu + BuE)) / quotient;
				e_Ptr[k][x_new(i, j, xDim)] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), smoothParameters.edge_detector_coefficient);

				// Calculation of coefficients in west direction
				ux = (uW - u) / smoothParameters.h;
				uy = ((uNW + uN) - (uSW + uS)) / quotient;
				uz = ((TuW + Tu) - (BuW + Bu)) / quotient;
				w_Ptr[k][x_new(i, j, xDim)] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), smoothParameters.edge_detector_coefficient);

				// Calculation of coefficients in north direction
				ux = ((uNE + uE) - (uNW + uW)) / quotient;
				uy = (uN - u) / smoothParameters.h;
				uz = ((TuN + Tu) - (BuN + Bu)) / quotient;
				n_Ptr[k][x_new(i, j, xDim)] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), smoothParameters.edge_detector_coefficient);

				// Calculation of coefficients in south direction
				ux = ((uE + uSE) - (uW + uSW)) / quotient;
				uy = (uS - u) / smoothParameters.h;
				uz = ((TuS + Tu) - (BuS + Bu)) / quotient;
				s_Ptr[k][x_new(i, j, xDim)] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), smoothParameters.edge_detector_coefficient);

				// Calculation of coefficients in top direction
				ux = ((TuE + uE) - (TuW + uW)) / quotient;
				uy = ((TuN + uN) - (TuS + uS)) / quotient;
				uz = (Tu - u) / smoothParameters.h;
				t_Ptr[k][x_new(i, j, xDim)] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), smoothParameters.edge_detector_coefficient);

				// Calculation of coefficients in bottom direction
				ux = ((BuW + uW) - (BuE + uE)) / quotient;
				uy = ((BuN + uN) - (BuS + uS)) / quotient;
				uz = (Bu - u) / smoothParameters.h;
				b_Ptr[k][x_new(i, j, xDim)] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), smoothParameters.edge_detector_coefficient);
			}
		}
	}

	Vtk_File_Info* vtkInfo = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	vtkInfo->spacing[0] = 1.0; vtkInfo->spacing[1] = 1.0; vtkInfo->spacing[2] = 1.0;
	vtkInfo->origin[0] = 0.0; vtkInfo->origin[1] = 0.0; vtkInfo->origin[2] = 0.0;
	vtkInfo->dimensions[0] = xDim; vtkInfo->dimensions[1] = yDim; vtkInfo->dimensions[2] = zDim;
	vtkInfo->vDataType = dta_Flt; vtkInfo->operation = copyTo;
	vtkDataForm dataForm = dta_binary;
	const char* pathSaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/_edge_g_top.vtk";

	vtkInfo->dataPointer = t_Ptr;
	storeVtkFile(pathSaveVTK, vtkInfo, dataForm);


	//statusArray has false everywhere in the beginning
	//false---> non-processed and true---> processed
	//segmentedImage has 0 everywhere
	for (k = 0; k < zDim; k++) {
		for (i = 0; i < xDim; i++) {
			for (j = 0; j < yDim; j++) {

				if (statusArray[k][x_new(i, j, xDim)] == false) {

					if (imageDataPtr[k][x_new(i, j, xDim)] >= thres_min && imageDataPtr[k][x_new(i, j, xDim)] <= thres_max) {

						//top neighbor
						if (k > 0 && statusArray[k - 1][x_new(i, j, xDim)] == false) {
							if (imageDataPtr[k - 1][x_new(i, j, xDim)] >= thres_min && imageDataPtr[k - 1][x_new(i, j, xDim)] <= thres_max) {
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
							if (imageDataPtr[k + 1][x_new(i, j, xDim)] >= thres_min && imageDataPtr[k + 1][x_new(i, j, xDim)] <= thres_max) {
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
							if (imageDataPtr[k][x_new(i - 1, j, xDim)] >= thres_min && imageDataPtr[k][x_new(i - 1, j, xDim)] <= thres_max) {
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
							if (imageDataPtr[k][x_new(i + 1, j, xDim)] >= thres_min && imageDataPtr[k][x_new(i + 1, j, xDim)] <= thres_max) {
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
							if (imageDataPtr[k][x_new(i, j - 1, xDim)] >= thres_min && imageDataPtr[k][x_new(i, j - 1, xDim)] <= thres_max) {
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
							if (imageDataPtr[k][x_new(i, j + 1, xDim)] >= thres_min && imageDataPtr[k][x_new(i, j + 1, xDim)] <= thres_max) {
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
						//If the current element isn't on edge we give its label, and update its status
						if (t_Ptr[k][x_new(i, j, xDim)] <= 0.95) {
							segmentedImage[k][x_new(i, j, xDim)] = label;
						}
						statusArray[k][x_new(i, j, xDim)] = true;

						//Now we check neighbors of its neighbors saved in the stacks
						//We start by the last added element in stacks
						//If there is no neighbor iStack.size() = jStack.size() = kStack.size() = 0 , and the while loop will no be ran
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
									//Not in region, update status and go to the next neighbor
									statusArray[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] - 1, xDim)] = true;
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
								}
							}

							//updating of processed element befor removal
							if (t_Ptr[kStack[kNew]][x_new(iStack[iNew], jStack[jNew], xDim)] <= 0.95) {
								segmentedImage[kStack[kNew]][x_new(iStack[iNew], jStack[jNew], xDim)] = label;
							}
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


	for (k = 0; k < zDim; k++) {
		free(n_Ptr[k]);
		free(s_Ptr[k]);
		free(w_Ptr[k]);
		free(e_Ptr[k]);
		free(b_Ptr[k]);
		free(t_Ptr[k]);
	}
	free(n_Ptr);
	free(s_Ptr);
	free(w_Ptr);
	free(e_Ptr);
	free(b_Ptr);
	free(t_Ptr);

	for (k = 0; k < zDim_ext; k++) {
		free(extendedImagePtr[k]);
	}
	free(extendedImagePtr);

	free(vtkInfo);

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
//=========================


