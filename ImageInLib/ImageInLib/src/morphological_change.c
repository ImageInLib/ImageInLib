#pragma warning(disable : 6011)
#pragma warning(disable : 4244)
#pragma warning(disable : 4477)
#pragma warning(disable : 4267)

#include<stdio.h>
#include<stdlib.h>
#include "../include/morphological_change.h"

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

//Erosion ---> Shrink object and remove boundaries pixels
bool erosion3D(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background) {
	if (imageDataPtr == NULL)
		return false;

	size_t i, j, k, cpt, dim2d = xDim * yDim;
	size_t xDim_ext = xDim + 2, yDim_ext = yDim + 2, zDim_ext = zDim + 2, dim2d_ext = xDim_ext * yDim_ext;

	dataType** maskImage = (dataType**)malloc(zDim * sizeof(dataType*));
	for (k = 0; k < zDim; k++)
		maskImage[k] = (dataType*)malloc(dim2d * sizeof(dataType));
	if (maskImage == NULL) return false;
	//Initialization
	for (k = 0; k < zDim; k++) {
		for (i = 0; i < xDim; i++) {
			for (j = 0; j < yDim; j++) {
				maskImage[k][x_new(i, j, xDim)] = background;
			}
		}
	}

	dataType** extendedImage = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	for (k = 0; k < zDim_ext; k++)
		extendedImage[k] = (dataType*)malloc(dim2d_ext * sizeof(dataType));
	if (extendedImage == NULL) return false;
	//Initialization
	for (k = 0; k < zDim_ext; k++) {
		for (i = 0; i < xDim_ext; i++) {
			for (j = 0; j < yDim_ext; j++) {
				extendedImage[k][x_new(i, j, xDim_ext)] = background;
			}
		}
	}
	copyDataToExtendedArea(imageDataPtr, extendedImage, zDim, xDim, yDim);
	reflection3D(extendedImage, zDim_ext, xDim_ext, yDim_ext);
	
	////Copy data to extended area
	//for (k = 0; k < zDim; k++) {
	//	for (i = 0; i < xDim; i++) {
	//		for (j = 0; j < yDim; j++) {
	//			reflectedImage[k + 1][x_new(i + 1, j + 1, xDim)] = imageDataPtr[k][x_new(i, j, xDim)];
	//		}
	//	}
	//}
	//// Y reflection
	//for (k = 0; k < zDim + 2; k++){
	//	for (i = 0; i < xDim + 2; i++){
	//		reflectedImage[k][x_new(i, 0, xDim + 2)] = reflectedImage[k][x_new(i, 1, xDim + 2)];
	//		reflectedImage[k][x_new(i, yDim + 1, xDim + 2)] = reflectedImage[k][x_new(i, yDim, xDim + 2)];
	//	}
	//}
	//// X Direction
	//for (k = 0; k < zDim + 2; k++){
	//	for (j = 0; j < yDim + 2; j++){
	//		reflectedImage[k][x_new(0, j, xDim + 2)] = reflectedImage[k][x_new(1, j, xDim + 2)];
	//		reflectedImage[k][x_new(xDim + 1, j, xDim + 2)] = reflectedImage[k][x_new(xDim, j, xDim + 2)];
	//	}
	//}
	//// Z Reflection
	//for (i = 0; i < xDim + 2; i++){
	//	for (j = 0; j < yDim + 2; j++) {
	//		reflectedImage[0][x_new(i, j, xDim + 2)] = reflectedImage[1][x_new(i, j, xDim + 2)];
	//		reflectedImage[zDim + 1][x_new(i, j, xDim + 2)] = reflectedImage[zDim][x_new(i, j, xDim + 2)];
	//	}
	//}

	//neigbors check
	for (k = 1; k < zDim_ext - 1; k++) {
		for (i = 1; i < xDim_ext - 1; i++) {
			for (j = 1; j < yDim_ext - 1; j++) {

				cpt = 0;

				if (extendedImage[k - 1][x_new(i - 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i - 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i - 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i + 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i + 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i + 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}

				if (extendedImage[k][x_new(i - 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i - 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i - 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i + 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i + 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i + 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}

				if (extendedImage[k + 1][x_new(i - 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i - 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i - 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i + 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i + 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i + 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}

				if (cpt == 26) {
					maskImage[k - 1][x_new(i - 1, j - 1, xDim)] = object;
				}
				else {
					maskImage[k - 1][x_new(i - 1, j - 1, xDim)] = background;
				}

			}
		}
	}

	//copy
	for (k = 0; k < zDim; k++) {
		for (i = 0; i < xDim; i++) {
			for (j = 0; j < yDim; j++) {
				imageDataPtr[k][x_new(i, j, xDim)] = maskImage[k][x_new(i, j, xDim)];
			}
		}
	}

	//clean memory
	for (k = 0; k < zDim; k++)
		free(maskImage[k]);
	free(maskImage);

	for (k = 0; k < zDim + 2; k++)
		free(extendedImage[k]);
	free(extendedImage);

	return true;
}

//Dilatation ---> Enlarge object and fill holes
bool dilatation3D(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background) {
	if (imageDataPtr == NULL)
		return false;

	size_t i, j, k, cpt, dim2d = xDim * yDim;
	size_t xDim_ext = xDim + 2, yDim_ext = yDim + 2, zDim_ext = zDim + 2, dim2d_ext = xDim_ext * yDim_ext;

	dataType** maskImage = (dataType**)malloc(zDim * sizeof(dataType*));
	for (k = 0; k < zDim; k++)
		maskImage[k] = (dataType*)malloc(dim2d * sizeof(dataType));
	if (maskImage == NULL) return false;
	//Initialization
	for (k = 0; k < zDim; k++) {
		for (i = 0; i < xDim; i++) {
			for (j = 0; j < yDim; j++) {
				maskImage[k][x_new(i, j, xDim)] = background;
			}
		}
	}

	dataType** extendedImage = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	for (k = 0; k < zDim_ext; k++)
		extendedImage[k] = (dataType*)malloc(dim2d_ext * sizeof(dataType));
	if (extendedImage == NULL) return false;
	//Initialization
	for (k = 0; k < zDim_ext; k++) {
		for (i = 0; i < xDim_ext; i++) {
			for (j = 0; j < yDim_ext; j++) {
				extendedImage[k][x_new(i, j, xDim_ext)] = background;
			}
		}
	}
	copyDataToExtendedArea(imageDataPtr, extendedImage, zDim, xDim, yDim);
	reflection3D(extendedImage, zDim_ext, xDim_ext, yDim_ext);

	////Copy data to extended area
	//for (k = 0; k < zDim; k++) {
	//	for (i = 0; i < xDim; i++) {
	//		for (j = 0; j < yDim; j++) {
	//			reflectedImage[k + 1][x_new(i + 1, j + 1, xDim)] = imageDataPtr[k][x_new(i, j, xDim)];
	//		}
	//	}
	//}
	//// Y reflection
	//for (k = 0; k < zDim + 2; k++) {
	//	for (i = 0; i < xDim + 2; i++) {
	//		reflectedImage[k][x_new(i, 0, xDim + 2)] = reflectedImage[k][x_new(i, 1, xDim + 2)];
	//		reflectedImage[k][x_new(i, yDim + 1, xDim + 2)] = reflectedImage[k][x_new(i, yDim, xDim + 2)];
	//	}
	//}
	//// X Direction
	//for (k = 0; k < zDim + 2; k++) {
	//	for (j = 0; j < yDim + 2; j++) {
	//		reflectedImage[k][x_new(0, j, xDim + 2)] = reflectedImage[k][x_new(1, j, xDim + 2)];
	//		reflectedImage[k][x_new(xDim + 1, j, xDim + 2)] = reflectedImage[k][x_new(xDim, j, xDim + 2)];
	//	}
	//}
	//// Z Reflection
	//for (i = 0; i < xDim + 2; i++) {
	//	for (j = 0; j < yDim + 2; j++) {
	//		reflectedImage[0][x_new(i, j, xDim + 2)] = reflectedImage[1][x_new(i, j, xDim + 2)];
	//		reflectedImage[zDim + 1][x_new(i, j, xDim + 2)] = reflectedImage[zDim][x_new(i, j, xDim + 2)];
	//	}
	//}

	//neigbors check
	for (k = 1; k < zDim_ext - 1; k++) {
		for (i = 1; i < xDim_ext - 1; i++) {
			for (j = 1; j < yDim_ext - 1; j++) {

				cpt = 0;

				if (extendedImage[k - 1][x_new(i - 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i - 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i - 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i + 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i + 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i + 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}

				if (extendedImage[k][x_new(i - 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i - 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i - 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i + 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i + 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i + 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}

				if (extendedImage[k + 1][x_new(i - 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i - 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i - 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i + 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i + 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i + 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}

				if (cpt > 0) {
					maskImage[k - 1][x_new(i - 1, j - 1, xDim)] = object;
				}
			}
		}
	}

	//copy
	for (k = 0; k < zDim; k++) {
		for (i = 0; i < xDim; i++) {
			for (j = 0; j < yDim; j++) {
				imageDataPtr[k][x_new(i, j, xDim)] = maskImage[k][x_new(i, j, xDim)];
			}
		}
	}

	//clean memory
	for (k = 0; k < zDim; k++)
		free(maskImage[k]);
	free(maskImage);

	for (k = 0; k < zDim + 2; k++)
		free(extendedImage[k]);
	free(extendedImage);

	return true;
}

