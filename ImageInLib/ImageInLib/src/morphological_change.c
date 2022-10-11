#pragma warning(disable : 6011)
#pragma warning(disable : 4244)
#pragma warning(disable : 4477)
#pragma warning(disable : 4267)

#include<stdio.h>
#include<stdlib.h>
#include "../include/morphological_change.h"

//Erosion removes floating pixels and thin lines so that only substantive objects remain
bool erosion3D(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background) {
	if (imageDataPtr == NULL)
		return false;

	size_t i, j, k, cpt, dim2d = xDim * yDim;
	size_t xDim_ext = xDim + 2, yDim_ext = yDim + 2, zDim_ext = zDim + 2, dim2d_ext = xDim_ext * yDim_ext;

	dataType** extendedImage = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	dataType** maskImage = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	for (k = 0; k < zDim_ext; k++) {
		extendedImage[k] = (dataType*)malloc(dim2d_ext * sizeof(dataType));
		maskImage[k] = (dataType*)malloc(dim2d_ext * sizeof(dataType));
	}
	if (extendedImage == NULL || maskImage == NULL) return false;

	//Initialization
	for (k = 0; k < zDim_ext; k++) {
		for (i = 0; i < xDim_ext; i++) {
			for (j = 0; j < yDim_ext; j++) {
				extendedImage[k][x_new(i, j, xDim_ext)] = background;
				maskImage[k][x_new(i, j, xDim_ext)] = background;
			}
		}
	}
	copyDataToExtendedArea(imageDataPtr, extendedImage, zDim, xDim, yDim);
	copyDataToExtendedArea(imageDataPtr, maskImage, zDim, xDim, yDim);
	reflection3D(extendedImage, zDim_ext, xDim_ext, yDim_ext);

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

				if (cpt < 26) {
					maskImage[k][x_new(i, j, xDim_ext)] = background;
				}
			}
		}
	}
	copyDataToReducedArea(imageDataPtr, maskImage, zDim, xDim, yDim);

	//clean memory
	for (k = 0; k < zDim_ext; k++) {
		free(extendedImage[k]); 
		free(maskImage[k]);
	}
	free(extendedImage); free(maskImage);

	return true;
}

bool erosion3dHeighteenNeigbours(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background) {
	if (imageDataPtr == NULL)
		return false;

	size_t i, j, k, cpt, dim2d = xDim * yDim;
	size_t xDim_ext = xDim + 2, yDim_ext = yDim + 2, zDim_ext = zDim + 2, dim2d_ext = xDim_ext * yDim_ext;

	dataType** extendedImage = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	dataType** maskImage = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	for (k = 0; k < zDim_ext; k++) {
		extendedImage[k] = (dataType*)malloc(dim2d_ext * sizeof(dataType));
		maskImage[k] = (dataType*)malloc(dim2d_ext * sizeof(dataType));
	}
	if (extendedImage == NULL || maskImage == NULL) return false;
	//Initialization
	for (k = 0; k < zDim_ext; k++) {
		for (i = 0; i < xDim_ext; i++) {
			for (j = 0; j < yDim_ext; j++) {
				extendedImage[k][x_new(i, j, xDim_ext)] = background;
				maskImage[k][x_new(i, j, xDim_ext)] = background;
			}
		}
	}
	copyDataToExtendedArea(imageDataPtr, extendedImage, zDim, xDim, yDim);
	copyDataToExtendedArea(imageDataPtr, maskImage, zDim, xDim, yDim);
	reflection3D(extendedImage, zDim_ext, xDim_ext, yDim_ext);

	//neigbors check ---> 18 neigbors
	for (k = 1; k < zDim_ext - 1; k++) {
		for (i = 1; i < xDim_ext - 1; i++) {
			for (j = 1; j < yDim_ext - 1; j++) {

				cpt = 0;

				if (extendedImage[k - 1][x_new(i - 1, j + 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i, j, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i, j + 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i + 1, j - 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i + 1, j, xDim_ext)] != background) {
					cpt++;
				}

				if (extendedImage[k][x_new(i - 1, j - 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k][x_new(i - 1, j, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k][x_new(i - 1, j + 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k][x_new(i, j - 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k][x_new(i, j + 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k][x_new(i + 1, j - 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k][x_new(i + 1, j, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k][x_new(i + 1, j + 1, xDim_ext)] != background) {
					cpt++;
				}

				if (extendedImage[k + 1][x_new(i - 1, j + 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i, j, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i, j + 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i + 1, j - 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i + 1, j, xDim_ext)] != background) {
					cpt++;
				}

				if (cpt < 18) {
					maskImage[k][x_new(i, j, xDim_ext)] = background;
				}

			}
		}
	}
	copyDataToReducedArea(imageDataPtr, maskImage, zDim, xDim, yDim);

	for (k = 0; k < zDim_ext; k++) {
		free(extendedImage[k]);
		free(maskImage[k]);
	}
	free(extendedImage);
	free(maskImage);

	return true;
}

//Dilation makes objects more visible and fills in small holes in objects.
bool dilatation3D(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background) {
	if (imageDataPtr == NULL)
		return false;

	size_t i, j, k, cpt, dim2d = xDim * yDim;
	size_t xDim_ext = xDim + 2, yDim_ext = yDim + 2, zDim_ext = zDim + 2, dim2d_ext = xDim_ext * yDim_ext;

	dataType** extendedImage = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	dataType** maskImage = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	for (k = 0; k < zDim_ext; k++) {
		extendedImage[k] = (dataType*)malloc(dim2d_ext * sizeof(dataType));
		maskImage[k] = (dataType*)malloc(dim2d_ext * sizeof(dataType));
	}
	if (extendedImage == NULL || maskImage == NULL) return false;

	//Initialization
	for (k = 0; k < zDim_ext; k++) {
		for (i = 0; i < xDim_ext; i++) {
			for (j = 0; j < yDim_ext; j++) {
				extendedImage[k][x_new(i, j, xDim_ext)] = background;
				maskImage[k][x_new(i, j, xDim_ext)] = background;
			}
		}
	}
	copyDataToExtendedArea(imageDataPtr, extendedImage, zDim, xDim, yDim);
	reflection3D(extendedImage, zDim_ext, xDim_ext, yDim_ext);
	copyDataToAnotherArray(extendedImage, maskImage, zDim_ext, xDim_ext, yDim_ext);

	//neigbors check ---> 26 neighbors
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
					maskImage[k][x_new(i, j, xDim_ext)] = object;
				}
			}
		}
	}
	copyDataToReducedArea(imageDataPtr, maskImage, zDim, xDim, yDim);

	//clean memory
	for (k = 0; k < zDim_ext; k++) {
		free(extendedImage[k]);
		free(maskImage[k]);
	}
	free(extendedImage); free(maskImage);

	return true;
}

bool dilatation3dHeighteenNeigbours(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background) {
	if (imageDataPtr == NULL)
		return false;

	size_t i, j, k, cpt, dim2d = xDim * yDim;
	size_t xDim_ext = xDim + 2, yDim_ext = yDim + 2, zDim_ext = zDim + 2, dim2d_ext = xDim_ext * yDim_ext;

	dataType** extendedImage = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	dataType** maskImage = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	for (k = 0; k < zDim_ext; k++) {
		extendedImage[k] = (dataType*)malloc(dim2d_ext * sizeof(dataType));
		maskImage[k] = (dataType*)malloc(dim2d_ext * sizeof(dataType));
	}
	if (extendedImage == NULL || maskImage == NULL) return false;
	//Initialization
	for (k = 0; k < zDim_ext; k++) {
		for (i = 0; i < xDim_ext; i++) {
			for (j = 0; j < yDim_ext; j++) {
				extendedImage[k][x_new(i, j, xDim_ext)] = background;
				maskImage[k][x_new(i, j, xDim_ext)] = background;
			}
		}
	}
	copyDataToExtendedArea(imageDataPtr, extendedImage, zDim, xDim, yDim);
	reflection3D(extendedImage, zDim_ext, xDim_ext, yDim_ext);
	copyDataToAnotherArray(extendedImage, maskImage, zDim_ext, xDim_ext, yDim_ext);

	//neigbors check ---> 18 neigbors
	for (k = 1; k < zDim_ext - 1; k++) {
		for (i = 1; i < xDim_ext - 1; i++) {
			for (j = 1; j < yDim_ext - 1; j++) {

				cpt = 0;

				if (extendedImage[k - 1][x_new(i - 1, j + 1, xDim_ext)] == object) {
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

				if (extendedImage[k + 1][x_new(i - 1, j + 1, xDim_ext)] == object) {
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

				if (cpt > 0) {
					maskImage[k][x_new(i, j, xDim_ext)] = object;
				}

			}
		}
	}
	copyDataToReducedArea(imageDataPtr, maskImage, zDim, xDim, yDim);

	for (k = 0; k < zDim_ext; k++) {
		free(extendedImage[k]);
		free(maskImage[k]);
	}
	free(extendedImage);
	free(maskImage);

	return true;
}
