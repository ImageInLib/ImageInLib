/*
* Author: Konan ALLALY
* Purpose: INFLANET project - Image Processing in Nuclear Medicine (2D/3D)
* Language:  C/C++
*/
#pragma warning(disable : 6011)
#pragma warning(disable : 4244)
#pragma warning(disable : 4477)
#pragma warning(disable : 4267)

#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include <common_vtk.h>
#include "heat_equation.h"
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

				if (extendedImage[k][x_new(i, j, xDim_ext)] != background) {

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

/* //Compute gradient
bool imageGradient(dataType ** imageDataPtr, const char* pathSavePtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType h) {
	
	if (imageDataPtr == NULL)
		return false;

	size_t i, j, k, dim2d = xDim * yDim;

	dataType ux, uy, uz, quotient = 4 * h;
	dataType u, uN, uS, uE, uW, uNW, uNE, uSE, uSW, Tu, TuN, TuS, TuE, TuW, TuNW, TuNE, TuSE, TuSW, Bu, BuN, BuS, BuE, BuW, BuNW, BuNE, BuSE, BuSW;

	dataType ** n_Ptr = (dataType**)malloc(zDim * sizeof(dataType*));
	dataType ** s_Ptr = (dataType**)malloc(zDim * sizeof(dataType*));
	dataType ** w_Ptr = (dataType**)malloc(zDim * sizeof(dataType*));
	dataType ** e_Ptr = (dataType**)malloc(zDim * sizeof(dataType*));
	dataType ** b_Ptr = (dataType**)malloc(zDim * sizeof(dataType*));
	dataType ** t_Ptr = (dataType**)malloc(zDim * sizeof(dataType*));

	for (k = 0; k < zDim; k++) {
		n_Ptr[k] = (dataType*)malloc(dim2d * sizeof(dataType));
		s_Ptr[k] = (dataType*)malloc(dim2d * sizeof(dataType));
		w_Ptr[k] = (dataType*)malloc(dim2d * sizeof(dataType));
		e_Ptr[k] = (dataType*)malloc(dim2d * sizeof(dataType));
		b_Ptr[k] = (dataType*)malloc(dim2d * sizeof(dataType));
		t_Ptr[k] = (dataType*)malloc(dim2d * sizeof(dataType));
	}
	if (n_Ptr == NULL || s_Ptr == NULL || w_Ptr == NULL || w_Ptr == NULL || b_Ptr == NULL || t_Ptr == NULL)
		return false;

	const size_t xDim_ext = xDim + 2, yDim_ext = yDim + 2, zDim_ext = zDim + 2;

	dataType ** extendedImagePtr = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	for (k = 0; k < zDim_ext; k++) {
		extendedImagePtr[k] = (dataType*)malloc(xDim_ext * yDim_ext * sizeof(dataType));
	}
	if (extendedImagePtr == NULL)
		return false;

	copyDataToExtendedArea(imageDataPtr, extendedImagePtr, zDim, xDim, yDim);
	reflection3D(extendedImagePtr, zDim_ext, xDim_ext, yDim_ext);

	size_t k_ext, i_ext, j_ext;

	for (k = 0, k_ext = 1; k < zDim; k++, k_ext++) {
		for (i = 0, i_ext = 1; i < xDim; i++, i_ext++) {
			for (j = 0, j_ext = 1; j < yDim; j++, j_ext++) {

				u = extendedImagePtr[k_ext][x_new(i_ext, j_ext, xDim_ext)];
				uN = extendedImagePtr[k_ext][x_new(i_ext, j_ext - 1, xDim_ext)];
				uS = extendedImagePtr[k_ext][x_new(i_ext, j_ext + 1, xDim_ext)];
				uE = extendedImagePtr[k_ext][x_new(i_ext + 1, j_ext, xDim_ext)];
				uW = extendedImagePtr[k_ext][x_new(i_ext - 1, j_ext, xDim_ext)];
				uNW = extendedImagePtr[k_ext][x_new(i_ext - 1 , j_ext - 1, xDim_ext)];
				uNE = extendedImagePtr[k_ext][x_new(i_ext + 1, j_ext - 1, xDim_ext)];
				uSE = extendedImagePtr[k_ext][x_new(i_ext + 1, j_ext + 1, xDim_ext)];
				uSW = extendedImagePtr[k_ext][x_new(i_ext - 1, j_ext + 1, xDim_ext)];

				Tu = extendedImagePtr[k_ext - 1][x_new(i_ext, j_ext, xDim_ext)];
				TuN = extendedImagePtr[k_ext - 1][x_new(i_ext, j_ext - 1, xDim_ext)];
				TuS = extendedImagePtr[k_ext - 1][x_new(i_ext, j_ext + 1, xDim_ext)];
				TuE = extendedImagePtr[k_ext - 1][x_new(i_ext + 1, j_ext, xDim_ext)];
				TuW = extendedImagePtr[k_ext - 1][x_new(i_ext - 1, j_ext, xDim_ext)];
				TuNW = extendedImagePtr[k_ext - 1][x_new(i_ext - 1, j_ext - 1, xDim_ext)];
				TuNE = extendedImagePtr[k_ext - 1][x_new(i_ext + 1, j_ext - 1, xDim_ext)];
				TuSE = extendedImagePtr[k_ext - 1][x_new(i_ext + 1, j_ext + 1, xDim_ext)];
				TuSW = extendedImagePtr[k_ext - 1][x_new(i_ext - 1, j_ext + 1, xDim_ext)];

				Bu = extendedImagePtr[k_ext + 1][x_new(i_ext, j_ext, xDim_ext)];
				BuN = extendedImagePtr[k_ext + 1][x_new(i_ext, j_ext - 1, xDim_ext)];
				BuS = extendedImagePtr[k_ext + 1][x_new(i_ext, j_ext + 1, xDim_ext)];
				BuE = extendedImagePtr[k_ext + 1][x_new(i_ext + 1, j_ext, xDim_ext)];
				BuW = extendedImagePtr[k_ext + 1][x_new(i_ext - 1, j_ext, xDim_ext)];
				BuNW = extendedImagePtr[k_ext + 1][x_new(i_ext - 1, j_ext - 1, xDim_ext)];
				BuNE = extendedImagePtr[k_ext + 1][x_new(i_ext + 1, j_ext - 1, xDim_ext)];
				BuSE = extendedImagePtr[k_ext + 1][x_new(i_ext + 1, j_ext + 1, xDim_ext)];
				BuSW = extendedImagePtr[k_ext + 1][x_new(i_ext - 1, j_ext + 1, xDim_ext)];

				//calculation of coefficients on the image data

				// Calculation of coefficients in east direction
				ux = (uE - u) / h;
				uy = ((uN + uNE) - (uS + uSE)) / quotient;
				uz = ((Tu + TuE) - (Bu + BuE))/ quotient;
				e_Ptr[k][x_new(i, j, xDim)] = sqrt((ux * ux) + (uy * uy) + (uz * uz));

				// Calculation of coefficients in west direction
				ux = (uW - u) / h;
				uy = ((uNW + uN) - (uSW + uS)) / quotient;
				uz = ((TuW + Tu) - (BuW + Bu)) / quotient;
				w_Ptr[k][x_new(i, j, xDim)] = sqrt( (ux * ux) + (uy * uy) + (uz * uz));

				// Calculation of coefficients in north direction
				ux = ((uNE + uE) - (uNW + uW)) / quotient;
				uy = (uN - u) / h;
				uz = ((TuN + Tu) - (BuN + Bu)) / quotient;
				n_Ptr[k][x_new(i, j, xDim)] = sqrt((ux * ux) + (uy * uy) + (uz * uz));

				// Calculation of coefficients in south direction
				ux = ((uE + uSE) - (uW + uSW)) / quotient;
				uy = (uS - u) / h;
				uz = ((TuS + Tu) - (BuS + Bu)) / quotient;
				s_Ptr[k][x_new(i, j, xDim)] = sqrt((ux * ux) + (uy * uy) + (uz * uz));

				// Calculation of coefficients in top direction
				ux = ((TuE + uE) - (TuW + uW)) / quotient;
				uy = ((TuN + uN) - (TuS + uS)) / quotient;
				uz = (Tu - u) / h;
				t_Ptr[k][x_new(i, j, xDim)] = sqrt((ux * ux) + (uy * uy) + (uz * uz));

				// Calculation of coefficients in bottom direction
				ux = ((BuW + uW) - (BuE + uE)) / quotient;
				uy = ((BuN + uN) - (BuS + uS)) / quotient;
				uz = (Bu - u) / h;
				b_Ptr[k][x_new(i, j, xDim)] = sqrt((ux * ux) + (uy * uy) + (uz * uz));
			}
		}
	}

	Vtk_File_Info* vtkInfo = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	vtkInfo->spacing[0] = 1.0; vtkInfo->spacing[1] = 1.0; vtkInfo->spacing[2] = 1.0;
	vtkInfo->origin[0] = 0.0; vtkInfo->origin[1] = 0.0; vtkInfo->origin[2] = 0.0;
	vtkInfo->dimensions[0] = xDim; vtkInfo->dimensions[1] = yDim; vtkInfo->dimensions[2] = zDim;
	vtkInfo->vDataType = dta_Flt; vtkInfo->operation = copyTo;
	vtkDataForm dataForm = dta_binary;

	//name construction
	unsigned char name[350];
	unsigned char name_ending[100];
	const char* pathSaveVTK;
	
	strcpy_s(name, sizeof name, pathSavePtr);
	sprintf_s(name_ending, sizeof(name_ending), "_g_north.vtk", i);
	strcat_s(name, sizeof(name), name_ending);
	pathSaveVTK = name;
	vtkInfo->dataPointer = n_Ptr;
	storeVtkFile(pathSaveVTK, vtkInfo, dataForm);

	strcpy_s(name, sizeof name, pathSavePtr);
	sprintf_s(name_ending, sizeof(name_ending), "_g_south.vtk", i);
	strcat_s(name, sizeof(name), name_ending);
	pathSaveVTK = name;
	vtkInfo->dataPointer = s_Ptr;
	storeVtkFile(pathSaveVTK, vtkInfo, dataForm);

	strcpy_s(name, sizeof name, pathSavePtr);
	sprintf_s(name_ending, sizeof(name_ending), "_g_west.vtk", i);
	strcat_s(name, sizeof(name), name_ending);
	pathSaveVTK = name;
	vtkInfo->dataPointer = w_Ptr;
	storeVtkFile(pathSaveVTK, vtkInfo, dataForm);

	strcpy_s(name, sizeof name, pathSavePtr);
	sprintf_s(name_ending, sizeof(name_ending), "_g_est.vtk", i);
	strcat_s(name, sizeof(name), name_ending);
	pathSaveVTK = name;
	vtkInfo->dataPointer = e_Ptr;
	storeVtkFile(pathSaveVTK, vtkInfo, dataForm);

	strcpy_s(name, sizeof name, pathSavePtr);
	sprintf_s(name_ending, sizeof(name_ending), "_g_top.vtk", i);
	strcat_s(name, sizeof(name), name_ending);
	pathSaveVTK = name;
	vtkInfo->dataPointer = t_Ptr;
	storeVtkFile(pathSaveVTK, vtkInfo, dataForm);

	strcpy_s(name, sizeof name, pathSavePtr);
	sprintf_s(name_ending, sizeof(name_ending), "_g_bottom.vtk", i);
	strcat_s(name, sizeof(name), name_ending);
	pathSaveVTK = name;
	vtkInfo->dataPointer = b_Ptr;
	storeVtkFile(pathSaveVTK, vtkInfo, dataForm);

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
*/

/*
bool edgesDetector(dataType** imageDataPtr, const char* pathSavePtr, const size_t xDim, const size_t yDim, const size_t zDim, Filter_Parameters smoothParameters) {

	if (imageDataPtr == NULL)
		return false;

	size_t i, j, k, dim2d = xDim * yDim;

	dataType ux, uy, uz, quotient = 4 * smoothParameters.h;
	dataType u, uN, uS, uE, uW, uNW, uNE, uSE, uSW, Tu, TuN, TuS, TuE, TuW, TuNW, TuNE, TuSE, TuSW, Bu, BuN, BuS, BuE, BuW, BuNW, BuNE, BuSE, BuSW;

	dataType** n_Ptr = (dataType**)malloc(zDim * sizeof(dataType*));
	dataType** s_Ptr = (dataType**)malloc(zDim * sizeof(dataType*));
	dataType** w_Ptr = (dataType**)malloc(zDim * sizeof(dataType*));
	dataType** e_Ptr = (dataType**)malloc(zDim * sizeof(dataType*));
	dataType** b_Ptr = (dataType**)malloc(zDim * sizeof(dataType*));
	dataType** t_Ptr = (dataType**)malloc(zDim * sizeof(dataType*));

	for (k = 0; k < zDim; k++) {
		n_Ptr[k] = (dataType*)malloc(dim2d * sizeof(dataType));
		s_Ptr[k] = (dataType*)malloc(dim2d * sizeof(dataType));
		w_Ptr[k] = (dataType*)malloc(dim2d * sizeof(dataType));
		e_Ptr[k] = (dataType*)malloc(dim2d * sizeof(dataType));
		b_Ptr[k] = (dataType*)malloc(dim2d * sizeof(dataType));
		t_Ptr[k] = (dataType*)malloc(dim2d * sizeof(dataType));
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

	presmoothData.imageDataPtr = extendedImagePtr;
	heatImplicitScheme(presmoothData, smoothParameters);

	size_t k_ext, i_ext, j_ext;

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

	//name construction
	unsigned char name[350];
	unsigned char name_ending[100];
	const char* pathSaveVTK;

	strcpy_s(name, sizeof name, pathSavePtr);
	sprintf_s(name_ending, sizeof(name_ending), "_ED_north.vtk", i);
	strcat_s(name, sizeof(name), name_ending);
	pathSaveVTK = name;
	vtkInfo->dataPointer = n_Ptr;
	storeVtkFile(pathSaveVTK, vtkInfo, dataForm);

	strcpy_s(name, sizeof name, pathSavePtr);
	sprintf_s(name_ending, sizeof(name_ending), "_ED_south.vtk", i);
	strcat_s(name, sizeof(name), name_ending);
	pathSaveVTK = name;
	vtkInfo->dataPointer = s_Ptr;
	storeVtkFile(pathSaveVTK, vtkInfo, dataForm);

	strcpy_s(name, sizeof name, pathSavePtr);
	sprintf_s(name_ending, sizeof(name_ending), "_ED_west.vtk", i);
	strcat_s(name, sizeof(name), name_ending);
	pathSaveVTK = name;
	vtkInfo->dataPointer = w_Ptr;
	storeVtkFile(pathSaveVTK, vtkInfo, dataForm);

	strcpy_s(name, sizeof name, pathSavePtr);
	sprintf_s(name_ending, sizeof(name_ending), "_ED_est.vtk", i);
	strcat_s(name, sizeof(name), name_ending);
	pathSaveVTK = name;
	vtkInfo->dataPointer = e_Ptr;
	storeVtkFile(pathSaveVTK, vtkInfo, dataForm);

	strcpy_s(name, sizeof name, pathSavePtr);
	sprintf_s(name_ending, sizeof(name_ending), "_ED_top.vtk", i);
	strcat_s(name, sizeof(name), name_ending);
	pathSaveVTK = name;
	vtkInfo->dataPointer = t_Ptr;
	storeVtkFile(pathSaveVTK, vtkInfo, dataForm);

	strcpy_s(name, sizeof name, pathSavePtr);
	sprintf_s(name_ending, sizeof(name_ending), "_ED_bottom.vtk", i);
	strcat_s(name, sizeof(name), name_ending);
	pathSaveVTK = name;
	vtkInfo->dataPointer = b_Ptr;
	storeVtkFile(pathSaveVTK, vtkInfo, dataForm);

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
*/