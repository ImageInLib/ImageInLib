#pragma warning(disable : 6011)
#pragma warning(disable : 4244)
#pragma warning(disable : 4477)
#pragma warning(disable : 4267)

#include <stdlib.h>
#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <omp.h>
#include <time.h>

#include "common_vtk.h"
#include "file.h"
#include "common_functions.h"
#include "../src/data_initialization.h"
#include "noiseGeneration.h"
#include "filtering.h"
#include "../src/data_load.h"
#include "../src/data_storage.h"
#include "../src/thresholding.h"
#include "Labelling.h"
#include "template_functions.h"
#include "../src/endianity_bl.h"
#include "../include/morphological_change.h"
#include "../src/shape_registration.h"
#include "../src/noise_generator.h"

#define originalMean 1095
#define offSet 0
#define margin 100
#define thresmin (originalMean + offSet - margin)
#define thresmax (originalMean + offSet + margin)
#define minimalSize 0


int main() {

	
	size_t i, j, k, xd;
	const size_t Length = 512;
	const size_t Width = 512;
	const size_t Height = 508;
	const size_t dim2D = Length * Width;

	short** image = (short**)malloc(Height * sizeof(short*));
	for (k = 0; k < Height; k++) {
		image[k] = (short*)malloc(dim2D * sizeof(short));
	}
	if (image == NULL) return false;

	//initialization
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				image[k][x_new(i, j, Length)] = 0;
			}
		}
	}

	//Loading
	load3dArrayRAW<short>(image, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/input/patient1b.raw");
	//store3dRawData<short>(image, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/loaded.raw");
	
	//Preparing Image container for float format
	Image_Data ImageData;
	ImageData.height = Height;
	ImageData.length = Length;
	ImageData.width = Width;
	ImageData.imageDataPtr = (dataType**)malloc(Height * sizeof(dataType*));
	dataType** inputImage = (dataType**)malloc(Height * sizeof(dataType*));
	for (k = 0;k < Height;k++) {
		ImageData.imageDataPtr[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		inputImage[k] = (dataType*)malloc(dim2D * sizeof(dataType));
	}
	if (ImageData.imageDataPtr == NULL || inputImage == NULL) return false;
	
	initialize3dArrayD(ImageData.imageDataPtr, Length, Width, Height, 0);
	initialize3dArrayD(inputImage, Length, Width, Height, 0);

	//copy input image in container
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				ImageData.imageDataPtr[k][x_new(i, j, Length)] = (dataType)image[k][x_new(i, j, Length)];
				inputImage[k][x_new(i, j, Length)] = (dataType)image[k][x_new(i, j, Length)];
			}
		}
	}
	
	////Find min and max values
	//dataType minData = 10000, maxData = -10000;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (ImageData.imageDataPtr[k][x_new(i, j, Length)] < minData)
	//				minData = ImageData.imageDataPtr[k][x_new(i, j, Length)];
	//			if (ImageData.imageDataPtr[k][x_new(i, j, Length)] > maxData)
	//				maxData = ImageData.imageDataPtr[k][x_new(i, j, Length)];
	//		}
	//	}
	//}

	//findClipBoxSingle(ImageData.imageDataPtr, Height, Length, Width);

	////Special histograms
	//dataType* HistoSliceBySlice = (dataType*)malloc(Height * sizeof(dataType));
	//for (k = 0; k < Height; k++) HistoSliceBySlice[k] = 0;
	//dataType nbElement;
	//for (k = 0; k < Height; k++) {
	//	nbElement = 0;
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (ImageData.imageDataPtr[k][x_new(i, j, Length)] >= thresmin && ImageData.imageDataPtr[k][x_new(i, j, Length)] <= thresmax) {
	//				nbElement++;
	//			}
	//		}
	//	}
	//	HistoSliceBySlice[k] = nbElement;
	//}

	//dataType* HistoRowByRow = (dataType*)malloc(Length * sizeof(dataType));
	//for (i = 0; i < Length; i++) HistoRowByRow[i] = 0;
	//for (i = 0; i < Length; i++) {
	//	nbElement = 0;
	//	for (k = 0; k < Height; k++) {
	//		for (j = 0; j < Width; j++) {
	//			if (ImageData.imageDataPtr[k][x_new(i, j, Length)] >= thresmin && ImageData.imageDataPtr[k][x_new(i, j, Length)] <= thresmax) {
	//				nbElement++;
	//			}
	//		}
	//	}
	//	HistoRowByRow[i] = nbElement;
	//}

	//dataType* HistoColumnByColumn = (dataType*)malloc(Length * sizeof(dataType));
	//for (i = 0; i < Length; i++) HistoColumnByColumn[i] = 0;
	//for (j = 0; j < Width; j++) {
	//	nbElement = 0;
	//	for (k = 0; k < Height; k++) {
	//		for (i = 0; i < Length; i++) {
	//			if (ImageData.imageDataPtr[k][x_new(i, j, Length)] >= thresmin && ImageData.imageDataPtr[k][x_new(i, j, Length)] <= thresmax) {
	//				nbElement++;
	//			}
	//		}
	//	}
	//	HistoColumnByColumn[j] = nbElement;
	//}

	//dataType maxHeight = 0, maxLenght = 0, maxWidth = 0;
	//size_t  k_peak = 0, i_peak = 0, j_peak = 0;

	//for (k = 0; k < Height; k++) {
	//	if (HistoSliceBySlice[k] > maxHeight) {
	//		maxHeight = HistoSliceBySlice[k];
	//		k_peak = k;
	//	}
	//}
	//for (i = 0; i < Length; i++) {
	//	if (HistoRowByRow[i] > maxLenght) {
	//		maxLenght = HistoRowByRow[i];
	//		i_peak = i;
	//	}
	//}
	//for (j = 0; j < Width; j++) {
	//	if (HistoColumnByColumn[j] > maxWidth) {
	//		maxWidth = HistoColumnByColumn[j];
	//		j_peak = j;
	//	}
	//}
	//const size_t height_new = 301, length_new = 301, width_new = 301;

	//dataType** croppedImage = (dataType**)malloc(height_new * sizeof(dataType*));
	//dataType** distanceMap = (dataType**)malloc(Height * sizeof(dataType*));
	//for (k = 0; k < height_new; k++) {
	//	croppedImage[k] = (dataType*)malloc(length_new * width_new * sizeof(dataType));
	//	distanceMap[k] = (dataType*)malloc(dim2D * sizeof(dataType));
	//}
	//if (croppedImage == NULL || distanceMap == NULL) return false;
	//initialize3dArrayD(croppedImage, length_new, width_new, height_new, 0);
	//initialize3dArrayD(distanceMap, length_new, width_new, height_new, 0);

	//size_t i_n = i_peak - 150, j_n = j_peak - 150, k_n = k_peak - 150;
	//for (k = 0; k < height_new; k++) {
	//	for (i = 0; i < length_new; i++) {
	//		for (j = 0; j < width_new; j++) {
	//			croppedImage[k][x_new(i, j, length_new)] = ImageData.imageDataPtr[k + k_n][x_new(i + i_n, j + j_n, Length)];
	//		}
	//	}
	//}
	//store3dRawData<dataType>(croppedImage, length_new, width_new, height_new, "C:/Users/Konan Allaly/Documents/Tests/output/cropped.raw");

	////Find min and max values
	//dataType minData = 10000, maxData = -10000;
	//for (k = 0; k < height_new; k++) {
	//	for (i = 0; i < length_new; i++) {
	//		for (j = 0; j < width_new; j++) {
	//			if (croppedImage[k][x_new(i, j, length_new)] < minData)
	//				minData = croppedImage[k][x_new(i, j, length_new)];
	//			if (croppedImage[k][x_new(i, j, length_new)] > maxData)
	//				maxData = croppedImage[k][x_new(i, j, length_new)];
	//		}
	//	}
	//}

	////Save Special Histogram
	//FILE* SpecialHisto;
	//if (fopen_s(&SpecialHisto, "C:/Users/Konan Allaly/Documents/Tests/output/Histogram.txt", "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//for (k = 0; k < Height; k++) {
	//	fprintf(SpecialHisto, "%d,", k);
	//	fprintf(SpecialHisto, "%f\n", HistoSliceBySlice[k]);
	//}
	//fclose(SpecialHisto);

	////Compute histogram
	//size_t totalClass = (size_t)(maxData - minData + 1);
	//size_t* histogram = (size_t*)malloc(totalClass * sizeof(size_t));
	//for (i = 0; i < totalClass; i++) histogram[i] = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			histogram[image[k][x_new(i, j, Length)]]++;
	//		}
	//	}
	//}

	////Save Histogram
	//FILE* histo;
	//if (fopen_s(&histo, "C:/Users/Konan Allaly/Documents/Tests/output/histoSlice304.txt", "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//for (i = 0; i < totalClass; i++) {
	//	fprintf(histo, "%d,", i);
	//	fprintf(histo, "%d\n", histogram[i]);
	//}
	//fclose(histo);

	////Compute probability
	//size_t numberOfCells = 0;
	//for (i = 0; i < totalClass; i++) {
	//	numberOfCells = numberOfCells + histogram[i];
	//}

	//dataType* Proba = (dataType*)malloc(totalClass * sizeof(dataType));
	//for (i = 0; i < totalClass; i++) Proba[i] = 0;

	//dataType sumProba = 0;
	//for (i = 0; i < totalClass; i++) {
	//	Proba[i] = (dataType)histogram[i]/numberOfCells;
	//	sumProba = sumProba + Proba[i]; // just for verification
	//}
	//
	//dataType* interClassVariance = (dataType*)malloc(totalClass * sizeof(dataType));
	//for (i = 0; i < totalClass; i++) interClassVariance[i] = 0;
	//size_t T, optimalThresholdValue = 0;
	//dataType sigma = 0, sum_weight;
	//dataType weightClass_b, meanClass_b, varClass_b;
	//dataType weightClass_f, meanClass_f, varClass_f;

	//for (T = 0; T < totalClass; T++) {
	//	weightClass_b = 0, meanClass_b = 0, varClass_b = 0;
	//	weightClass_f = 0, meanClass_f = 0, varClass_f = 0;
	//	//compute class weight
	//	for (i = 0; i <= T; i++) {
	//		weightClass_b = weightClass_b + Proba[i];
	//	}
	//	for (i = T + 1; i < totalClass; i++) {
	//		weightClass_f = weightClass_f + Proba[i];
	//	}
	//	sum_weight = weightClass_b + weightClass_f;
	//	//compute class mean
	//	for (i = 0; i <= T; i++) {
	//		meanClass_b = meanClass_b + i * Proba[i];
	//	}
	//	meanClass_b = meanClass_b / weightClass_b;
	//	for (i = T + 1; i < totalClass; i++) {
	//		meanClass_f = meanClass_f + i * Proba[i];
	//	}
	//	meanClass_f = meanClass_f / weightClass_f;
	//	//compute class variance
	//	for (i = 0; i <= T; i++) {
	//		varClass_b = varClass_b + (i - meanClass_b) * (i - meanClass_b) * Proba[i];
	//	}
	//	varClass_b = varClass_b / weightClass_b;
	//	for (i = T + 1; i < totalClass; i++) {
	//		varClass_f = varClass_f + (i - meanClass_f) * (i - meanClass_f) * Proba[i];
	//	}
	//	varClass_f = varClass_f / weightClass_f;
	//
	//	//compute inter-class variance
	//	sigma = weightClass_b * varClass_b + weightClass_f * varClass_f;
	//	//sigma = sqrt(sigma);
	//	interClassVariance[T] = sigma;
	//}
	//dataType minVariance = 1e10;
	//for (T = 0; T < totalClass; T++) {
	//	if (interClassVariance[T] < minVariance) {
	//		minVariance = interClassVariance[T];
	//	}
	//}
	//for (T = 0; T < totalClass; T++) {
	//	if (interClassVariance[T] == minVariance) {
	//		optimalThresholdValue = T;
	//	}
	//}
	//printf("optimal threshold value = %d \n", optimalThresholdValue);

	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (ImageData.imageDataPtr[k][x_new(i, j, Length)] < optimalThresholdValue) {
	//				ImageData.imageDataPtr[k][x_new(i, j, Length)] = optimalThresholdValue;
	//			}
	//			/*else {
	//				ImageData.imageDataPtr[k][x_new(i, j, Length)] = optimalThresholdValue;
	//			}*/
	//		}
	//	}
	//}
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/OTSUKV.raw");

	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			image[k][x_new(i, j, Length)] = (short)ImageData.imageDataPtr[k][x_new(i, j, Length)];
	//		}
	//	}
	//}
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			image[k][x_new(i, j, Length)] = (short)ImageData.imageDataPtr[k][x_new(i, j, Length)];
	//		}
	//	}
	//}
	//for (i = 0; i < totalClass; i++) {
	//	numberOfCells = numberOfCells + histogram[i];
	//}
	//for (i = 0; i < totalClass; i++) {
	//	Proba[i] = (dataType)histogram[i] / numberOfCells;
	//}

	//for (T = 0; T < totalClass; T++) {
	//	weightClass_b = 0, meanClass_b = 0, varClass_b = 0;
	//	weightClass_f = 0, meanClass_f = 0, varClass_f = 0;
	//	//compute class weight
	//	for (i = 0; i <= T; i++) {
	//		weightClass_b = weightClass_b + Proba[i];
	//	}
	//	for (i = T + 1; i < totalClass; i++) {
	//		weightClass_f = weightClass_f + Proba[i];
	//	}
	//	sum_weight = weightClass_b + weightClass_f;
	//	//compute class mean
	//	for (i = 0; i <= T; i++) {
	//		meanClass_b = meanClass_b + i * Proba[i];
	//	}
	//	meanClass_b = meanClass_b / weightClass_b;
	//	for (i = T + 1; i < totalClass; i++) {
	//		meanClass_f = meanClass_f + i * Proba[i];
	//	}
	//	meanClass_f = meanClass_f / weightClass_f;
	//	//compute class variance
	//	for (i = 0; i <= T; i++) {
	//		varClass_b = varClass_b + (i - meanClass_b) * (i - meanClass_b) * Proba[i];
	//	}
	//	varClass_b = varClass_b / weightClass_b;
	//	for (i = T + 1; i < totalClass; i++) {
	//		varClass_f = varClass_f + (i - meanClass_f) * (i - meanClass_f) * Proba[i];
	//	}
	//	varClass_f = varClass_f / weightClass_f;

	//	//compute inter-class variance
	//	sigma = weightClass_b * varClass_b + weightClass_f * varClass_f;
	//	//sigma = sqrt(sigma);
	//	interClassVariance[T] = sigma;
	//}
	//minVariance = 1e10;
	//for (T = 0; T < totalClass; T++) {
	//	if (interClassVariance[T] < minVariance) {
	//		minVariance = interClassVariance[T];
	//	}
	//}
	//for (T = 0; T < totalClass; T++) {
	//	if (interClassVariance[T] == minVariance) {
	//		optimalThresholdValue = T;
	//	}
	//}
	//printf("optimal threshold value = %d \n", optimalThresholdValue);

	//rescaleNewRange(ImageData.imageDataPtr, Length, Width, Height, 0, 1);
	//NoiseParameters Nparameters = { 0.05, 0.1, 10, 0, 1 };
	//const NoiseType Ntype = SALT_AND_PEPPER;
	//addNoiseToImage(ImageData.imageDataPtr, Length, Width, Height, Nparameters, Ntype);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "output/noisyImage.raw");

	////Filtering by geodesic mean curvature filter
	//const size_t maxNumberOfSolverIteration = 1000;
	//dataType  coef = 0.01, eps2 = 1e-4;
	//size_t numberOfTimeStep = 10;
	//Filter_Parameters GMC_filterParameters;
	//const FilterMethod methodFiltering = GEODESIC_MEAN_CURVATURE_FILTER;
	//dataType timeStepSize = 0.012, h = 1, sigma = 0.1, K = 0.018, omega_c = 1.5, tolerance = 5e-4;
	//size_t p = 1, timeStepsNum = 10, maxNumberOftimeSteps = 10;
	//GMC_filterParameters = { timeStepSize, h, sigma, K, omega_c, tolerance, eps2, p, timeStepsNum, maxNumberOfSolverIteration, maxNumberOftimeSteps };
	//rescaleNewRange(ImageData.imageDataPtr, Length, Width, Height, 0, 1);
	//filterImage(ImageData, GMC_filterParameters, maxNumberOfSolverIteration, coef, eps2, numberOfTimeStep, methodFiltering);
	//rescaleNewRange(ImageData.imageDataPtr, Length, Width, Height, minData, maxData);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/filtered_GMC.raw");

	////Filtering by mean curvature filter
	//const size_t maxNumberOfSolverIteration = 1000;
	//dataType  coef = 0.01, eps2 = 1e-4;
	//size_t numberOfTimeStep = 10;
	//Filter_Parameters MC_filterParameters;
	//const FilterMethod methodFiltering = MEAN_CURVATURE_FILTER;
	//dataType timeStepSize = 5, h = 1, sigma = 0.1, K = 0.018, omega_c = 1.5, tolerance = 5e-4;
	//size_t p = 1, timeStepsNum = 10, maxNumberOftimeSteps = 100;
	//MC_filterParameters = { timeStepSize, h, sigma, K, omega_c, tolerance, eps2, p, timeStepsNum, maxNumberOfSolverIteration, maxNumberOftimeSteps };
	//rescaleNewRange(ImageData.imageDataPtr, Length, Width, Height, 0, 1);
	//filterImage(ImageData, MC_filterParameters, maxNumberOfSolverIteration, coef, eps2, numberOfTimeStep, methodFiltering);
	//rescaleNewRange(ImageData.imageDataPtr, Length, Width, Height, minData, maxData);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "output/filtered_MC_10.raw");

	////Filtering Linear Heat Implicit
	//const size_t maxNumberOfSolverIteration = 0;
	//dataType  coef = 0, eps2 = 0;
	//size_t numberOfTimeStep = 0;
	//const FilterMethod methodFiltering = LINEAR_HEATEQUATION_IMPLICIT;
	//dataType timeStepSize = 1.2, h = 1, sigma = 0, K = 0, omega_c = 1.5, tolerance = 5e-4;
	//size_t p = 1, timeStepsNum = 12, maxNumberOftimeSteps = 0;
	//Filter_Parameters LHI_filterParameters{};
	//LHI_filterParameters = { timeStepSize, h, sigma, K, omega_c, tolerance, eps2, p, timeStepsNum, maxNumberOfSolverIteration, maxNumberOftimeSteps };
	//rescaleNewRange(ImageData.imageDataPtr, Length, Width, Height, 0, 1);
	//filterImage(ImageData, LHI_filterParameters, maxNumberOfSolverIteration, coef, eps2, numberOfTimeStep, methodFiltering);
	//rescaleNewRange(ImageData.imageDataPtr, Length, Width, Height, minData, maxData);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "output/filtered_LHI_12.raw");

	////Filtering Peronna Malick (Non linear implicit heat equation)
	//const size_t maxNumberOfSolverIteration = 0;
	//dataType  coef = 0, eps2 = 0;
	//size_t numberOfTimeStep = 0;
	//const FilterMethod methodFiltering = NONLINEAR_HEATEQUATION_IMPLICIT;
	//dataType timeStepSize = 1.2, h = 1, sigma = 0, K = 0.018, omega_c = 1.5, tolerance = 5e-4;
	//size_t p = 1, timeStepsNum = 10, maxNumberOftimeSteps = 0;
	//Filter_Parameters NLHI_filterParameters{};
	//NLHI_filterParameters = { timeStepSize, h, sigma, K, omega_c, tolerance, eps2, p, timeStepsNum, maxNumberOfSolverIteration, maxNumberOftimeSteps };
	//rescaleNewRange(ImageData.imageDataPtr, Length, Width, Height, 0, 1);
	//filterImage(ImageData, NLHI_filterParameters, maxNumberOfSolverIteration, coef, eps2, numberOfTimeStep, methodFiltering);
	//rescaleNewRange(ImageData.imageDataPtr, Length, Width, Height, minData, maxData);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "output/filtered_NLHI_10.raw");
	
	////Remove fat 
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (ImageData.imageDataPtr[k][x_new(i, j, Length)] >= 900 && ImageData.imageDataPtr[k][x_new(i, j, Length)] <= 950) {
	//				ImageData.imageDataPtr[k][x_new(i, j, Length)] = maxData;
	//			}
	//		}
	//	}
	//}
	// 
	//Thresholding

	//dataType vmin = thresmin; char pathSave[300]; int n;
	//for ( n = 0; n < 150; n++) {
	//	thresholding3dFunctionN(ImageData.imageDataPtr, Length, Width, Height, vmin, thresmax, minData, maxData);
	//	vmin++;
	//	sprintf_s(pathSave, "C:/Users/Konan Allaly/Documents/Tests/output/threshold3d00%d.raw", n);
	//	store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, pathSave);
	//	for (k = 0; k < Height; k++) {
	//		for (i = 0; i < Length; i++) {
	//			for (j = 0; j < Width; j++) {
	//				ImageData.imageDataPtr[k][x_new(i, j, Length)] = (dataType)image[k][x_new(i, j, Length)];
	//			}
	//		}
	//	}
	//}
	//dataType vmin = thresmin + 100;

	//thresholding3dFunctionN(ImageData.imageDataPtr, Length, Width, Height, thresmin, thresmax, minData, maxData);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/thres100.raw");
	//dilatation3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/dilate100.raw");
	//erosion3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/ErosionAfterd.raw");

	//for (i = 0; i < 15; i++) {
	//	erosion3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//}
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/Erosion20Times.raw");

	//thresholding3dFunctionN(croppedImage, length_new, width_new, height_new, thresmin, thresmax, minData, maxData);
	//bruteForceFunction_3D(distanceMap, croppedImage, length_new, width_new, height_new, 10000000, minData);
	//rouyTourinFunction_3D(distanceMap, ImageData.imageDataPtr, 0.5, Length, Width, Height, 0.4, 1.0);
	//store3dRawData<dataType>(distanceMap, length_new, width_new, height_new, "C:/Users/Konan Allaly/Documents/Tests/output/distanceCropped.raw");

	//char pathSave[300]; dataType dist_max = 30;
	//do {
	//	for (k = 0; k < Height; k++) {
	//		for (i = 0; i < Length; i++) {
	//			for (j = 0; j < Width; j++) {
	//				if (distanceMap[k][x_new(i, j, Length)] < dist_max && (inputImage[k][x_new(i, j, Length)] >= thresmin && inputImage[k][x_new(i, j, Length)] <= thresmax)) {
	//					inputImage[k][x_new(i, j, Length)] = minData;
	//				}
	//				else {
	//					inputImage[k][x_new(i, j, Length)] = maxData;
	//				}
	//			}
	//		}
	//	}
	//	sprintf_s(pathSave, "C:/Users/Konan Allaly/Documents/Tests/output/segmented3d0%d.raw", (int)dist_max);
	//	store3dRawData<dataType>(inputImage, Length, Width, Height, pathSave);
	//	for (k = 0; k < Height; k++) {
	//		for (i = 0; i < Length; i++) {
	//			for (j = 0; j < Width; j++) {
	//				inputImage[k][x_new(i, j, Length)] = (dataType)image[k][x_new(i, j, Length)];
	//			}
	//		}
	//	}
	//	dist_max = dist_max + 1;
	//} while (dist_max < 101);
	
	
	////Loop for Erosion
	//char pathSaveErosion[200]; int n;
	//for (n = 0; n < 11; n++) {
	//	erosion3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//	//dilatation3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//	sprintf_s(pathSaveErosion, "C:/Users/Konan Allaly/Documents/Tests/output/erosion0%d.raw", n);
	//	store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, pathSaveErosion);
	//}

	////Loop for Dilatation
	//char pathSaveDilatation[200];
	//for (n = 1; n < 6; n++) {
	//	dilatation3D(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//	sprintf_s(pathSaveErosion, "C:/Users/Konan Allaly/Documents/Tests/output/F_dilatation0%d.raw", n);
	//	store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, pathSaveErosion);
	//}

	////find the Centroid
	//dataType* cenTroid = (dataType*)malloc(3 * sizeof(dataType));
	//centroidImage(imageData, cenTroid, Height, Length, Width, 0);
	//printf("Coordinates of centroid  : ");
	//for (i = 0; i < 3; i++) {
	//	printf("%f, ", cenTroid[i]);
	//}
	//printf("\n");

	//dilatation3D(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//dilatation3D(ImageData.imageDataPtr, Length, Width, Height, minData, maxData);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "output/DilatationAfterErosion30.raw");

	////Labelling
	//int** labelArray = (int**)malloc(Height * sizeof(int*));
	//bool** status = (bool**)malloc(Height * sizeof(bool*));
	//for (k = 0; k < Height; k++) {
	//	labelArray[k] = (int*)malloc(dim2D * sizeof(int));
	//	status[k] = (bool*)malloc(dim2D * sizeof(bool));
	//}
	//if (labelArray == NULL || status == NULL) return false;

	////initialization
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			labelArray[k][x_new(i, j, Length)] = 0;
	//			status[k][x_new(i, j, Length)] = false;
	//		}
	//	}
	//}
	//double start = clock();
	//labelling3D(ImageData.imageDataPtr, labelArray, status, Length, Width, Height, maxData);
	//double finish = clock();
	//printf("Execution time for the new code : %.3lf \n", (finish - start) / CLOCKS_PER_SEC);

	////number of region cells
	//int numberOfRegionsCells = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (ImageData.imageDataPtr[k][x_new(i, j, Length)] == maxData) {
	//				numberOfRegionsCells++;
	//			}
	//		}
	//	}
	//}
	//printf("Number of Regions Cells : %d \n", numberOfRegionsCells);

	////Counting
	//int* countingArray = (int*)malloc(numberOfRegionsCells * sizeof(int));
	//if (countingArray == NULL) return false;
	//for (i = 0; i < numberOfRegionsCells; i++) countingArray[i] = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (labelArray[k][x_new(i, j, Length)] > 0) {
	//				countingArray[labelArray[k][x_new(i, j, Length)]]++;
	//			}
	//		}
	//	}
	//}

	////Remove small regions 
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (countingArray[labelArray[k][x_new(i, j, Length)]] < minimalSize) {
	//				countingArray[labelArray[k][x_new(i, j, Length)]] = 0;
	//			}
	//		}
	//	}
	//}

	////Verification
	//int numberOfVoxelsSegmented = 0;
	//for (i = 0; i < numberOfRegionsCells; i++) {
	//	if (countingArray[i] > 0) {
	//		numberOfVoxelsSegmented = numberOfVoxelsSegmented + countingArray[i];
	//	}
	//}
	//printf("Number of Regions Cells : %d \n", numberOfRegionsCells);
	//printf("Number of Voxels segmented : %d \n", numberOfVoxelsSegmented);

	////Number of regions
	//int numberOfRegions = 0;
	//for (i = 0; i < numberOfRegionsCells; i++) {
	//	if (countingArray[i] > 0) {
	//		numberOfRegions++;
	//	}
	//}
	//printf("Number of regions = %d \n", numberOfRegions);

	////Statistics
	//FILE* file;
	//if (fopen_s(&file, "C:/Users/Konan Allaly/Documents/Tests/output/Statistics.txt", "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(file, " 3D real image with %d Slices \n", Height);
	//fprintf(file, " Labelling function takes = %.5lf seconds\n", (finish - start) / CLOCKS_PER_SEC);
	//fprintf(file, " Number of Objects Cells = %d \n", numberOfRegionsCells);
	//fprintf(file, " Number of Regions  = %d \n", numberOfRegions);
	//fprintf(file, "#############################################################\n\n");
	//fprintf(file, " Regions with more than = %d voxels\n", minimalSize);
	//fprintf(file, "\n\n");
	//fprintf(file, "Region Label        Count");
	//for (i = 0; i < numberOfRegionsCells; i++) {
	//	if (countingArray[i] >= minimalSize) {
	//		fprintf(file, "\n");
	//		fprintf(file, "   %d                 %d", i, countingArray[i]);
	//	}
	//}
	//fclose(file);

	////Remove
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (countingArray[labelArray[k][x_new(i, j, Length)]] < minimalSize) {
	//				labelArray[k][x_new(i, j, Length)] = 0;
	//			}
	//		}
	//	}
	//}

	////Saving with template function
	//store3dRawData<int>(labelArray, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/segmented.raw");

	////Keep one region 
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (labelArray[k][x_new(i, j, Length)] != 93) {
	//				labelArray[k][x_new(i, j, Length)] = 0;
	//			}
	//		}
	//	}
	//}
	//store3dRawData<int>(labelArray, Length, Width, Height, "output/region93.raw");

	////Find the centroid of keepped region
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			ImageData.imageDataPtr[k][x_new(i, j, Length)] = (dataType)labelArray[k][x_new(i, j, Length)];
	//		}
	//	}
	//}

	//dataType* cenTroid = (dataType*)malloc(3 * sizeof(dataType));
	//centroidImage(imageData, cenTroid, Height, Length, Width, 0);
	//printf("Coordinates of centroid  : ");
	//for (i = 0; i < 3; i++) {
	//	printf("%f, ", cenTroid[i]);
	//}
	//printf("\n");

	////Delatation on region 
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			imageData[k][x_new(i, j, Length)] = (dataType)labelArray[k][x_new(i, j, Length)];
	//		}
	//	}
	//}
	//dilatation3D(imageData, Length, Width, Height, 5185, background);
	//dilatation3D(imageData, Length, Width, Height, 5185, background);
	//dilatation3D(imageData, Length, Width, Height, 5185, background);
	//store3dRawData<dataType>(imageData, Length, Width, Height, "output/region5185WithDilatation.raw");

	//free memory
	for (k = 0; k < Height; k++) {
		free(image[k]); 
		free(ImageData.imageDataPtr[k]);
		free(inputImage[k]);
		//free(labelArray[k]); free(status[k]);
	}
	free(image); free(ImageData.imageDataPtr); free(inputImage);
	//free(HistoSliceBySlice); free(HistoRowByRow); free(HistoColumnByColumn);

	//for (k = 0; k < height_new; k++) {
	//	free(distanceMap[k]); free(croppedImage);
	//}
	//free(croppedImage); free(distanceMap);

	//free(Proba); free(histogram); free(interClassVariance);
	//free(labelArray); free(status); free(countingArray);
	//free(cenTroid);
	
	return EXIT_SUCCESS;
}
