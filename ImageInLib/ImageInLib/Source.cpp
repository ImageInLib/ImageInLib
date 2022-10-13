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
#define minimalSize 2000


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
	dataType** distanceMap = (dataType**)malloc(Height * sizeof(dataType*));
	for (k = 0;k < Height;k++) {
		ImageData.imageDataPtr[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		inputImage[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		distanceMap[k] = (dataType*)malloc(dim2D * sizeof(dataType));
	}
	if (ImageData.imageDataPtr == NULL || inputImage == NULL || distanceMap == NULL) return false;
	
	initialize3dArrayD(ImageData.imageDataPtr, Length, Width, Height, 0);
	initialize3dArrayD(inputImage, Length, Width, Height, 0);
	initialize3dArrayD(distanceMap, Length, Width, Height, 0);

	//copy input image in container
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				ImageData.imageDataPtr[k][x_new(i, j, Length)] = (dataType)image[k][x_new(i, j, Length)];
				inputImage[k][x_new(i, j, Length)] = (dataType)image[k][x_new(i, j, Length)];
			}
		}
	}
	
	//Find min and max values
	dataType minData = 10000, maxData = -10000;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (ImageData.imageDataPtr[k][x_new(i, j, Length)] < minData)
					minData = ImageData.imageDataPtr[k][x_new(i, j, Length)];
				if (ImageData.imageDataPtr[k][x_new(i, j, Length)] > maxData)
					maxData = ImageData.imageDataPtr[k][x_new(i, j, Length)];
			}
		}
	}

	//Filtering by geodesic mean curvature filter
	const size_t maxNumberOfSolverIteration = 1000;
	dataType  coef = 0.01, eps2 = 1e-4;
	size_t numberOfTimeStep = 100;
	Filter_Parameters GMC_filterParameters;
	const FilterMethod methodFiltering = GEODESIC_MEAN_CURVATURE_FILTER;
	dataType timeStepSize = 1e-4, h = 10, sigma = 0.1, K = 1000, omega_c = 1.4, tolerance = 1e-10;
	size_t p = 1, timeStepsNum = 100, maxNumberOftimeSteps = 100;
	GMC_filterParameters = { timeStepSize, h, sigma, K, omega_c, tolerance, eps2, p, timeStepsNum, maxNumberOfSolverIteration, maxNumberOftimeSteps };
	rescaleNewRange(ImageData.imageDataPtr, Length, Width, Height, 0, 1);
	filterImage(ImageData, GMC_filterParameters, maxNumberOfSolverIteration, coef, eps2, numberOfTimeStep, methodFiltering);
	rescaleNewRange(ImageData.imageDataPtr, Length, Width, Height, minData, maxData);
	store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/filtered_GMC.raw");


	////Remove fat tissues, bones , water, kidney, white matter
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			//Remove fat tissues
	//			if (ImageData.imageDataPtr[k][x_new(i, j, Length)] >= 850 && ImageData.imageDataPtr[k][x_new(i, j, Length)] <= 1000) {
	//				ImageData.imageDataPtr[k][x_new(i, j, Length)] = minData;
	//			}
	//			//bones
	//			if (ImageData.imageDataPtr[k][x_new(i, j, Length)] >= 2000) {
	//				ImageData.imageDataPtr[k][x_new(i, j, Length)] = minData;
	//			}
	//			//water
	//			if (ImageData.imageDataPtr[k][x_new(i, j, Length)] == 1000) {
	//				ImageData.imageDataPtr[k][x_new(i, j, Length)] = minData;
	//			}
	//			//kidney
	//			if (ImageData.imageDataPtr[k][x_new(i, j, Length)] == 1030) {
	//				ImageData.imageDataPtr[k][x_new(i, j, Length)] = minData;
	//			}
	//			//white matter
	//			if (ImageData.imageDataPtr[k][x_new(i, j, Length)] >= 970 && ImageData.imageDataPtr[k][x_new(i, j, Length)] <= 980) {
	//				ImageData.imageDataPtr[k][x_new(i, j, Length)] = minData;
	//			}
	//		}
	//	}
	//}
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/RemovedSeveralElementsPM.raw");
	
	////segment a ball in the liver's centre of gravity
	//const size_t x_liver = 258, y_liver = 157, z_liver = 288;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (sqrt((x_liver - i) * (x_liver - i) + (y_liver - j) * (y_liver - j) + (z_liver - k) * (z_liver - k)) > 20) {
	//				image[k][x_new(i, j, Length)] = minData;
	//			}
	//		}
	//	}
	//}
	//store3dRawData<short>(image, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/Ball_Liver.raw");

	//erosion3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/ErosionAfterRemoving.raw");

	//thresholding3dFunctionN(ImageData.imageDataPtr, Length, Width, Height, thresmin, thresmax, minData, maxData);
	//for (i = 0; i < 6; i++) {
	//	erosion3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//}
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/ErosionAfterThreshAndErosion.raw");

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
	//if (fopen_s(&file, "C:/Users/Konan Allaly/Documents/Tests/output/StatsE3.txt", "w") != 0) {
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

	////Save bigest region 
	//int maxVolume = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (countingArray[labelArray[k][x_new(i, j, Length)]] > maxVolume) {
	//				maxVolume = countingArray[labelArray[k][x_new(i, j, Length)]];
	//			}
	//		}
	//	}
	//}
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (countingArray[labelArray[k][x_new(i, j, Length)]] != maxVolume) {
	//				ImageData.imageDataPtr[k][x_new(i, j, Length)] = minData;
	//			}
	//		}
	//	}
	//}
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/bigestRegionErosion7.raw");

	//char pathSaveDil[300];
	//for (int i = 0; i < 7; i++) {
	//	dilatation3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//	sprintf_s(pathSaveDil, "C:/Users/Konan Allaly/Documents/Tests/output/dilatation0%d.raw", i);
	//	store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, pathSaveDil);
	//}
	//fastSweepingFunction_3D(distanceMap, ImageData.imageDataPtr, Length, Width, Height, 1, 100000, minData);
	////store3dRawData<dataType>(distanceMap, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/distanceMap.raw");
	//dataType distanceMax = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (distanceMap[k][x_new(i, j, Length)] > distanceMax) {
	//				distanceMax = distanceMap[k][x_new(i, j, Length)];
	//			}
	//		}
	//	}
	//}

	//dataType cptMax = 0;
	//int x_max = 0, y_max = 0, z_max = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (distanceMap[k][x_new(i, j, Length)] == distanceMax) {
	//				cptMax++; x_max = i; y_max = j; z_max = k;
	//			}
	//		}
	//	}
	//}
	//printf("The maximal distance is = %f, and %f points have that distance \n", distanceMax, cptMax);
	//printf("The coordinates are : (%d, %d, %d)", x_max, y_max, z_max);

	//free memory
	for (k = 0; k < Height; k++) {
		free(image[k]); 
		free(ImageData.imageDataPtr[k]);
		free(inputImage[k]); free(distanceMap[k]);
		//free(labelArray[k]); free(status[k]);
	}
	free(image); free(ImageData.imageDataPtr); free(inputImage); free(distanceMap);
	//free(labelArray); free(status); free(countingArray);
	
	//free(HistoSliceBySlice); free(HistoRowByRow); free(HistoColumnByColumn);
	//for (k = 0; k < height_new; k++) {
	//	free(distanceMap[k]); free(croppedImage);
	//}
	//free(croppedImage); free(distanceMap);
	//free(Proba); free(histogram); free(interClassVariance);
	//free(cenTroid);
	
	return EXIT_SUCCESS;
}
