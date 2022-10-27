#pragma warning(disable : 6011)
#pragma warning(disable : 4244)
#pragma warning(disable : 4477)
#pragma warning(disable : 4267)
#pragma warning(disable : 4700)

#include <stdlib.h>
#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <time.h>

#include "common_vtk.h"
#include "file.h"
#include "common_functions.h"
#include "../src/data_initialization.h"
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
#include "../src/segmentation3D_subsurf.h"

#define originalMean 71.1245
#define offSet 1024
#define standarDeviation 22.001
#define thresmin 995//(originalMean + offSet - 2*standarDeviation)
#define thresmax 1213//(originalMean + offSet + 2*standarDeviation)
#define minimalSize 2000


int main() {

	size_t i, j, k;
	const size_t Length = 512;
	const size_t Width = 512;
	const size_t Height = 406;
	const size_t dim2D = Length * Width;

	short** image = (short**)malloc(Height * sizeof(short*));
	for (k = 0; k < Height; k++) {
		image[k] = (short*)malloc(dim2D * sizeof(short));
		//imageDataVTK[k] = (int*)malloc(dim2D * sizeof(int));
	}
	if (image == NULL ) return false;

	//initialization
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				image[k][x_new(i, j, Length)] = 0;
			}
		}
	}

	//Loading
	load3dArrayRAW<short>(image, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/input/patient2.raw");
	//store3dRawData<short>(image, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/loaded.raw");

	//Preparing Image container for float format
	Image_Data ImageData; ImageData.height = Height; ImageData.length = Length; ImageData.width = Width;
	ImageData.imageDataPtr = (dataType**)malloc(Height * sizeof(dataType*));
	dataType** distanceMap = (dataType**)malloc(Height * sizeof(dataType*));
	for (k = 0;k < Height;k++) {
		ImageData.imageDataPtr[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		distanceMap[k] = (dataType*)malloc(dim2D * sizeof(dataType));
	}
	if (ImageData.imageDataPtr == NULL || distanceMap == NULL) return false;
	
	initialize3dArrayD(ImageData.imageDataPtr, Length, Width, Height, 0);
	initialize3dArrayD(distanceMap, Length, Width, Height, 0);

	//copy input image in container
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				ImageData.imageDataPtr[k][x_new(i, j, Length)] = (dataType)image[k][x_new(i, j, Length)];
			}
		}
	}

	//Storage_Flags flags = { false, false };
	//unsigned char name[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/loaded.vtk";
	//store3dDataVtkD(ImageData.imageDataPtr, Length, Width, Height, name, 1, flags);
	//store3dDataArrayD(ImageData.imageDataPtr, Length, Width, Height, name, flags);

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

	//thresholdingOTSU(ImageData.imageDataPtr, Length, Width, Height, minData, maxData);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/thresh.raw");

	////Filtering by linear heat equation filter
	//dataType timeStepSize = 0.5, h = 1, tol = 1e-4, omega_c = 1.5;
	//size_t p = 1, timeStepsNum = 5, iter_max = 100;
	//Filter_Parameters LH_filterParameters;
	//const FilterMethod methodFiltering = LINEAR_HEATEQUATION_EXPLICIT;
	//LH_filterParameters = { timeStepSize, h, 0, 0, omega_c, tol, 0, 0, p, timeStepsNum, iter_max};
	//rescaleNewRange(ImageData.imageDataPtr, Length, Width, Height, 0, 1);
	//filterImage(ImageData, LH_filterParameters, methodFiltering);
	//rescaleNewRange(ImageData.imageDataPtr, Length, Width, Height, minData, maxData);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/filtered_ex_M.raw");

	//dataType* Histogram = (dataType*)malloc((int)(maxData - minData + 1) * sizeof(dataType));
	//if (Histogram == NULL) return false;
	//for (i = 0; i < maxData - minData + 1; i++) Histogram[i] = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			Histogram[image[k][x_new(i, j, Length)]]++;
	//		}
	//	}
	//}
	//dataType meanLiver = 0, maxFrequence = 0;
	//for (i = 600; i < maxData - minData + 1; i++) {
	//	if (Histogram[i] > maxFrequence) {
	//		maxFrequence = Histogram[i];
	//		meanLiver = i;
	//	}
	//}
	//printf("The mean of the liver is = %f", meanLiver);

	////Remove fat tissues, bones , water, kidney, white matter
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			//Remove fat tissues
	//			if (ImageData.imageDataPtr[k][x_new(i, j, Length)] >= 850 && ImageData.imageDataPtr[k][x_new(i, j, Length)] <= 1000) {
	//				ImageData.imageDataPtr[k][x_new(i, j, Length)] = minData;
	//			}
	//			//bones
	//			if (ImageData.imageDataPtr[k][x_new(i, j, Length)] >= 1150 && ImageData.imageDataPtr[k][x_new(i, j, Length)] <= 1250) {
	//				ImageData.imageDataPtr[k][x_new(i, j, Length)] = minData;
	//			}
	//			////water
	//			//if (ImageData.imageDataPtr[k][x_new(i, j, Length)] == 1000) {
	//			//	ImageData.imageDataPtr[k][x_new(i, j, Length)] = minData;
	//			//}
	//			////kidney
	//			//if (ImageData.imageDataPtr[k][x_new(i, j, Length)] == 1030) {
	//			//	ImageData.imageDataPtr[k][x_new(i, j, Length)] = minData;
	//			//}
	//			////white matter
	//			//if (ImageData.imageDataPtr[k][x_new(i, j, Length)] >= 970 && ImageData.imageDataPtr[k][x_new(i, j, Length)] <= 980) {
	//			//	ImageData.imageDataPtr[k][x_new(i, j, Length)] = minData;
	//			//}
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

	thresholding3dFunctionN(ImageData.imageDataPtr, Length, Width, Height, thresmin, thresmax, minData, maxData);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/ThreshP1b.raw");
	//for (i = 0; i < 6; i++) {
	//	erosion3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//}
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/Erosion.raw");

	//erosion3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//erosion3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//erosion3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//erosion3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//dilatation3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//erosion3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//erosion3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//dilatation3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//dilatation3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//erosion3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//erosion3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//dilatation3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//erosion3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/ErosionP3.raw");

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
	////store3dRawData<int>(labelArray, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/segmented.raw");
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
	////Save biggest region 
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
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/bigestRegionP3V2.raw");

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

	//char pathSaveDil[300];
	//for (int i = 0; i < 7; i++) {
	//	dilatation3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//	sprintf_s(pathSaveDil, "C:/Users/Konan Allaly/Documents/Tests/output/dilatation0%d.raw", i);
	//	store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, pathSaveDil);
	//}

	fastSweepingFunction_3D(distanceMap, ImageData.imageDataPtr, Length, Width, Height, 1, 100000000, minData);
	//store3dRawData<dataType>(distanceMap, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/distanceMap.raw");
	dataType distanceMax = -1;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (distanceMap[k][x_new(i, j, Length)] >= distanceMax) {
					distanceMax = distanceMap[k][x_new(i, j, Length)];
				}
			}
		}
	}
	dataType cptMax = 0;
	int i_max = 0, j_max = 0, k_max = 0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (distanceMap[k][x_new(i, j, Length)] == distanceMax) {
					cptMax++; 
					i_max = (int)i; j_max = (int)j; k_max = (int)k;
				}
			}
		}
	}
	printf("The maximal distance is = %f, and %f points have that distance \n", distanceMax, cptMax);
	printf("The coordinates are : (%d, %d, %d)", j_max, i_max, k_max);

	////save according distance
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (distanceMap[k][x_new(i, j, Length)] > 2.5) {
	//				ImageData.imageDataPtr[k][x_new(i, j, Length)] = maxData;
	//			}
	//			else {
	//				ImageData.imageDataPtr[k][x_new(i, j, Length)] = minData;
	//			}
	//		}
	//	}
	//}
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/segment2.raw");

	//dilatation3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//dilatation3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//dilatation3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//dilatation3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/dilateBigestV2.raw");

	//// Centroid of the segmented Liver
	//dataType* cenTroid = (dataType*)malloc(3 * sizeof(dataType));
	//centroidImage(ImageData.imageDataPtr, cenTroid, Height, Length, Width, 0);
	//printf("Coordinates of centroid  : ");
	//for (i = 0; i < 3; i++) {
	//	printf("%f, ", cenTroid[i]);
	//}
	//printf("\n");

	////Ball arround the highest distance point
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (sqrt( (i_max - i)* (i_max - i) + (j_max - j) * (j_max - j) + (k_max - k) * (k_max - k) ) < 15 ) {
	//				ImageData.imageDataPtr[k][x_new(i, j, Length)] = maxData;
	//			}
	//			else {
	//				ImageData.imageDataPtr[k][x_new(i, j, Length)] = minData;
	//			}
	//		}
	//	}
	//}
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/ball.raw");

	////Saving
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (ImageData.imageDataPtr[k][x_new(i, j, Length)] == maxData) {
	//				inputImage[k][x_new(i, j, Length)] = minData;
	//			}
	//			if (inputImage[k][x_new(i, j, Length)] != maxData) {
	//				inputImage[k][x_new(i, j, Length)] = maxData;
	//			}
	//		}
	//	}
	//}
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/newSegment.raw");

	//copy input image in container
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				ImageData.imageDataPtr[k][x_new(i, j, Length)] = (dataType)image[k][x_new(i, j, Length)];
			}
		}
	}

	Segmentation_Parameters segment_parameters;
	segment_parameters.maxNoGSIteration = 100; segment_parameters.coef = 1000; segment_parameters.eps2 = 0.1;
	segment_parameters.gauss_seidelTolerance = 1e-6; segment_parameters.h = 1; segment_parameters.numberOfTimeStep = 100;
	segment_parameters.tau = 0.6; segment_parameters.omega_c = 1.0; segment_parameters.mod = 1; segment_parameters.maxNoOfTimeSteps = 100;
	Point3D * center_segment = (Point3D*)malloc(sizeof(Point3D)); center_segment->x = i_max; center_segment->y = j_max; center_segment->z = k_max;
	size_t number_of_centers = 1;
	Filter_Parameters filtering_parameters;
	filtering_parameters.timeStepSize = 0.5; filtering_parameters.edge_detector_coefficient = 1000; filtering_parameters.maxNumberOfSolverIteration = 100;
	filtering_parameters.eps2 = 0.01; filtering_parameters.omega_c = 1.5; filtering_parameters.timeStepsNum = 1; filtering_parameters.tolerance = 1e-4;
	filtering_parameters.h = 1.0; filtering_parameters.p = 1; filtering_parameters.sigma = 1.2; filtering_parameters.coef = 1000;
	unsigned char outputPath[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/";
	//rescaleNewRange(ImageData.imageDataPtr, Length, Width, Height, 0, 1);
	subsurfSegmentation(ImageData, segment_parameters, filtering_parameters, center_segment, number_of_centers, outputPath);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/distanceMap.raw");

	//free memory
	for (k = 0; k < Height; k++) {
		free(image[k]); 
		free(ImageData.imageDataPtr[k]);
		free(distanceMap[k]);
		//free(labelArray[k]); free(status[k]);
	}
	free(image); free(ImageData.imageDataPtr); free(distanceMap);
	//free(labelArray); free(status); free(countingArray);
	//free(Histogram);
	
	//free(HistoSliceBySlice); free(HistoRowByRow); free(HistoColumnByColumn);
	//for (k = 0; k < height_new; k++) {
	//	free(distanceMap[k]); free(croppedImage);
	//}
	//free(croppedImage); free(distanceMap);
	//free(Proba); free(histogram); free(interClassVariance);
	//free(cenTroid);
	
	return EXIT_SUCCESS;
}
