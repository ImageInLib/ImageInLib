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
#include "../src/segmentation3D_subsurf.h"

#define originalMean 71.1245
#define offSet 1024
#define standarDeviation 22.001
#define thresmin 995
#define thresmax 1213
#define minimalSize 2000


int main() {

	size_t i, j, k;
	const size_t Width = 512;
	const size_t Length = 512;
	const size_t Height = 406;
	const size_t dim2D = Width * Length;

	short** image = (short**)malloc(Height * sizeof(short*));
	dataType** imageData = (dataType**)malloc(Height * sizeof(dataType*));
	for (k = 0; k < Height; k++) {
		image[k] = (short*)malloc(dim2D * sizeof(short));
		imageData[k] = (dataType*)malloc(dim2D * sizeof(dataType));
	}
	if (image == NULL || imageData == NULL) return false;

	//initialization
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				image[k][x_new(i, j, Length)] = 0;
				imageData[k][x_new(i, j, Length)] = 0;
			}
		}
	}

	//Loading
	load3dArrayRAW<short>(image, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/input/patient2.raw");

	//copy on data container
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				imageData[k][x_new(i, j, Length)] = (dataType)image[k][x_new(i, j, Length)];
			}
		}
	}
	//store3dRawData<dataType>(imageData, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/loaded.raw");

	////Preparing Image container for float format
	//Image_Data ImageData; ImageData.height = Height; ImageData.length = Length; ImageData.width = Width;
	//ImageData.imageDataPtr = (dataType**)malloc(Height * sizeof(dataType*));
	//for (k = 0;k < Height;k++) {
	//	ImageData.imageDataPtr[k] = (dataType*)malloc(dim2D * sizeof(dataType));
	//}
	//if (ImageData.imageDataPtr == NULL) return false;
	//initialize3dArrayD(ImageData.imageDataPtr, Length, Width, Height, 0);
	////copy input image in container
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			ImageData.imageDataPtr[k][x_new(i, j, Length)] = (dataType)image[k][x_new(i, j, Length)];
	//		}
	//	}
	//}

	//Storage_Flags flags = { false, false };
	//unsigned char name[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/loaded.vtk";
	//store3dDataVtkD(ImageData.imageDataPtr, Length, Width, Height, name, 1, flags);
	//store3dDataArrayD(ImageData.imageDataPtr, Length, Width, Height, name, flags);

	//Find min and max values
	dataType minData = 10000, maxData = -10000;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (imageData[k][x_new(i, j, Length)] < minData)
					minData = imageData[k][x_new(i, j, Length)];
				if (imageData[k][x_new(i, j, Length)] > maxData)
					maxData = imageData[k][x_new(i, j, Length)];
			}
		}
	}
	printf("Min data = %1.lf and Max data = %.1lf\n", minData, maxData);

	//for (k = 0; k < zDim; k++) {
	//	for (i = 0; i < xDim; i++) {
	//		for (j = 0; j < yDim; j++) {
	//			if (imageData[k][x_new(i, j, xDim)] == maxData) {
	//				if (k < kMinCrop) kMinCrop = k; if (k > kMaxCrop) kMaxCrop = k;
	//				if (i < iMinCrop) iMinCrop = i; if (i > iMaxCrop) iMaxCrop = i;
	//				if (j < jMinCrop) jMinCrop = j; if (j > jMaxCrop) jMaxCrop = j;
	//			}
	//		}
	//	}
	//}
	//printf("Coodinates of the cropped volume : %d , %d, %d \n", iMinCrop, jMinCrop, kMinCrop);
	//printf("Coodinates of the cropped volume : %d , %d, %d \n", iMaxCrop, jMaxCrop, kMaxCrop);

	//const size_t Length = iMaxCrop - iMinCrop + 1, Width = jMaxCrop - jMinCrop + 1, Height = kMaxCrop - kMinCrop + 1;
	//printf("Dimensions of cropped volume : %d x %d x %d", Length, Width, Height);
	//dataType** croppedVolume = (dataType**)malloc(Height * sizeof(dataType*));
	//dataType** distanceMap = (dataType**)malloc(Height * sizeof(dataType*));
	//for (k = 0; k < Height; k++) {
	//	croppedVolume[k] = (dataType*)malloc(Length * Width * sizeof(dataType));
	//	distanceMap[k] = (dataType*)malloc(Length * Width * sizeof(dataType));
	//}
	//if (croppedVolume == NULL || distanceMap == NULL) return false;

	//initialize3dArrayD(croppedVolume, Length, Width, Height, 0);
	//initialize3dArrayD(distanceMap, Length, Width, Height, 0);

	//size_t i_c, j_c, k_c;
	//for (k_c = 0, k = kMinCrop; k_c < Height; k_c++, k++) {
	//	for (i_c = 0, i = iMinCrop; i_c < Length; i_c++, i++) {
	//		for (j_c = 0, j = jMinCrop; j_c < Width; j_c++, j++) {
	//			croppedVolume[k_c][x_new(i_c, j_c, Length)] = (dataType)image[k][x_new(i, j, xDim)];
	//		}
	//	}
	//}
	//printf("\nDimensions of the cropped volume : %d x %d x %d \n", Width, Length, Height);
	//store3dRawData<dataType>(croppedVolume, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/newVolume.raw");

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
	//thresholding3dFunctionN(imageData, Length, Width, Height, thresmin, thresmax, minData, maxData);

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

	//// Centroid of the segmented Liver
	//dataType* cenTroid = (dataType*)malloc(3 * sizeof(dataType));
	//centroidImage(ImageData.imageDataPtr, cenTroid, Height, Length, Width, 0);
	//printf("Coordinates of centroid  : ");
	//for (i = 0; i < 3; i++) {
	//	printf("%f, ", cenTroid[i]);
	//}
	//printf("\n");

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

	//char pathSaveDil[300];
	//for (int i = 0; i < 7; i++) {
	//	dilatation3dHeighteenNeigbours(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//	sprintf_s(pathSaveDil, "C:/Users/Konan Allaly/Documents/Tests/output/dilatation0%d.raw", i);
	//	store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, pathSaveDil);
	//}


	//Manual Cropping
	size_t iMinCrop = 170, iMaxCrop = 360, jMinCrop = 130, jMaxCrop = 310, kMinCrop = 157, kMaxCrop = 231;
	size_t lengthNew = jMaxCrop - jMinCrop + 1, widthNew = iMaxCrop - iMinCrop + 1, heightNew = kMaxCrop - kMinCrop + 1;

	Image_Data segment;
	segment.height = heightNew; segment.length = lengthNew; segment.width = widthNew;
	segment.imageDataPtr = (dataType**)malloc(heightNew * sizeof(dataType*));
	dataType** distanceMap = (dataType**)malloc(heightNew * sizeof(dataType*));
	dataType** maskThresh = (dataType**)malloc(heightNew * sizeof(dataType*));
	for (k = 0; k < heightNew; k++) {
		segment.imageDataPtr[k] = (dataType*)malloc(lengthNew * widthNew * sizeof(dataType));
		distanceMap[k] = (dataType*)malloc(lengthNew * widthNew * sizeof(dataType));
		maskThresh[k] = (dataType*)malloc(lengthNew * widthNew * sizeof(dataType));
	}
	if (segment.imageDataPtr == NULL || distanceMap == NULL || maskThresh == NULL) return false;
	initialize3dArrayD(segment.imageDataPtr, lengthNew, widthNew, heightNew, 0);
	initialize3dArrayD(distanceMap, lengthNew, widthNew, heightNew, 0);
	initialize3dArrayD(maskThresh, lengthNew, widthNew, heightNew, 0);
	//copy data (Cropping)
	for (k = 0; k < heightNew; k++) {
		for (i = 0; i < lengthNew; i++) {
			for (j = 0; j < widthNew; j++) {
				segment.imageDataPtr[k][x_new(i, j, lengthNew)] = imageData[k + kMinCrop][x_new(i + iMinCrop, j + jMinCrop, Length)];
				maskThresh[k][x_new(i, j, lengthNew)] = imageData[k + kMinCrop][x_new(i + iMinCrop, j + jMinCrop, Length)];
			}
		}
	}
	printf("New dimensions %d x %d x %d\n", widthNew, lengthNew, heightNew);
	store3dRawData<dataType>(segment.imageDataPtr, lengthNew, widthNew, heightNew, "C:/Users/Konan Allaly/Documents/Tests/output/newVolume.raw");

	//Filtering
	dataType timeStepSize = 1.2, h = 1.0, sigma = 1e-3, K = 200, omega_c = 1.1, tol = 1e-7, coef = 1e-3, eps2 = 1e-3;
	size_t p = 1, timeStepsNum = 1, iter_max = 1000;
	Filter_Parameters PM_filterParameters;
	const FilterMethod methodFiltering = GEODESIC_MEAN_CURVATURE_FILTER;
	PM_filterParameters = {timeStepSize, h, sigma, K, omega_c, tol, eps2, coef, p, timeStepsNum, iter_max};
	rescaleNewRange(segment.imageDataPtr, lengthNew, widthNew, heightNew, 0, 1);
	filterImage(segment, PM_filterParameters, methodFiltering);
	//rescaleNewRange(ImageData.imageDataPtr, Length, Width, Height, minData, maxData);
	store3dRawData<dataType>(segment.imageDataPtr, lengthNew, widthNew, heightNew, "C:/Users/Konan Allaly/Documents/Tests/output/filtered_3GMC.raw");

	thresholding3dFunctionN(maskThresh, lengthNew, widthNew, heightNew, thresmin, thresmax, minData, maxData);
	store3dRawData<dataType>(maskThresh, lengthNew, widthNew, heightNew, "C:/Users/Konan Allaly/Documents/Tests/output/Thresh.raw");

	//Fast sweeping for cropped volume
	fastSweepingFunction_3D(distanceMap, maskThresh, lengthNew, widthNew, heightNew, 1, 100000000, minData);
	store3dRawData<dataType>(distanceMap, lengthNew, widthNew, heightNew, "C:/Users/Konan Allaly/Documents/Tests/output/distanceMap.raw");
	dataType distanceMax = -1;
	for (k = 0; k < heightNew; k++) {
		for (i = 0; i < lengthNew; i++) {
			for (j = 0; j < widthNew; j++) {
				if (distanceMap[k][x_new(i, j, lengthNew)] >= distanceMax) {
					distanceMax = distanceMap[k][x_new(i, j, lengthNew)];
				}
			}
		}
	}
	dataType cptMax = 0;
	int i_max = 0, j_max = 0, k_max = 0;
	for (k = 0; k < heightNew; k++) {
		for (i = 0; i < lengthNew; i++) {
			for (j = 0; j < widthNew; j++) {
				if (distanceMap[k][x_new(i, j, lengthNew)] == distanceMax) {
					cptMax++;
					i_max = (int)i; j_max = (int)j; k_max = (int)k;
				}
			}
		}
	}
	printf("\nThe maximal distance is = %f, and %f points have that distance \n", distanceMax, cptMax);
	printf("Coordinates of the point with highest distance : (%d, %d, %d)\n\n", j_max, i_max, k_max);

	//Segmentation
	Segmentation_Parameters segment_parameters;
	segment_parameters.maxNoGSIteration = 100; segment_parameters.coef = 200; segment_parameters.eps2 = 1e-10;
	segment_parameters.numberOfTimeStep = 300; segment_parameters.mod = 1; segment_parameters.maxNoOfTimeSteps = 300;
	segment_parameters.tau = 4.0; segment_parameters.h = 1; segment_parameters.omega_c = 1.1;
	segment_parameters.gauss_seidelTolerance = 1e-6; segment_parameters.segTolerance = 1e-4;
	Point3D* center_segment = (Point3D*)malloc(sizeof(Point3D)); //center_segment->x = 50; center_segment->y = 55; center_segment->z = 20;
	center_segment->x = i_max; center_segment->y = j_max; center_segment->z = k_max;
	size_t number_of_centers = 1;
	Filter_Parameters filtering_parameters;
	filtering_parameters.timeStepSize = 1.0/6.0; filtering_parameters.edge_detector_coefficient = 200; filtering_parameters.maxNumberOfSolverIteration = 300;
	filtering_parameters.eps2 = 0.001; filtering_parameters.omega_c = 1.2; filtering_parameters.timeStepsNum = 1; filtering_parameters.tolerance = 1e-4;
	filtering_parameters.h = 1.0; filtering_parameters.p = 1; filtering_parameters.sigma = 1.2; filtering_parameters.coef = 1000;
	unsigned char outputPath[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/";
	//rescaleNewRange(segment.imageDataPtr, lengthNew, widthNew, heightNew, 0, 1);
	////store3dRawData<dataType>(segment.imageDataPtr, xDim, yDim, zDim, "C:/Users/Konan Allaly/Documents/Tests/output/rescalled.raw");
	subsurfSegmentation(segment, segment_parameters, filtering_parameters, center_segment, number_of_centers, outputPath);

	////Fast sweeping for whole volume
	//dataType** distanceMap = (dataType**)malloc(Height * sizeof(dataType*));
	//dataType** ballCode = (dataType**)malloc(Height * sizeof(dataType*));
	//dataType** ballSlicer = (dataType**)malloc(Height * sizeof(dataType*));
	//for (k = 0; k < Height; k++) {
	//	distanceMap[k] = (dataType*)malloc(dim2D * sizeof(dataType));
	//	ballCode[k] = (dataType*)malloc(dim2D * sizeof(dataType));
	//	ballSlicer[k] = (dataType*)malloc(dim2D * sizeof(dataType));
	//}
	//if (distanceMap == NULL || ballCode == NULL || ballSlicer == NULL) return false;
	////initialization
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			distanceMap[k][x_new(i, j, Length)] = 0;
	//			ballCode[k][x_new(i, j, Length)] = 0;
	//			ballSlicer[k][x_new(i, j, Length)] = 0;
	//		}
	//	}
	//}
	////initialize3dArrayD(distanceMap, Height, Width, Length, 0);
	////initialize3dArrayD(ballCode, Height, Width, Length, 0);
	////initialize3dArrayD(ballSlicer, Height, Width, Length, 0);
	//fastSweepingFunction_3D(distanceMap, imageData, Length, Width, Height, 1, 100000000, minData);
	//store3dRawData<dataType>(distanceMap, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/distanceMapP7.raw");
	//dataType distanceMax = -1;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (distanceMap[k][x_new(i, j, Length)] >= distanceMax) {
	//				distanceMax = distanceMap[k][x_new(i, j, Length)];
	//			}
	//		}
	//	}
	//}
	//dataType cptMax = 0;
	//int i_max = 0, j_max = 0, k_max = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (distanceMap[k][x_new(i, j, Length)] == distanceMax) {
	//				cptMax++;
	//				i_max = (int)i; j_max = (int)j; k_max = (int)k;
	//			}
	//		}
	//	}
	//}
	//printf("\nThe maximal distance is = %f, and %f points have that distance \n", distanceMax, cptMax);
	//printf("Coordinates of the point with highest distance : (%d, %d, %d)\n", j_max, i_max, k_max);
	//float i_s = 189.022, j_s = 308.971, k_s = 351.193;
	//printf("Centroid from Slicer : (%.3f, %.3f, %.3f)\n", i_s, j_s, k_s);
	//printf("The distance between the two points is : %.3f", sqrt((j_max - i_s)*((j_max - i_s)) + (i_max - j_s)*(i_max - j_s) + (k_max - k_s)* (k_max - k_s)) );
	////Ball arround the highest distance point and centroid from Slicer
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			//ball around the highest distance
	//			if (sqrt( (i_max - i)* (i_max - i) + (j_max - j) * (j_max - j) + (k_max - k) * (k_max - k) ) < 10 ) {
	//				ballCode[k][x_new(i, j, Length)] = maxData;
	//			}
	//			else {
	//				ballCode[k][x_new(i, j, Length)] = minData;
	//			}
	//			//ball around centroid
	//			if (sqrt((i_s - j) * (i_s - j) + (j_s - i) * (j_s - i) + (k_s - k) * (k_s - k)) < 10) {
	//				ballSlicer[k][x_new(i, j, Length)] = maxData;
	//			}
	//			else {
	//				ballSlicer[k][x_new(i, j, Length)] = minData;
	//			}
	//		}
	//	}
	//}
	//store3dRawData<dataType>(ballCode, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/ballCodeP7.raw");
	//store3dRawData<dataType>(ballSlicer, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/ballSlicerP7.raw");

	//free memory
	for (k = 0; k < Height; k++) {
		free(image[k]);
		free(imageData[k]);
		//free(distanceMap[k]); free(ballCode[k]); free(ballSlicer[k]);
		//free(labelArray[k]); free(status[k]);
	}
	free(image); free(imageData); 
	//free(distanceMap); free(ballCode); free(ballSlicer);
	//free(labelArray); free(status); free(countingArray);

	for (k = 0; k < heightNew; k++) {
		free(segment.imageDataPtr[k]);
		free(distanceMap[k]);
		free(maskThresh[k]);
	}
	free(segment.imageDataPtr); 
	free(distanceMap); 
	free(maskThresh);
	free(center_segment);

	return EXIT_SUCCESS;
}
