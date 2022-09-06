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
#include "../src/data_initialization.h"
#include "noiseGeneration.h"
#include "filtering.h"
#include "../src/data_load.h"
#include "../src/data_storage.h"
#include "../src/thresholding.h"
#include "../src/generate_3d_shapes.h"
#include "Labelling.h"
#include "template_functions.h"
#include "../src/endianity_bl.h"
#include "../include/morphological_change.h"
#include "../src/shape_registration.h"



#define object 0
#define background SHRT_MAX

#define originalMean 71
#define offSet 1024
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

	dataType** imageData = (dataType**)malloc(Height * sizeof(dataType*));
	short** image = (short**)malloc(Height * sizeof(short*));
	int** labelArray = (int**)malloc(Height * sizeof(int*));
	bool** status = (bool**)malloc(Height * sizeof(bool*));
	for (k = 0; k < Height; k++) {
		imageData[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		image[k] = (short*)malloc(dim2D * sizeof(short));
		labelArray[k] = (int*)malloc(dim2D * sizeof(int));
		status[k] = (bool*)malloc(dim2D * sizeof(bool));
	}
	if (imageData == NULL) return false;
	if (image == NULL) return false;
	if (labelArray == NULL) return false;
	if (status == NULL) return false;

	//initialization
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				imageData[k][x_new(i, j, Length)] = 0;
				image[k][x_new(i, j, Length)] = 0;
				labelArray[k][x_new(i, j, Length)] = 0;
				status[k][x_new(i, j, Length)] = false;
			}
		}
	}

	//Loading
	//unsigned char pathPtr[] = "data/3D/patient1b.raw";
	//OperationType operation = LOAD_DATA_RAW;
	//LoadDataType dType = BINARY_DATA;
	//Storage_Flags flags = { false, false };
	//manageFile(image, Length, Width, Height, pathPtr, operation, dType, flags);
	load3dArrayRAW<short>(image, Length, Width, Height, "data/3D/patient1b.raw");
	//store3dRawData<short>(image, Length, Width, Height, "output/loadedImage.raw");

	////Create artificial image
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (sqrt((128 - i) * (128 - i) + (384 - j) * (384 - j) + (254 - k) * (254 - k)) < 50) {
	//				imageData[k][x_new(i, j, Length)] = 1;
	//			}
	//		}
	//	}
	//}
	//char pathSave[100];
	//store3dRawData<dataType>(imageData, Length, Width, Height, "output/ballInRightCorner.raw");

	////loop for erosion
	//for (k = 1; k < 28; k++) {
	//	erosion3D(imageData, Length, Width, Height, object, background);
	//	sprintf_s(pathSave, "output/Erosion/erosion%d.raw", k);
	//	store3dRawData<float>(imageData, Length, Width, Height, pathSave);
	//}
	////Loop for dilatation
	//for (k = 1; k <= 28; k++) {
	//	dilatation3D(imageData, Length, Width, Height, object, background);
	//	sprintf_s(pathSave, "output/Dilatation/dilatation%d.raw", k);
	//	store3dRawData<float>(imageData, Length, Width, Height, pathSave);
	//}

	////Copy of input image in container
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			imageData[k][x_new(i, j, Length)] = image[k][x_new(i, j, Length)];
	//		}
	//	}
	//}

	//Preparing Image container for float format
	Image_Data ImageData;
	ImageData.height = Height;
	ImageData.length = Length;
	ImageData.width = Width;
	ImageData.imageDataPtr = (dataType**)malloc(Height * sizeof(dataType*));
	for (k = 0;k < Height;k++) {
		ImageData.imageDataPtr[k] = (dataType*)malloc(dim2D * sizeof(dataType));
	}
	if (ImageData.imageDataPtr == NULL) return false;

	//copy input image in container
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				ImageData.imageDataPtr[k][x_new(i, j, Length)] = (float)image[k][x_new(i, j, Length)];
			}
		}
	}

	//Find min and max values
	short minData = 10000, maxData = -10000;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (image[k][x_new(i, j, Length)] < minData)
					minData = image[k][x_new(i, j, Length)];
				if (image[k][x_new(i, j, Length)] > maxData)
					maxData = image[k][x_new(i, j, Length)];
			}
		}
	}

	////Compute and save histogram
	//short totalClass = maxData - minData + 1;
	//dataType* histogram = (dataType*)malloc(totalClass * sizeof(dataType));
	//for (i = 0; i < totalClass; i++) histogram[i] = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			histogram[image[k][x_new(i, j, Length)]]++;
	//		}
	//	}
	//}
	//FILE* histo;
	//if (fopen_s(&histo, "output/histogram.csv", "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//for (i = 0; i < totalClass; i++) {
	//	fprintf(histo, "%f \n", histogram[i]);
	//}
	//fclose(histo);
	//float numberOfCells = 0;
	//for (i = 0; i < totalClass; i++) {
	//	numberOfCells = numberOfCells + histogram[i];
	//}
	//float q1 = 0, mu1 = 0, sigma1 = 0, q2 = 0, mu2 = 0, sigma2 = 0, sigmab = 0;
	//short T = 2000;
	//float* P = (float*)malloc(totalClass * sizeof(float));
	//for (i = 0; i < totalClass; i++) P[i] = 0;
	//for (i = 0; i <= T; i++) {
	//	P[i] = histogram[i] / numberOfCells;
	//	q1 = q1 + P[i];
	//}
	//for (i = T+1 ; i < totalClass ; i++) {
	//	P[i] = histogram[i] / numberOfCells;
	//	q2 = q2 + P[i];
	//}
	//for (i = 0; i <= T; i++) {
	//	mu1 = mu1 + i*P[i];
	//}
	//mu1 = mu1 / q1;
	//for (i = T + 1; i < totalClass; i++) {
	//	mu2 = mu2 + i * P[i];
	//}
	//mu2 = mu2 / q2;
	//for (i = 0; i <= T; i++) {
	//	sigma1 = sigma1 + (i - mu1)*(i - mu1)*P[i];
	//}
	//sigma1 = sigma1 / q1;
	//for (i = T + 1; i < totalClass; i++) {
	//	sigma2 = sigma2 + (i - mu2) * (i - mu2) * P[i];
	//}
	//sigma2 = sigma2 / q2;
	//sigmab = q1 * sigma1 + q2 * sigma1;

	//rescaleNewRange(ImageData.imageDataPtr, Length, Width, Height, 0, 1);

	NoiseParameters Nparameters = { 0.05, 0.1, 10, minData, maxData };
	const NoiseType Ntype = SALT_AND_PEPPER;
	addNoiseToImage(ImageData.imageDataPtr, Length, Width, Height, Nparameters, Ntype);
	store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "output/noisyImage.raw");

	////Filtering
	//const size_t maxNumberOfSolverIteration = 100;
	//float  coef = 0.01;
	//float  eps2 = 1e-2;
	//size_t numberOfTimeStep = 10;
	//Filter_Parameters MC_filterParameters;
	//MC_filterParameters = {5, 1, 0.1, 0.018, 1.5, 5e-4, 1e-2, 1, 10, 100, 10};
	//const FilterMethod methodFiltering = MEAN_CURVATURE_FILTER;
	////Time step = 10
	//filterImage(ImageData, MC_filterParameters, maxNumberOfSolverIteration, coef, eps2, numberOfTimeStep, methodFiltering);
	//rescaleNewRange(ImageData.imageDataPtr, Length, Width, Height, minData, maxData);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "output/filteredImage10timeStep.raw");

	//Thresholding
	//thresholding3dFunctionN(ImageData.imageDataPtr, Length, Width, Height, 1823, 1944, minData, maxData);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "output/thres.raw");

	////find the Centroid
	//dataType* cenTroid = (dataType*)malloc(3 * sizeof(dataType));
	//centroidImage(imageData, cenTroid, Height, Length, Width, 0);
	//printf("Coordinates of centroid  : ");
	//for (i = 0; i < 3; i++) {
	//	printf("%f, ", cenTroid[i]);
	//}
	//printf("\n");

	//erosion3D(ImageData.imageDataPtr, Length, Width, Height, minData, maxData);
	//erosion3D(ImageData.imageDataPtr, Length, Width, Height, minData, maxData);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "output/Erosion.raw");

	//dilatation3D(ImageData.imageDataPtr, Length, Width, Height, object, background);
	//dilatation3D(ImageData.imageDataPtr, Length, Width, Height, object, background);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "output/DilatationAfterErosion.raw");

	//double start = clock();
	//labelling3D(ImageData.imageDataPtr, labelArray, status, Length, Width, Height, minData);
	//double finish = clock();
	//printf("Execution time for the new code : %.3lf \n", (finish - start) / CLOCKS_PER_SEC);

	////number of region cells
	//int numberOfRegionsCells = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (ImageData.imageDataPtr[k][x_new(i, j, Length)] == minData) {
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
	//if (fopen_s(&file, "output/Statistics.txt", "w") != 0) {
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
	//				//imageData[k][x_new(i, j, Length)] = background;
	//			}
	//		}
	//	}
	//}

	//Saving with template function
	//store3dRawData<int>(labelArray, Length, Width, Height, "output/segmented.raw");

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
		free(labelArray[k]); free(imageData[k]); free(status[k]);
		//free(ImageData.imageDataPtr[k]);
	}
	free(imageData);
	free(labelArray); free(image); free(status); //free(countingArray);
	//free(ImageData.imageDataPtr); //
	//free(cenTroid);
	//free(P); free(histogram);

	return EXIT_SUCCESS;
}
