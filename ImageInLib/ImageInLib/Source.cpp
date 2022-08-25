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



#define object 0
#define background SHRT_MAX

#define originalMean 71
#define offSet 1024
#define margin 100
#define thresmin (originalMean + offSet - margin)
#define thresmax (originalMean + offSet + margin)
#define minimalSize 1


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
	//			if (sqrt((90 - i) * (90 - i) + (90 - j) * (90 - j) + (200 - k) * (200 - k)) < 90) {
	//				imageData[k][x_new(i, j, Length)] = object;
	//			}
	//			if (sqrt((410 - i) * (410 - i) + (90 - j) * (90 - j) + (200 - k) * (200 - k)) < 90) {
	//				imageData[k][x_new(i, j, Length)] = object;
	//			}
	//			if (sqrt((90 - i) * (90 - i) + (410 - j) * (410 - j) + (200 - k) * (200 - k)) < 90) {
	//				imageData[k][x_new(i, j, Length)] = object;
	//			}
	//			if (sqrt((410 - i) * (410 - i) + (410 - j) * (410 - j) + (200 - k) * (200 - k)) < 90) {
	//				imageData[k][x_new(i, j, Length)] = object;
	//			}
	//			if (sqrt((250 - i) * (250 - i) + (250 - j) * (250 - j) + (200 - k) * (200 - k)) < 90) {
	//				imageData[k][x_new(i, j, Length)] = object;
	//			}
	//		}
	//	}
	//}
	//char pathSave[100];
	//store3dRawData<dataType>(image, Length, Width, Height, "output/artificialImage.raw");

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

	// Preparing Image container for float format
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

	//find min and max values of data
	float minData = 10000, maxData = 0;
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

	//Rescalling
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				ImageData.imageDataPtr[k][x_new(i, j, Length)] = ImageData.imageDataPtr[k][x_new(i, j, Length)] / maxData;
			}
		}
	}

	//Filtering
	const size_t maxNumberOfSolverIteration = 1000;
	float  coef = 0.01;
	float  eps2 = 1e-4;
	size_t numberOfTimeStep = 1000;
	Filter_Parameters MC_filterParameters;
	MC_filterParameters = { 1.2, 10, 0.1, 1000, 1.5, 1e-7, 1e-4, 1, 100, 1000, 1000};
	const FilterMethod methodFiltering = MEAN_CURVATURE_FILTER;
	filterImage(ImageData, MC_filterParameters, maxNumberOfSolverIteration, coef, eps2, numberOfTimeStep, methodFiltering);

	//Rescalling
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				ImageData.imageDataPtr[k][x_new(i, j, Length)] = ImageData.imageDataPtr[k][x_new(i, j, Length)] * maxData + 0.5;
			}
		}
	}
	store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "output/filteredImage.raw");

	//Thresholding
	thresholding3dFunctionN(ImageData.imageDataPtr, Length, Width, Height, thresmin, thresmax, object, background);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "output/thres.raw");

	erosion3D(ImageData.imageDataPtr, Length, Width, Height, object, background);
	erosion3D(ImageData.imageDataPtr, Length, Width, Height, object, background);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "output/TwoErosionWF.raw");

	//dilatation3D(ImageData.imageDataPtr, Length, Width, Height, object, background);
	//dilatation3D(ImageData.imageDataPtr, Length, Width, Height, object, background);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "output/DilatationAfterTwoErosions.raw");

	double start = clock();
	labelling3D(ImageData.imageDataPtr, labelArray, status, Length, Width, Height, object);
	double finish = clock();
	printf("Execution time for the new code : %.3lf \n", (finish - start) / CLOCKS_PER_SEC);

	//number of region cells
	int numberOfRegionsCells = 0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (ImageData.imageDataPtr[k][x_new(i, j, Length)] == object) {
					numberOfRegionsCells++;
				}
			}
		}
	}
	printf("Number of Regions Cells : %d \n", numberOfRegionsCells);

	//Counting
	int* countingArray = (int*)malloc(numberOfRegionsCells * sizeof(int));
	if (countingArray == NULL) return false;
	for (i = 0; i < numberOfRegionsCells; i++) countingArray[i] = 0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (labelArray[k][x_new(i, j, Length)] > 0) {
					countingArray[labelArray[k][x_new(i, j, Length)]]++;
				}
			}
		}
	}

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

	//Number of regions
	int numberOfRegions = 0;
	for (i = 0; i < numberOfRegionsCells; i++) {
		if (countingArray[i] > 0) {
			numberOfRegions++;
		}
	}
	printf("Number of regions = %d \n", numberOfRegions);

	//Statistics
	FILE* file;
	if (fopen_s(&file, "output/Statistics.txt", "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(file, " 3D real image with %d Slices \n", Height);
	fprintf(file, " Labelling function takes = %.5lf seconds\n", (finish - start) / CLOCKS_PER_SEC);
	fprintf(file, " Number of Objects Cells = %d \n", numberOfRegionsCells);
	fprintf(file, " Number of Regions  = %d \n", numberOfRegions);
	fprintf(file, "#############################################################\n\n");
	fprintf(file, " Regions with more than = %d voxels\n", minimalSize);
	fprintf(file, "\n\n");
	fprintf(file, "Region Label        Count");
	for (i = 0; i < numberOfRegionsCells; i++) {
		if (countingArray[i] >= minimalSize) {
			fprintf(file, "\n");
			fprintf(file, "   %d                 %d", i, countingArray[i]);
		}
	}
	fclose(file);

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
	store3dRawData<int>(labelArray, Length, Width, Height, "output/segmentedN.raw");

	////Keep one region 
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (labelArray[k][x_new(i, j, Length)] != 383) {
	//				labelArray[k][x_new(i, j, Length)] = 0;
	//			}
	//		}
	//	}
	//}
	//store3dRawData<int>(labelArray, Length, Width, Height, "output/region383.raw");

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
	}
	free(imageData);
	free(labelArray); free(image); free(status); //free(countingArray);
	return EXIT_SUCCESS;
}
