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


#define object 0
#define background SHRT_MAX

#define originalMean 71
#define offSet 1024
#define margin 40
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
	dataType** image = (dataType**)malloc(Height * sizeof(dataType*));
	int** labelArray = (int**)malloc(Height * sizeof(int*));
	bool** status = (bool**)malloc(Height * sizeof(bool*));
	for (k = 0; k < Height; k++) {
		imageData[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		image[k] = (dataType*)malloc(dim2D * sizeof(dataType));
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
	unsigned char pathPtr[] = "data/3D/patient1b.raw";
	OperationType operation = LOAD_DATA_RAW;
	LoadDataType dType = BINARY_DATA;
	Storage_Flags flags = { false, false };
	manageFile(image, Length, Width, Height, pathPtr, operation, dType, flags);

	//Copy of input image in container
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				imageData[k][x_new(i, j, Length)] = image[k][x_new(i, j, Length)];
			}
		}
	}

	//Thresholding and saving of thresholding image
	thresholding3dFunctionN(imageData, Length, Width, Height, thresmin, thresmax, object, background);
	store3dRawData<dataType>(imageData, Length, Width, Height, "output/thresbeforeErosion.raw");

	erosion3D(imageData, Length, Width, Height, object, background);

	store3dRawData<dataType>(imageData, Length, Width, Height, "output/thresAfterErosionFirstTime.raw");

	erosion3D(imageData, Length, Width, Height, object, background);

	store3dRawData<dataType>(imageData, Length, Width, Height, "output/thresAfterErosionSecondTime.raw");

	double start = clock();
	labelling3D(imageData, labelArray, status, Length, Width, Height, object);
	double finish = clock();
	printf("Execution time for the new code : %.3lf \n", (finish - start) / CLOCKS_PER_SEC);

	//number of region cells
	int numberOfRegionsCells = 0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (imageData[k][x_new(i, j, Length)] == object) {
					numberOfRegionsCells++;
				}
			}
		}
	}
	printf("Number of Regions Cells : %d \n", numberOfRegionsCells);

	int* countingArray = (int*)malloc(numberOfRegionsCells * sizeof(int));
	if (countingArray == NULL) return false;
	for (i = 0; i < numberOfRegionsCells; i++) countingArray[i] = 0;

	//Counting
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (labelArray[k][x_new(i, j, Length)] > 0) {
					countingArray[labelArray[k][x_new(i, j, Length)]]++;
				}
			}
		}
	}

	//Remove small regions 
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (countingArray[labelArray[k][x_new(i, j, Length)]] < minimalSize) {
					countingArray[labelArray[k][x_new(i, j, Length)]] = 0;
				}
			}
		}
	}

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
	if (fopen_s(&file, "output/Labels_Count_New.txt", "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(file, " 3D real image with %d Slices \n", Height);
	fprintf(file, " New labelling function takes = %.5lf seconds\n", (finish - start) / CLOCKS_PER_SEC);
	fprintf(file, " Number of Objects Cells = %d \n", numberOfRegionsCells);
	fprintf(file, " Number of Regions with new code = %d \n", numberOfRegions);
	fprintf(file, "#############################################################\n");
	fprintf(file, "\n\n");
	fprintf(file, "Region Label        Count");
	for (i = 0; i < numberOfRegionsCells; i++) {
		if (countingArray[i] >= minimalSize) {
			fprintf(file, "\n");
			fprintf(file, "   %d                 %d", i, countingArray[i]);
		}
	}
	fclose(file);


	//change the backrgound
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (countingArray[labelArray[k][x_new(i, j, Length)]] < minimalSize) {
					labelArray[k][x_new(i, j, Length)] = 0;
				}
			}
		}
	}

	//Saving with template function
	store3dRawData<int>(labelArray, Length, Width, Height, "output/segment.raw");

	//free memory
	for (k = 0; k < Height; k++) {
		free(image[k]);
		free(labelArray[k]); free(imageData[k]); free(status[k]);
	}
	free(imageData);
	free(labelArray); free(image); free(status); free(countingArray);
	return EXIT_SUCCESS;
}
