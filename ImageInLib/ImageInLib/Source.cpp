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
//#include "../src/thresholding.h"
#include "../src/generate_3d_shapes.h"
#include "Labelling.h"
#include "template_functions.h"


#define object 0
#define background SHRT_MAX

#define thresmin 930
#define thresmax 1150


int main() {

	size_t i, j, k, xd;
	const size_t Length = 512;
	const size_t Width = 512;
	const size_t Height = 50;
	const size_t dim2D = Length * Width;
	const size_t dim3D = Height * Length * Width;

	dataType** image = (dataType**)malloc(Height * sizeof(dataType*));
	for (k = 0; k < Height; k++) {
		image[k] = (dataType*)malloc(dim2D * sizeof(dataType));
	}
	if (image == NULL)
		return false;

	//initialization
	initialize3dArrayD(image, Length, Width, Height, 0);

	//Loading
	//unsigned char pathPtr[] = "data/3D/artificial2.raw";
	unsigned char pathPtr[] = "data/3D/patient2_50Slices.raw";
	OperationType operation = LOAD_DATA_RAW;
	LoadDataType dType = BINARY_DATA;
	Storage_Flags flags = { false, false };
	manageFile(image, Length, Width, Height, pathPtr, operation, dType, flags);

	//Saving of loding image
	operation = STORE_DATA_RAW;
	//unsigned char loadPathPtr[] = "output/loaded.raw";
	//manageFile(image, Length, Width, Height, loadPathPtr, operation, dType, flags);

	//Thresholding and saving of thresholding image
	//thresholding3dFunctionN(image, Length, Width, Height, thresmin, thresmax, background, object);  //930-1150 / 650-1100 // 950-1100
	//unsigned char threshPathPtr[] = "output/thresh.raw";
	//manageFile(image, Length, Width, Height, threshPathPtr, operation, dType, flags);



	//######################################################################################
	//For new code
	dataType** segmented = (dataType**)malloc(Height * sizeof(dataType*));
	for (k = 0; k < Height; k++) {
		segmented[k] = (dataType*)malloc(dim2D * sizeof(dataType));
	}
	if (segmented == NULL)
		return false;

	dataType* imageNew = (dataType*)malloc(dim3D * sizeof(dataType));
	dataType* labelArrayNew = (dataType*)malloc(dim3D * sizeof(dataType));
	bool* status = (bool*)malloc(dim3D * sizeof(bool));

	if (imageNew == NULL)
		return false;
	if (labelArrayNew == NULL)
		return false;
	if (status == NULL)
		return false;

	//Initialization
	initialize3dArrayD(segmented, Length, Width, Height, 0);
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				imageNew[x_flat(i, j, k, Length, Width)] = 0;
				labelArrayNew[x_flat(i, j, k, Length, Width)] = 0;
				status[x_flat(i, j, k, Length, Width)] = false;
			}
		}
	}
	//Copy of input image in container
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				imageNew[x_flat(i, j, k, Length, Width)] = image[k][x_new(i, j, Length)];
			}
		}
	}

	int height = (int)Height;
	int length = (int)Length;
	int width = (int)Width;

	double startNew = clock();
	labelling3D(imageNew, labelArrayNew, status, length, width, height, object);
	double finishNew = clock();
	printf("Execution time for the new code : %.3lf \n", (finishNew - startNew) / CLOCKS_PER_SEC);

	int numberOfRegionsCellsNew = 0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (imageNew[x_flat(i, j, k, Length, Width)] == object) {
					numberOfRegionsCellsNew++;
				}
			}
		}
	}
	printf("Number of Regions Cells : %d\n", numberOfRegionsCellsNew);

	int* countingArrayNew = (int*)malloc(numberOfRegionsCellsNew * sizeof(int));
	if (countingArrayNew == NULL)
		return false;
	for (i = 0; i < numberOfRegionsCellsNew; i++)
		countingArrayNew[i] = 0;

	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (labelArrayNew[x_flat(i, j, k, Length, Width)] > 0) {
					countingArrayNew[labelArrayNew[x_flat(i, j, k, Length, Width)]]++;
				}
			}
		}
	}
	int numberOfRegionsNew = 0;
	for (i = 0; i < numberOfRegionsCellsNew; i++) {
		if (countingArrayNew[i] > 0) {
			numberOfRegionsNew++;
		}
	}
	printf("Number of regions = %d\n", numberOfRegionsNew);


	FILE* file;
	if (fopen_s(&file, "output/Labels_Count_New.txt", "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(file, " 3D real image with %d Slices \n", Height);
	fprintf(file, " Labelling function takes = %.5lf seconds\n", (finishNew - startNew) / CLOCKS_PER_SEC);
	fprintf(file, " Number of Objects Cells = %d \n", numberOfRegionsCellsNew);
	fprintf(file, " Number of Regions = %d \n", numberOfRegionsNew);
	fprintf(file, "#############################################################\n");
	fprintf(file, "\n\n");
	fprintf(file, "Region Label        Count");
	for (i = 0; i < numberOfRegionsCellsNew; i++) {
		if (countingArrayNew[i] > 0) {
			fprintf(file, "\n");
			fprintf(file, "   %d                 %d", i, countingArrayNew[i]);
		}
	}
	fclose(file);

	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (labelArrayNew[x_flat(i, j, k, Length, Width)] == 0) {
					segmented[k][x_new(i, j, Length)] = background;
				}
				else {
					segmented[k][x_new(i, j, Length)] = labelArrayNew[x_flat(i, j, k, Length, Width)];
				}
			}
		}
	}

	//unsigned char segmentedPathPtr[] = "output/segNew.raw";
	//manageFile(segmented, Length, Width, Height, segmentedPathPtr, operation, dType, flags);
	store3dRawData<dataType>(labelArrayNew, Length, Width, Height, "output/segmentN.raw");

	printf("\n##############################################\n\n");

	//###################################################################################
	//For Old code
	dataType** imageOld = (dataType**)malloc(Height * sizeof(dataType*));
	dataType** labelArrayOld = (dataType**)malloc(Height * sizeof(dataType*));
	for (k = 0; k < Height; k++) {
		imageOld[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		labelArrayOld[k] = (dataType*)malloc(dim2D * sizeof(dataType));
	}

	if (imageOld == NULL)
		return false;
	if (labelArrayOld == NULL)
		return false;

	//Initialization
	initialize3dArrayD(imageOld, Length, Width, Height, 0);
	initialize3dArrayD(labelArrayOld, Length, Width, Height, 0);

	//Copy of input image in container
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				imageOld[k][x_new(i, j, Length)] = image[k][x_new(i, j, Length)];
			}
		}
	}

	double startOld = clock();
	regionLabelling3D(imageOld, labelArrayOld, Length, Width, Height, background, object, false, 16);
	double finishOld = clock();
	printf("Execution time for the old code = %.3lf \n", (finishOld - startOld) / CLOCKS_PER_SEC);

	int numberOfRegionsCellsOld = 0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (image[k][x_new(i, j, Length)] == object) {
					numberOfRegionsCellsOld++;
				}
			}
		}
	}
	printf("Number of Regions Cells Old: %d\n", numberOfRegionsCellsOld);

	int* countingArrayOld = (int*)malloc(numberOfRegionsCellsOld * sizeof(int));
	if (countingArrayOld == NULL)
		return false;
	for (i = 0; i < numberOfRegionsCellsOld; i++)
		countingArrayOld[i] = 0;

	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (labelArrayOld[k][x_new(i, j, Length)] > 0) {
					countingArrayOld[labelArrayOld[k][x_new(i, j, Length)]]++;
				}
			}
		}
	}
	int numberOfRegionsOld = 0;
	for (i = 0; i < numberOfRegionsCellsOld; i++) {
		if (countingArrayOld[i] > 0) {
			numberOfRegionsOld++;
		}
	}
	printf("Number of regions = %d\n", numberOfRegionsOld);

	FILE* fileOld;
	if (fopen_s(&fileOld, "output/Labels_Count_Old.txt", "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(fileOld, " 3D real image with %d Slices \n", Height);
	fprintf(fileOld, " Labelling function takes = %.5lf seconds\n", (finishOld - startOld) / CLOCKS_PER_SEC);
	fprintf(fileOld, " Number of Objects Cells = %d \n", numberOfRegionsCellsOld);
	fprintf(fileOld, " Number of Regions = %d \n", numberOfRegionsOld);
	fprintf(fileOld, "#############################################################\n");
	fprintf(fileOld, "\n\n");
	fprintf(fileOld, "Region Label        Count");
	for (i = 0; i < numberOfRegionsCellsOld; i++) {
		if (countingArrayOld[i] > 0) {
			fprintf(fileOld, "\n");
			fprintf(fileOld, "   %d                 %d", i, countingArrayOld[i]);
		}
	}
	fclose(fileOld);

	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (labelArrayOld[k][x_new(i, j, Length)] == 0) {
					labelArrayOld[k][x_new(i, j, Length)] = background;
				}
			}
		}
	}

	unsigned char segOld[] = "output/segOld.raw";
	manageFile(labelArrayOld, Length, Width, Height, segOld, operation, dType, flags);

	for (k = 0; k < Height; k++) {
		free(imageOld[k]); free(labelArrayOld[k]);
	}
	free(imageOld); free(labelArrayOld); free(countingArrayOld);

	for (k = 0; k < Height; k++) {
		free(segmented[k]);
	}
	free(segmented); free(labelArrayNew); free(imageNew); free(status); free(countingArrayNew);


	//#####################################################################################
	for (k = 0; k < Height; k++) {
		free(image[k]);
	}
	free(image);

	return EXIT_SUCCESS;
}
