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
#include "../src/generate_3d_shapes.h"
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

	dataType** imageData = (dataType**)malloc(Height * sizeof(dataType*));
	short** image = (short**)malloc(Height * sizeof(short*));
	for (k = 0; k < Height; k++) {
		imageData[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		image[k] = (short*)malloc(dim2D * sizeof(short));
	}
	if (imageData == NULL || image == NULL) return false;

	//initialization
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				imageData[k][x_new(i, j, Length)] = 0;
				image[k][x_new(i, j, Length)] = 0;
			}
		}
	}

	//Loading
	unsigned char pathPtr[] = "data/3D/patient1b.raw";
	OperationType operation = LOAD_DATA_RAW;
	LoadDataType dType = BINARY_DATA;
	Storage_Flags flags = { false, false };
	manageFile(imageData, Length, Width, Height, pathPtr, operation, dType, flags);
	load3dArrayRAW<short>(image, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/input/patient1b.raw");
	//store3dRawData<short>(image, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/loaded.raw");

	////Create artificial image
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (sqrt((256 - i) * (256 - i) + (256 - j) * (256 - j) + (254 - k) * (254 - k)) < 250) {
	//				imageData[k][x_new(i, j, Length)] = 1;
	//			}
	//		}
	//	}
	//}
	//store3dRawData<dataType>(imageData, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/bigBall.raw");

	////loop for erosion
	//char pathSave[200]; int n;
	//for (n = 1; n < 40; n++) {
	//	erosion3D(imageData, Length, Width, Height, 1, 0);
	//	sprintf_s(pathSave, "C:/Users/Konan Allaly/Documents/Tests/output/Erosion/erosion0%d.raw", n);
	//	store3dRawData<dataType>(imageData, Length, Width, Height, pathSave);
	//}
	////Loop for dilatation
	//for (n = 1; n <= 40; n++) {
	//	dilatation3D(imageData, Length, Width, Height, 1, 0);
	//	sprintf_s(pathSave, "C:/Users/Konan Allaly/Documents/Tests/output/Dilatation/dilatation0%d.raw", n);
	//	store3dRawData<dataType>(imageData, Length, Width, Height, pathSave);
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
				ImageData.imageDataPtr[k][x_new(i, j, Length)] = (dataType)image[k][x_new(i, j, Length)];
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

	/*size_t x = Length / 2, y = Width / 2;
	ImageData.imageDataPtr[0][x_new(x - 1, y - 1, Length)] = 0;
	ImageData.imageDataPtr[0][x_new(x - 1, y, Length)] = 0;
	ImageData.imageDataPtr[0][x_new(x - 1, y + 1, Length)] = 0;
	ImageData.imageDataPtr[0][x_new(x, y - 1, Length)] = 0;
	ImageData.imageDataPtr[0][x_new(x, y, Length)] = 0;
	ImageData.imageDataPtr[0][x_new(x, y + 1, Length)] = 0;
	ImageData.imageDataPtr[0][x_new(x + 1, y - 1, Length)] = 0;
	ImageData.imageDataPtr[0][x_new(x + 1, y, Length)] = 0;
	ImageData.imageDataPtr[0][x_new(x + 1, y + 1, Length)] = 0;
	store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "output/loadedImage.raw");*/

	/*dataType* arraytest = (dataType*) malloc(dim2D * sizeof(dataType));
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				arraytest[x_new(i, j, Length)] = ImageData.imageDataPtr[k][x_new(i, j, Length)];
			}
		}
	}*/

	//NoiseParameters Nparameters = { 0.05, 0, 0, minData, maxData };
	//const NoiseType Ntype = SALT_AND_PEPPER;
	////addNoiseToImage(ImageData.imageDataPtr, Length, Width, Height, Nparameters, Ntype);
	//saltAndPepper2dNoise_D(arraytest, Length, Width, 0.05, maxData);
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			ImageData.imageDataPtr[k][x_new(i, j, Length)] = arraytest[x_new(i, j, Length)];
	//		}
	//	}
	//}
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "output/noisyImage2D.raw");

	//Test of effect of rescalling
	rescaleNewRange(ImageData.imageDataPtr, Length, Width, Height, 0, 1);
	rescaleNewRange(ImageData.imageDataPtr, Length, Width, Height, minData, maxData);
	store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/rescalled.raw");

	////Compute and save histogram
	//short totalClass = (short)(maxData - minData + 1);
	//dataType* histogram = (dataType*)malloc(totalClass * sizeof(dataType));
	//for (i = 0; i < totalClass; i++) histogram[i] = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			histogram[image[k][x_new(i, j, Length)]]++;
	//		}
	//	}
	//}
	/*FILE* histo;
	if (fopen_s(&histo, "output/histogram.csv", "w") != 0) {
		printf("Enable to open");
		return false;
	}
	for (i = 0; i < totalClass; i++) {
		fprintf(histo, "%f \n", histogram[i]);
	}
	fclose(histo);*/

	//dataType numberOfCells = 0;
	//for (i = 0; i < totalClass; i++) {
	//	numberOfCells = numberOfCells + histogram[i];
	//}

	////Compute probability
	//dataType* Proba = (dataType*)malloc(totalClass * sizeof(dataType));
	//dataType sumProba = 0;
	//for (i = 0; i < totalClass; i++) Proba[i] = 0;
	//for (i = 0; i < totalClass; i++) {
	//	Proba[i] = histogram[i] / numberOfCells;
	//	sumProba = sumProba + Proba[i];
	//}
	//
	//dataType tol = 0.1, sigma = 0; short T = (short)minData;
	//dataType weightClass_b, meanClass_b, weightClass_f, meanClass_f;
	//do {
	//	weightClass_b = 0, meanClass_b = 0, weightClass_f = 0, meanClass_f = 0;
	//	//compute class weight
	//	for (i = 0; i <= T; i++) {
	//		weightClass_b = weightClass_b + Proba[i];
	//	}
	//	for (i = T + 1; i < totalClass; i++) {
	//		weightClass_f = weightClass_f + Proba[i];
	//	}
	//	//compute class mean
	//	for (i = 0; i <= T; i++) {
	//		meanClass_b = meanClass_b + i * Proba[i];
	//	}
	//	meanClass_b = meanClass_b / weightClass_b;
	//	for (i = T + 1; i < totalClass; i++) {
	//		meanClass_f = meanClass_f + i * Proba[i];
	//	}
	//	////compute class variance
	//	//for (i = 0; i <= T; i++) {
	//	//	varClass_b = Class_b + i * Proba[i];
	//	//}
	//	meanClass_f = meanClass_f / weightClass_b;
	//	//compute inter-class variance
	//	sigma = weightClass_b * meanClass_b + weightClass_f * meanClass_f;

	//	T = T + 1;
	//} while (T < (short)maxData && abs(sigma - T) > tol);	

	//rescaleNewRange(ImageData.imageDataPtr, Length, Width, Height, 0, 1);
	//NoiseParameters Nparameters = { 0.05, 0.1, 10, 0, 1 };
	//const NoiseType Ntype = SALT_AND_PEPPER;
	//addNoiseToImage(ImageData.imageDataPtr, Length, Width, Height, Nparameters, Ntype);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "output/noisyImage.raw");

	////Filtering by geodesic mean curvature filter
	//const size_t maxNumberOfSolverIteration = 1000;
	//dataType  coef = 0.01, eps2 = 1e-4;
	//size_t numberOfTimeStep = 1;
	//Filter_Parameters GMC_filterParameters;
	//const FilterMethod methodFiltering = GEODESIC_MEAN_CURVATURE_FILTER;
	//dataType timeStepSize = 0.002, h = 1, sigma = 0.01, K = 0.018, omega_c = 1.5, tolerance = 10;
	//size_t p = 1, timeStepsNum = 1, maxNumberOftimeSteps = 100;
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
	
	//Thresholding
	//thresholding3dFunctionN(ImageData.imageDataPtr, Length, Width, Height, thresmin, thresmax, minData, maxData);
	//thresholding3dFunctionN(ImageData.imageDataPtr, Length, Width, Height, 1279, 1379, minData, maxData);
	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/thres.raw");

	////Loop for Erosion
	//char pathSaveErosion[200]; int n;
	//for (n = 1; n < 11; n++) {
	//	erosion3D(ImageData.imageDataPtr, Length, Width, Height, maxData, minData);
	//	//sprintf_s(pathSaveErosion, "C:/Users/Konan Allaly/Documents/Tests/output/F_erosion0%d.raw", n);
	//	//store3dRawData<dataType>(ImageData.imageDataPtr, Length, Width, Height, pathSaveErosion);
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
		free(image[k]); free(imageData[k]);
		free(ImageData.imageDataPtr[k]);
		//free(labelArray[k]); free(status[k]);
	}
	free(imageData); free(image);
	free(ImageData.imageDataPtr);
	//free(labelArray); free(status); free(countingArray);
	//free(cenTroid);
	//free(P); free(histogram);
	//free(arraytest);

	return EXIT_SUCCESS;
}
