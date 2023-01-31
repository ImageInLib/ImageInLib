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
#include "../src/segmentation3d_gsubsurf.h"
#include "../src/image_difference.h"
#include "../src/noise_generator.h"
#include "../src/image_norm.h"
#include "../src/imageInterpolation.h"

#define originalMean 71.1245
#define offSet 1024
#define standarDeviation 22.001
#define thresmin 995
#define thresmax 1213
#define minimalSize 2000


int main() {

	size_t i, j, k, x;

	//image Dimensions
	const size_t Width = 512;
	const size_t Length = 512;
	const size_t Height = 406;
	const size_t dim2D = Width * Length;
	/*The Width and Lenght are the same for all the patients
	Height : patient1b = patient3 = 508, patient2 = 406, patient6 = 880, patient7 = 607*/
	
	//Preparation of image data pointers
	dataType** imageData = (dataType**)malloc(Height * sizeof(dataType*));
	short** image = (short**)malloc(Height * sizeof(short*));
	for (k = 0; k < Height; k++) {
		imageData[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		image[k] = (short*)malloc(dim2D * sizeof(short));
	}

	if (imageData == NULL || image == NULL) 
		return false;

	std::string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/";
	std::string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";
	std::string inputImagePath = inputPath +  "patient2.raw"; // "Slices304.raw";

	if (load3dArrayRAW<short>(image, Length, Width, Height, inputImagePath.c_str()) == false)
	{
		printf("inputImagePath does not exist\n");
	}

	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				x = x_new(i, j, Length);
				imageData[k][x] = (dataType)image[k][x];
			}
		}
	}

	const size_t HeightNew = 866;
	dataType ** resampledImage = (dataType**)malloc(HeightNew * sizeof(dataType*));
	dataType** distanceMap = (dataType**)malloc(HeightNew * sizeof(dataType*));
	dataType** maskThresh = (dataType**)malloc(HeightNew * sizeof(dataType*));
	for (k = 0; k < HeightNew; k++) {
		resampledImage[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		distanceMap[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		maskThresh[k] = (dataType*)malloc(dim2D * sizeof(dataType));
	}
	if (resampledImage == NULL || distanceMap == NULL || maskThresh == NULL)
		return false;

	dataType k_spacingOld = 2.5, k_spacingNew = 1.171875;
	linear2dInterpolation(imageData, resampledImage, Length, Width, Height, k_spacingOld, k_spacingNew);

	for (k = 0; k < HeightNew; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				x = x_new(i, j, Length);
				maskThresh[k][x] = resampledImage[k][x];
				distanceMap[k][x] = 0;
			}
		}
	}

	//std::string loadedImagePath = outputPath + "resampled.raw";
	//store3dRawData<dataType>(resampledImage, Length, Width, HeightNew, loadedImagePath.c_str());

	thresholding3dFunctionN(maskThresh, Length, Width, HeightNew, thresmin, thresmax, 0, 1);
	//std::string thresholdedImagePath = outputPath + "thresholded.raw";
	//store3dRawData<dataType>(maskThresh, Length, Width, Height, thresholdedImagePath.c_str());

	//erosion3dHeighteenNeigbours(maskThresh, Length, Width, HeightNew, 1, 0);
	//erosion3dHeighteenNeigbours(maskThresh, Length, Width, Height, 1, 0);
	//std::string erodedImagePath = outputPath + "erodedV2.raw";
	//store3dRawData<dataType>(imageData, Length, Width, Height, erodedImagePath.c_str());

	for (k = 0; k < HeightNew; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				x = x_new(i, j, Length);
				if (maskThresh[k][x] == 0) {
					resampledImage[k][x] = 0;
				}
			}
		}
	}

	fastSweepingFunction_3D(distanceMap, maskThresh, Length, Width, HeightNew, 1, 100000000, 0);
	//std::string distanceImagePath = outputPath + "distanceMap.raw";
	//store3dRawData<dataType>(distanceMap, Length, Width, HeightNew, distanceImagePath.c_str());

	//finding of a point with the highest distance
	dataType distanceMax = -1;
	int i_max = 0, j_max = 0, k_max = 0;

	for (k = 0; k < HeightNew; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				x = x_new(i, j, Length);
				if (distanceMap[k][x] >= distanceMax) {
					distanceMax = distanceMap[k][x];
					i_max = (int)i; j_max = (int)j; k_max = (int)k;
				}
			}
		}
	}

	printf("Maximal distance : %f \n", distanceMax);
	printf("Coordinates of the highest distance x = %d, y = %d and z = %d \n", i_max, j_max, k_max);

	//Region growing
	dataType** segmented = (dataType**)malloc(HeightNew * sizeof(dataType*));
	bool** status = (bool**)malloc(HeightNew * sizeof(bool*));
	for (k = 0; k < HeightNew; k++) {
		segmented[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		status[k] = (bool*)malloc(dim2D * sizeof(bool));
	}
	if (segmented == NULL || status == NULL)
		return false;

	for (k = 0; k < HeightNew; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				x = x_new(i, j, Length);
				segmented[k][x] = 0;
				status[k][x] = false;
			}
		}
	}

	//printf("seed voxel intensity : %f \n", imageData[k_max][x_new(i_max, j_max, Length)]);
	//printf("top voxel intensity : %f \n", imageData[k_max - 1][x_new(i_max, j_max, Length)]);
	//printf("bottom voxel intensity : %f \n", imageData[k_max + 1][x_new(i_max, j_max, Length)]);
	//printf("... voxel intensity : %f \n", imageData[k_max][x_new(i_max - 1, j_max, Length)]);
	//printf("... voxel intensity : %f \n", imageData[k_max][x_new(i_max + 1, j_max, Length)]);
	//printf("... voxel intensity : %f \n", imageData[k_max][x_new(i_max, j_max - 1, Length)]);
	//printf("... voxel intensity : %f \n", imageData[k_max][x_new(i_max, j_max + 1, Length)]);

	Point3D* seed = (Point3D*)malloc(sizeof(Point3D));
	seed->x = i_max; seed->y = j_max; seed->z = k_max;
	regionGrowing(resampledImage, segmented, status, Length, Width, HeightNew, thresmin, thresmax, seed);

	//std::string segmentedImagePath = outputPath + "segmentedImageV3.raw";
	//store3dRawData<dataType>(segmented, Length, Width, HeightNew, segmentedImagePath.c_str());

	//------------------------------------------------------------------------------------------------------

	//Manual croppping

	//kMin depends on the interpolated image
	size_t kMin = 332, iMin = 126, jMin = 180, kn, in, jn;
	const size_t heightNew = 170, lengthNew = 180, widthNew = 180, dim2dNew = lengthNew * widthNew;

	dataType** croppedImage = (dataType**)malloc(sizeof(dataType*) * heightNew);
	dataType** maskThresholdCropped = (dataType**)malloc(sizeof(dataType*) * heightNew);
	dataType** distanceMapCropped = (dataType**)malloc(sizeof(dataType*) * heightNew);
	dataType** initialSegment = (dataType**)malloc(sizeof(dataType*) * heightNew);
	for (k = 0; k < heightNew; k++) {
		croppedImage[k] = (dataType*)malloc(sizeof(dataType) * dim2dNew);
		distanceMapCropped[k] = (dataType*)malloc(sizeof(dataType) * dim2dNew);
		maskThresholdCropped[k] = (dataType*)malloc(sizeof(dataType) * dim2dNew);
		initialSegment[k] = (dataType*)malloc(sizeof(dataType) * dim2dNew);
	}
	if (croppedImage == NULL || distanceMapCropped == NULL || maskThresholdCropped == NULL || initialSegment == NULL) {
		return false;
	}

	const dataType initialSegmentationValue = 1.0;

	for (k = 0, kn = kMin; k < heightNew; k++, kn++) {
		for (i = 0, in = iMin; i < lengthNew; i++, in++) {
			for (j = 0, jn = jMin; j < widthNew; j++, jn++) {
				x = x_new(i, j, lengthNew);
				croppedImage[k][x] = resampledImage[kn][x_new(in, jn, Length)];
				maskThresholdCropped[k][x] = croppedImage[k][x];
				initialSegment[k][x] = initialSegmentationValue;
			}
		}
	}

	//std::string croppedVolumePath = outputPath + "cropped.raw";
	//store3dRawData<dataType>(croppedImage, lengthNew, widthNew, heightNew, croppedVolumePath.c_str());

	//------------------------------------------------------------------------------------------------------------

	//Fast sweeping to find the point with the higest distance
	thresholding3dFunctionN(maskThresholdCropped, lengthNew, widthNew, heightNew, thresmin, thresmax, 0, 1);
	fastSweepingFunction_3D(distanceMapCropped, maskThresholdCropped, lengthNew, widthNew, heightNew, 1, 100000000, 0);

	//finding of the point with the highest distance
	distanceMax = -1;

	for (k = 0; k < heightNew; k++) {
		for (i = 0; i < lengthNew; i++) {
			for (j = 0; j < widthNew; j++) {
				if (distanceMapCropped[k][x_new(i, j, lengthNew)] >= distanceMax) {
					distanceMax = distanceMapCropped[k][x_new(i, j, lengthNew)];
					i_max = (int)i; j_max = (int)j; k_max = (int)k;
				}
			}
		}
	}
	
	printf("Maximal distance for the cropped volume : %f \n", distanceMax);
	printf("Coordinates of the highest distance (Cropped Volume) : x = %d, y = %d and z = %d \n", i_max, j_max, k_max);

	//------------------------------------------------------------------------------------------------------------
	 
	//Segmentation parameters for real image
	size_t numb_centers = 1; Point3D * centerSeg = (Point3D*)malloc(sizeof(Point3D) * numb_centers);
	centerSeg->x = i_max; centerSeg->y = j_max; centerSeg->z = k_max; //---> used for one center
	//used for multiple centers
	//centerSeg[0].x = i_max; centerSeg[0].y = j_max; centerSeg[0].z = k_max;
	//centerSeg[1].x = i_max + 50; centerSeg[1].y = j_max - 70; centerSeg[1].z = k_max;

	//If we want to start with the segmentatation function originally implemented in the library
	generateInitialSegmentationFunctionForMultipleCentres(initialSegment, lengthNew, widthNew, heightNew, centerSeg, 0.5, 15, numb_centers);

	std::string segmFolderPath = outputPath + "segmentation/";
	store3dRawData<dataType>(initialSegment, lengthNew, widthNew, heightNew, (segmFolderPath + std::string("_seg_func_000.raw")).c_str());

	Image_Data segment; segment.height = heightNew; segment.length = lengthNew; segment.width = widthNew; segment.imageDataPtr = croppedImage;
	rescaleNewRange(segment.imageDataPtr, lengthNew, widthNew, heightNew, 0, 1);
	Segmentation_Parameters segmentParameters; segmentParameters.coef = 10000; segmentParameters.eps2 = 1e-6; segmentParameters.gauss_seidelTolerance = 1e-3;
	segmentParameters.h = k_spacingNew; segmentParameters.maxNoGSIteration = 100; segmentParameters.maxNoOfTimeSteps = 300; segmentParameters.mod = 1;
	segmentParameters.numberOfTimeStep = 300; segmentParameters.omega_c = 1.5; segmentParameters.segTolerance = 1e-4; segmentParameters.tau = 4;
	Filter_Parameters filterParameters; filterParameters.coef = 1e-6; filterParameters.edge_detector_coefficient = 100; filterParameters.eps2 = 1e-6;
	filterParameters.h = k_spacingNew; filterParameters.maxNumberOfSolverIteration = 100; filterParameters.omega_c = 1.5; filterParameters.p = 1;
	filterParameters.sigma = 1e-3; filterParameters.timeStepSize = 1.2; filterParameters.timeStepsNum = 1; filterParameters.tolerance = 1e-3;

	unsigned char outputPathPtr[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/";
	subsurfSegmentation(segment, initialSegment, segmentParameters, filterParameters, centerSeg, numb_centers, outputPathPtr);
	//generalizedSubsurfSegmentation(segment, initialSegment, segmentParameters, filterParameters, centerSeg, numb_centers, outputPathPtr, -1.0, 1.0);

	//------------------------------------------------------------------------------------------------

	//free memory
	for (k = 0; k < Height; k++) {
		free(imageData[k]); free(image[k]);
	}
	free(imageData); free(image);

	for (k = 0; k < HeightNew; k++) {
		free(resampledImage[k]);
		free(distanceMap[k]); free(maskThresh[k]);
		free(segmented[k]); free(status[k]);
	}
	free(resampledImage);
	free(distanceMap); free(maskThresh);
	free(segmented); free(status);

	free(seed);

	//for (k = 0; k < zDim; k++) {
	//	free(resampledImage[k]); free(resampledLiver[k]);
	//}
	//free(resampledImage); free(resampledLiver);

	for (k = 0; k < heightNew; k++) {
		free(croppedImage[k]); free(distanceMapCropped[k]); free(maskThresholdCropped[k]); free(initialSegment[k]);
	}
	free(croppedImage); free(distanceMapCropped); free(maskThresholdCropped); free(initialSegment);
	 
	free(centerSeg);

	return EXIT_SUCCESS;
}
