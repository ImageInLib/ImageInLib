#pragma warning(disable : 6011)
#pragma warning(disable : 4244)
#pragma warning(disable : 4477)
#pragma warning(disable : 4267)
#pragma warning(disable : 4700)

#include <iostream>
#include<vector>
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
#include "distanceForPathFinding.h"
#include "distanceMaps.h"
#include "../src/trajectories.h"

#define thresmin 995
#define thresmax 1213


int main() {

	size_t i, j, k, xd, m, n;

	//image Dimensions
	const size_t Width = 512;
	const size_t Length = 512;
	const size_t Height = 508; /*607*/ /*508*/
	const size_t dim2D = Width * Length;

	//-------------Real 3D image -------------------------
	
	//Preparation of image data pointers
	dataType** imageData = new dataType * [Height];
	short** image = new short * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType [dim2D];
		image[k] = new short [dim2D];
	}
	if (imageData == NULL || image == NULL) {
		return false;
	}
	std::string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/";
	std::string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";

	//std::string inputImagePath = inputPath + "patient2.raw";
	////std::string inputImagePath = inputPath + "patient2_filtered.raw";
	////std::string inputImagePath = inputPath + "filteredK100.raw";
	////std::string inputImagePath = inputPath + "filteredK1.raw"; filteredHeatEQ
	//std::string inputImagePath = inputPath + "filteredHeatEQ.raw";
	std::string inputImagePath = inputPath + "patient1b.raw";
	//std::string inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/output/filteredP1bMC.raw";
	//std::string inputImagePath = inputPath + "patient7.raw";
	//std::string inputImagePath = inputPath + "Slices304.raw";
	//std::string inputImagePath = inputPath + "filteredMC.raw";

	if (load3dArrayRAW<short>(image, Length, Width, Height, inputImagePath.c_str()) == false)
	{
		printf("inputImagePath does not exist\n");
	}
	
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				xd = x_new(i, j, Length);
				imageData[k][xd] = (dataType)image[k][xd];
			}
		}
	}

	//std::string loadedImagePath = outputPath + "loaded.raw";
	//store3dRawData<dataType>(imageData, Length, Width, Height, loadedImagePath.c_str());

	const size_t HeightNew = 866;
	dataType ** resampledImage = (dataType**)malloc(HeightNew * sizeof(dataType*));

	for (k = 0; k < HeightNew; k++) {
		resampledImage[k] = (dataType*)malloc(dim2D * sizeof(dataType));
	}
	if (resampledImage == NULL)
		return false;

	dataType k_spacingOld = 2.5, k_spacingNew = 1.171875;
	//linear2dInterpolation(imageData, resampledImage, Length, Width, Height, k_spacingOld, k_spacingNew);
	//std::string interpolatedImagePath = outputPath + "interpolated.raw";
	//store3dRawData<dataType>(imageData, Length, Width, Height, interpolatedImagePath.c_str());

	//------------- Filtering -------------------------------------------
	 
	//rescaleNewRange(imageData, Length, Width, Height, 0, 1);
	//FilterMethod method = MEAN_CURVATURE_FILTER; //NONLINEAR_HEATEQUATION_IMPLICIT; // MEAN_CURVATURE_FILTER; // LINEAR_HEATEQUATION_IMPLICIT; // 
	//Filter_Parameters filterParm; //filterParm = { 1.2, 1.0, 0.1, 1000, 1.5, 0.0004, 0.0001, 0.01, 1, 1, 1000 };
	//filterParm.timeStepSize = 1.2; filterParm.h = 1.0; filterParm.sigma = 0.1; filterParm.edge_detector_coefficient = 1;
	//filterParm.omega_c = 1.5; filterParm.tolerance = 0.0004; filterParm.eps2 = 0.0001; filterParm.coef = 0.01;
	//filterParm.timeStepsNum = 1; filterParm.p = 1; filterParm.maxNumberOfSolverIteration = 1000;
	//Image_Data toBeFiltered;
	//toBeFiltered.imageDataPtr = imageData; toBeFiltered.height = Height; toBeFiltered.length = Length; toBeFiltered.width = Width;
	//filterImage(toBeFiltered, filterParm, method);
	//rescaleNewRange(imageData, Length, Width, Height, 0, 4000);
	//thresholding3dFunctionN(imageData, Length, Width, Height, thresmin, thresmax, 0, 1);
	//erosion3dHeighteenNeigbours(imageData, Length,Width, Height, 1, 0);
	//erosion3dHeighteenNeigbours(imageData, Length, Width, Height, 1, 0);
	//erosion3dHeighteenNeigbours(imageData, Length, Width, Height, 1, 0);
	//std::string filteredImagePath = outputPath + "filteredThresshErod.raw";
	//store3dRawData<dataType>(imageData, Length, Width, Height, filteredImagePath.c_str());
	 
	//-----------------------------------------------------------------------------------------------------
	 
	////Region growing
	//dataType** segmented = (dataType**)malloc(HeightNew * sizeof(dataType*));
	//bool** status = (bool**)malloc(HeightNew * sizeof(bool*));
	//for (k = 0; k < HeightNew; k++) {
	//	segmented[k] = (dataType*)malloc(dim2D * sizeof(dataType));
	//	status[k] = (bool*)malloc(dim2D * sizeof(bool));
	//}
	//if (segmented == NULL || status == NULL)
	//	return false;

	//for (k = 0; k < HeightNew; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			x = x_new(i, j, Length);
	//			segmented[k][x] = 0;
	//			status[k][x] = false;
	//		}
	//	}
	//}

	//printf("seed voxel intensity : %f \n", imageData[k_max][x_new(i_max, j_max, Length)]);
	//printf("top voxel intensity : %f \n", imageData[k_max - 1][x_new(i_max, j_max, Length)]);
	//printf("bottom voxel intensity : %f \n", imageData[k_max + 1][x_new(i_max, j_max, Length)]);
	//printf("... voxel intensity : %f \n", imageData[k_max][x_new(i_max - 1, j_max, Length)]);
	//printf("... voxel intensity : %f \n", imageData[k_max][x_new(i_max + 1, j_max, Length)]);
	//printf("... voxel intensity : %f \n", imageData[k_max][x_new(i_max, j_max - 1, Length)]);
	//printf("... voxel intensity : %f \n", imageData[k_max][x_new(i_max, j_max + 1, Length)]);

	//Point3D* seed = (Point3D*)malloc(sizeof(Point3D));
	//seed->x = i_max; seed->y = j_max; seed->z = k_max;
	//regionGrowing(resampledImage, segmented, status, Length, Width, HeightNew, thresmin, thresmax, seed);

	//std::string segmentedImagePath = outputPath + "segmentedThreeTimesErosion.raw";
	//store3dRawData<dataType>(segmented, Length, Width, HeightNew, segmentedImagePath.c_str());

	//------------------------------------------------------------------------------------------------------

	//Manual croppping

	//kMin depends on the interpolated image
	size_t kMin = 332, iMin = 126, jMin = 180, kn, in, jn;
	const size_t heightNew = 170, lengthNew = 180, widthNew = 180, dim2dNew = lengthNew * widthNew;

	dataType** croppedImage = (dataType**)malloc(sizeof(dataType*) * heightNew);
	dataType** maskThreshold = (dataType**)malloc(sizeof(dataType*) * heightNew);
	dataType** distanceMap = (dataType**)malloc(sizeof(dataType*) * heightNew);
	dataType** initialSegment = (dataType**)malloc(sizeof(dataType*) * heightNew);
	for (k = 0; k < heightNew; k++) {
		croppedImage[k] = (dataType*)malloc(sizeof(dataType) * dim2dNew);
		distanceMap[k] = (dataType*)malloc(sizeof(dataType) * dim2dNew);
		maskThreshold[k] = (dataType*)malloc(sizeof(dataType) * dim2dNew);
		initialSegment[k] = (dataType*)malloc(sizeof(dataType) * dim2dNew);
	}
	if (croppedImage == NULL || distanceMap == NULL || maskThreshold == NULL || initialSegment == NULL) {
		return false;
	}

	const dataType initialSegmentationValue = 1.0;

	for (k = 0, kn = kMin; k < heightNew; k++, kn++) {
		for (i = 0, in = iMin; i < lengthNew; i++, in++) {
			for (j = 0, jn = jMin; j < widthNew; j++, jn++) {
				x = x_new(i, j, lengthNew);
				croppedImage[k][x] = resampledImage[kn][x_new(in, jn, Length)];
				maskThreshold[k][x] = croppedImage[k][x];
				initialSegment[k][x] = initialSegmentationValue;
			}
		}
	}

	//std::string croppedVolumePath = outputPath + "croppedImage.raw";
	//store3dRawData<dataType>(croppedImage, lengthNew, widthNew, heightNew, croppedVolumePath.c_str());

	Distance_Map_Params distParameters = {0.4, 1.0, 0.0, 100000, 0.001};
	DistanceMapMethod method = FAST_SWEEP;
	thresholding3dFunctionN(imageThresh, length_new, width_new, height_new, thresmin, thresmax, 0.0, 1.0);
	computeDistanceMap(distanceMap, imageThresh, length_new, width_new, height_new, distParameters, method);

	//Fast sweeping to find the point with the higest distance
	//thresholding3dFunctionN(maskThreshold, lengthNew, widthNew, heightNew, thresmin, thresmax, 0, 1);
	//fastSweepingFunction_3D(distanceMap, maskThreshold, lengthNew, widthNew, heightNew, 1, 100000000, 0);

	//finding of the point with the highest distance
	dataType distanceMax = -1;
	int i_max, j_max, k_max;

	for (k = 0; k < heightNew; k++) {
		for (i = 0; i < lengthNew; i++) {
			for (j = 0; j < widthNew; j++) {
				if (distanceMap[k][x_new(i, j, lengthNew)] >= distanceMax) {
					distanceMax = distanceMap[k][x_new(i, j, lengthNew)];
					i_max = (int)i; j_max = (int)j; k_max = (int)k;
				}
			}
		}
	}
	
	//printf("Maximal distance for the cropped volume : %f \n", distanceMax);
	//printf("Coordinates of the highest distance (Cropped Volume) : x = %d, y = %d and z = %d \n", i_max, j_max, k_max);

	//------------------------------------------------------------------------------------------------------------
	 
	//Segmentation parameters for real image
	size_t numb_centers = 1; Point3D* centerSeg = (Point3D*)malloc(sizeof(Point3D) * numb_centers);
	centerSeg->x = i_max; centerSeg->y = j_max; centerSeg->z = k_max; //---> used for one center
	//used for multiple centers
	//centerSeg[0].x = i_max; centerSeg[0].y = j_max; centerSeg[0].z = k_max;
	//centerSeg[1].x = i_max + 50; centerSeg[1].y = j_max - 70; centerSeg[1].z = k_max;

	////If we want to start with the segmentatation function originally implemented in the library
	//generateInitialSegmentationFunctionForMultipleCentres(initialSegment, Length, Width, Height, centerSeg, 0.5, 30, nb_centers);

	std::string segmFolderPath = outputPath + "segmentation/";
	//store3dRawData<dataType>(initialSegment, lengthNew, widthNew, heightNew, (segmFolderPath + std::string("_seg_func_000.raw")).c_str());

	Image_Data segment; segment.height = heightNew; segment.length = lengthNew; segment.width = widthNew; segment.imageDataPtr = croppedImage;
	rescaleNewRange(segment.imageDataPtr, lengthNew, widthNew, heightNew, 0, 1);
	Segmentation_Parameters segmentParameters; segmentParameters.coef = 10000; segmentParameters.eps2 = 1e-6; segmentParameters.gauss_seidelTolerance = 1e-3;
	segmentParameters.h = 1.0; segmentParameters.maxNoGSIteration = 100; segmentParameters.maxNoOfTimeSteps = 300; segmentParameters.mod = 1;
	segmentParameters.numberOfTimeStep = 300; segmentParameters.omega_c = 1.5; segmentParameters.segTolerance = 1e-4; segmentParameters.tau = 8;
	
	Filter_Parameters filterParameters; filterParameters.coef = 1e-2; filterParameters.edge_detector_coefficient = 1; filterParameters.eps2 = 1e-4;
	filterParameters.h = 1.0; filterParameters.maxNumberOfSolverIteration = 100; filterParameters.omega_c = 1.5; filterParameters.p = 1;
	filterParameters.sigma = 1e-3; filterParameters.timeStepSize = 1.2; filterParameters.timeStepsNum = 1; filterParameters.tolerance = 4*1e-4;

	Filter_Parameters filterParameters; filterParameters.coef = 1e-2; filterParameters.edge_detector_coefficient = 1; filterParameters.eps2 = 1e-6;
	filterParameters.h = 1.0; filterParameters.maxNumberOfSolverIteration = 100; filterParameters.omega_c = 1.5; filterParameters.p = 1;
	filterParameters.sigma = 1e-3; filterParameters.timeStepSize = 1.0; filterParameters.timeStepsNum = 1; filterParameters.tolerance = 4 * 1e-4;

	unsigned char outputPathPtr[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/";
	//subsurfSegmentation(segment, initialSegment, segmentParameters, filterParameters, centerSeg, nb_centers, outputPathPtr);
	//generalizedSubsurfSegmentation(segment, initialSegment, segmentParameters, filterParameters, centerSeg, nb_centers, outputPathPtr, 1.0, 1.0);

	//delete[] centerSeg;
	for (k = 0; k < height_new; k++) {
		//delete[] imageData[k];
		delete[] imageThresh[k];
		delete[] distanceMap[k];
		//delete[] initialSegment[k];
		delete[] croppedImage[k];
	}
	//delete[] imageData;
	delete[] imageThresh;
	delete[] distanceMap;
	delete[] croppedImage;
	//delete[] initialSegment;

	return EXIT_SUCCESS;
}
