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

	///////////////////////
	//size_t kMin = 332, iMin = 126, jMin = 180, kn, in, jn;
	//const size_t heightNew = 170, lengthNew = 180, widthNew = 180, dim2dNew = lengthNew * widthNew;

	//dataType** croppedImage = (dataType**)malloc(sizeof(dataType*) * heightNew);
	//dataType** maskThreshold = (dataType**)malloc(sizeof(dataType*) * heightNew);
	//dataType** distanceMap = (dataType**)malloc(sizeof(dataType*) * heightNew);
	//dataType** initialSegment = (dataType**)malloc(sizeof(dataType*) * heightNew);
	//for (k = 0; k < heightNew; k++) {
	//	croppedImage[k] = (dataType*)malloc(sizeof(dataType) * dim2dNew);
	//	distanceMap[k] = (dataType*)malloc(sizeof(dataType) * dim2dNew);
	//	maskThreshold[k] = (dataType*)malloc(sizeof(dataType) * dim2dNew);
	//	initialSegment[k] = (dataType*)malloc(sizeof(dataType) * dim2dNew);
	//}
	//if (croppedImage == NULL || distanceMap == NULL || maskThreshold == NULL || initialSegment == NULL) {
	//	return false;
	//}
	//const dataType initialSegmentationValue = 1.0;

	//for (k = 0, kn = kMin; k < heightNew; k++, kn++) {
	//	for (i = 0, in = iMin; i < lengthNew; i++, in++) {
	//		for (j = 0, jn = jMin; j < widthNew; j++, jn++) {
	//			x = x_new(i, j, lengthNew);
	//			croppedImage[k][x] = resampledImage[kn][x_new(in, jn, Length)];
	//			maskThreshold[k][x] = croppedImage[k][x];
	//			initialSegment[k][x] = initialSegmentationValue;
	//		}
	//	}
	//}
	///////////////////////

	
	size_t length_new = 180, width_new = 180, height_new = 170;
	dataType** croppedImage = new dataType * [height_new];
	dataType** imageThresh = new dataType * [height_new];
	dataType** distanceMap = new dataType * [height_new];
	for (k = 0; k < height_new; k++) {
		croppedImage[k] = new dataType[length_new * width_new];
		imageThresh[k] = new dataType[length_new * width_new];
		distanceMap[k] = new dataType[length_new * width_new];
	}
	if (croppedImage == NULL || imageThresh == NULL || distanceMap == NULL) {
		return false;
	}

	//cropping
	size_t i_ext, j_ext, k_ext;
	size_t kMin = 332, iMin = 126, jMin = 180;
	for (k = 0, k_ext = kMin; k < height_new; k++, k_ext++) {
		for (i = 0, i_ext = iMin; i < length_new; i++, i_ext++) {
			for (j = 0, j_ext = jMin; j < width_new; j++, j_ext++) {
				xd = x_new(i, j, length_new);
				croppedImage[k][xd] = imageData[k_ext][x_new(i_ext, j_ext, Length)];
				imageThresh[k][xd] = imageData[k_ext][x_new(i_ext, j_ext, Length)];
				distanceMap[k][xd] = 0.0;
			}
		}
	}

	//load3dArrayRAW<dataType>(imageData, Length, Width, Height, inputImagePath.c_str());
	std::string loaded3D = outputPath + "croppedImage.raw";
	store3dRawData<dataType>(croppedImage, length_new, width_new, height_new, loaded3D.c_str());

	//Freeing image pointer 
	for (k = 0; k < Height; k++) {
		delete[] image[k];
		delete[] imageData[k];
	}
	delete[] image;
	delete[] imageData;

	size_t nb_centers = 1;
	Point3D* centerSeg = new Point3D[nb_centers];

	//-------------Filtering-------------------------

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
	//std::string filteredImagePath = outputPath + "filtered1bMC.raw";
	//store3dRawData<dataType>(imageData, Length, Width, Height, filteredImagePath.c_str());

	//--------------- Generate artificial Image -----

	//// Artificial Sphere
	//size_t x = Length / 2, y = Width / 2, z = Height / 2;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (sqrt((x - i)*(x - i) + (y - j)*(y - j) + (z - k)*(z - k)) <= 50) {
	//				imageData[k][x_new(i, j, Length)] = 1;
	//			}
	//			else {
	//				imageData[k][x_new(i, j, Length)] = 0;
	//			}
	//		}
	//	}
	//}

	//Point3D c = { x,y,z };
	//dataType radius = 20;
	//dataType small_radius = 6;
	//size_t pitch_ball = 30;
	//dataType fill_value = 1.0;
	//std::string artificialImage = outputPath + "ballsWithHoles.raw";

	//generateSphere(imageData, c, Length, Width, Height, radius, fill_value, artificialImage.c_str());
	//generateSphereWithSixHoles(imageData, c, Length, Width, Height, radius, small_radius, fill_value, artificialImage.c_str());
	//ballsOnHelix(imageData, pitch_ball, Length, Width, Height, radius, small_radius, fill_value);
	//store3dRawData<dataType>(imageData, Length, Width, Height, artificialImage.c_str());

	//--------------- Distance Map

	Distance_Map_Params distParameters = {0.4, 1.0, 0.0, 100000, 0.001};
	DistanceMapMethod method = FAST_SWEEP;
	thresholding3dFunctionN(imageThresh, length_new, width_new, height_new, thresmin, thresmax, 0.0, 1.0);
	computeDistanceMap(distanceMap, imageThresh, length_new, width_new, height_new, distParameters, method);

	loaded3D = outputPath + "distanceMap.raw";
	store3dRawData<dataType>(distanceMap, length_new, width_new, height_new, loaded3D.c_str());

	dataType distance_max = 0.0;

	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			xd = x_new(i, j, Length);
	//			if (distanceMap[k][xd] > distance_max) {
	//				distance_max = distanceMap[k][xd];
	//				centerSeg->x = i;
	//				centerSeg->y = j;
	//				centerSeg->z = k;
	//			}
	//		}
	//	}
	//}

	//--------------- Test segmentation -------------

	//dataType** initialSegment = new dataType * [Height];
	//for (k = 0; k < Height; k++) {
	//	initialSegment[k] = new dataType[Width * Length];
	//}

	////If we want to start with the segmentatation function originally implemented in the library
	//generateInitialSegmentationFunctionForMultipleCentres(initialSegment, Length, Width, Height, centerSeg, 0.5, 30, nb_centers);

	//std::string segmFolderPath = outputPath + "segmentation/";
	//store3dRawData<dataType>(initialSegment, Length, Width, Height, (segmFolderPath + std::string("_seg_func_000.raw")).c_str());

	//Image_Data segment; segment.height = Height; segment.length = Length; segment.width = Width; segment.imageDataPtr = imageData;
	//rescaleNewRange(segment.imageDataPtr, Length, Width, Height, 0.0, 1.0);
	Segmentation_Parameters segmentParameters; segmentParameters.coef = 10000; segmentParameters.eps2 = 1e-6; segmentParameters.gauss_seidelTolerance = 1e-3;
	segmentParameters.h = 1.0; segmentParameters.maxNoGSIteration = 100; segmentParameters.maxNoOfTimeSteps = 100; segmentParameters.mod = 1;
	segmentParameters.numberOfTimeStep = 100; segmentParameters.omega_c = 1.5; segmentParameters.segTolerance = 1e-8; segmentParameters.tau = 1.0;

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
