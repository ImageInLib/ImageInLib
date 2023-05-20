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
#include "../src/segmentation3D_atlas.h"
#include "segmentation.h"

#define thresmin 995
#define thresmax 1213


int main() {

	size_t i, j, k, xd, m, n, x;

	//image Dimensions
	const size_t Width = 128; //256; //512; //415; //512;
	const size_t Length = 128; //256; //512; //279; //512;
	const size_t Height = 217; //433; //406; /*607*/ /*508*/
	const size_t dim2D = Width * Length;

	//-------------Real 3D image -------------------------
	
	//Preparation of image data pointers
	dataType** imageData = new dataType * [Height];
	//dataType** distanceMap = new dataType * [Height];
	short** image = new short * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType [dim2D];
		//distanceMap[k] = new dataType[dim2D];
		image[k] = new short [dim2D];
	}
	if (imageData == NULL || image == NULL) {
		return false;
	}

	std::string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/";
	std::string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";

	//std::string inputImagePath = inputPath + "filteredHeatEQ.raw";
	//std::string inputImagePath = inputPath + "patient2_downsampled.raw";
	std::string inputImagePath = inputPath + "patient2_down_down.raw";

	//Operation operation = LOAD_DATA;
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, inputImagePath.c_str(), operation, false);

	//std::string loadedImagePath = outputPath + "loadedFloat.raw";
	//operation = STORE_DATA;
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, loadedImagePath.c_str(), operation, false);

	//for (k = 0; k < Height; k++) {
	//	delete[] image[k];
	//	delete[] imageData[k];
	//}
	//delete[] image;

	////std::string inputImagePath = inputPath + "patient2_filtered.raw";
	////std::string inputImagePath = inputPath + "filteredK100.raw";
	////std::string inputImagePath = inputPath + "filteredK1.raw"; filteredHeatEQ
	//std::string inputImagePath = inputPath + "filteredHeatEQ.raw";
	//std::string inputImagePath = inputPath + "patient1b.raw";
	//std::string inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/output/filteredP1bMC.raw";
	//std::string inputImagePath = inputPath + "patient7.raw";
	//std::string inputImagePath = inputPath + "Slices304.raw";
	//std::string inputImagePath = inputPath + "filteredMC.raw";
	//std::string inputImagePath = inputPath + "slice175.raw";

	if (load3dArrayRAW<short>(image, Length, Width, Height, inputImagePath.c_str(), false) == false)
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

	std::string loadedImagePath = outputPath + "loaded.raw";
	store3dRawData<dataType>(imageData, Length, Width, Height, loadedImagePath.c_str());

	for (k = 0; k < Height; k++) {
		delete[] image[k];
		//delete[] imageData[k];
	}
	delete[] image;
	//delete[] imageData;

	//------------------- Filtering --------------------------------

	//Filter_Parameters filterParameters; filterParameters.coef = 1e-4; filterParameters.edge_detector_coefficient = 1000; filterParameters.eps2 = 1e-4;
	//filterParameters.h = 1.0; filterParameters.maxNumberOfSolverIteration = 100; filterParameters.omega_c = 1.5; filterParameters.p = 1;
	//filterParameters.sigma = 1e-3; filterParameters.timeStepSize = 0.25; filterParameters.timeStepsNum = 1; filterParameters.tolerance = 1e-3;
	//rescaleNewRange(imageData, Length, Width, Height, 0.0, 1.0);
	//Image_Data inputImageData; inputImageData.imageDataPtr = imageData; inputImageData.height = Height; inputImageData.length = Length; inputImageData.width = Width;
	//const FilterMethod F_method = GEODESIC_MEAN_CURVATURE_FILTER;
	//filterImage(inputImageData, filterParameters, F_method);
	//std::string filteredImagePath = outputPath + "filteredGMC.raw";
	//store3dRawData<dataType>(imageData, Length, Width, Height, filteredImagePath.c_str());
	//rescaleNewRange(imageData, Length, Width, Height, 0.0, 4000.0);

	//for (k = 0; k < Height; k++) {
	//	delete[] imageData[k];
	//}
	//delete[] imageData;

	//-------------------- Resampling for to have the same spacing in all the directions -------------------

	//dataType k_spacingOld = 2.5, k_spacingNew = 1.171875;
	//const size_t HeightInterp = (size_t)((k_spacingOld / k_spacingNew) * Height);
	//dataType ** resampledImage = new dataType* [HeightInterp];
	////dataType** distanceMap = new dataType * [HeightInterp];
	//for (k = 0; k < HeightInterp; k++) {
	//	resampledImage[k] = new dataType[dim2D];
	//	//distanceMap[k] = new dataType[dim2D];
	//}
	//if (resampledImage == NULL)
	//	return false;

	//linear2dInterpolation(imageData, resampledImage, Length, Width, Height, k_spacingOld, k_spacingNew);
	//std::string interpolatedImagePath = outputPath + "interpolated.raw";
	//store3dRawData<dataType>(resampledImage, Length, Width, HeightInterp, interpolatedImagePath.c_str());

	//for (k = 0; k < Height; k++) {
	//	delete[] imageData[k];
	//}
	//delete[] imageData;

	//------------------- Fast Marching and Minimal path --------------------------------------------------
	 
	//dataType** distanceFunc = new dataType * [Height];
	//dataType** potentialFunc = new dataType * [Height];
	//dataType** resultedPath = new dataType * [Height];
	//dataType** maskThreshold = new dataType * [Height];
	//dataType** distanceMap = new dataType * [Height];
	//for (k = 0; k < Height; k++) {
	//	distanceFunc[k] = new dataType[dim2D];
	//	potentialFunc[k] = new dataType[dim2D];
	//	resultedPath[k] = new dataType[dim2D];
	//	maskThreshold[k] = new dataType[dim2D];
	//	distanceMap[k] = new dataType[dim2D];
	//}
	//if (distanceFunc == NULL || potentialFunc == NULL || resultedPath == NULL || maskThreshold == NULL || distanceMap == NULL)
	//	return false;

	//point3d* seed = new point3d[2];

	////////Patient 2
	//////seed[0].x = 261; seed[0].y = 257; seed[0].z = 151;
	//////seed[1].x = 295; seed[1].y = 317; seed[1].z = 261;

	////////Patient 1b
	//////seed[0].x = 288; seed[0].y = 308; seed[0].z = 364;
	//////seed[1].x = 259; seed[1].y = 256; seed[1].z = 244;

	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			xd = x_new(j, i, Width);
	//			distanceFunc[k][xd] = 0.0;
	//			potentialFunc[k][xd] = 0.0;
	//			resultedPath[k][xd] = 0.0;
	//			maskThreshold[k][xd] = imageData[k][xd];
	//			distanceMap[k][xd] = 0.0;

	//			/*if (sqrt(pow(seed[0].x - j, 2) + pow(seed[0].y - i, 2) + pow(seed[0].z - k, 2)) <= 5) {
	//				resultedPath[k][xd] = 1.0;
	//			}
	//			if (sqrt(pow(seed[1].x - j, 2) + pow(seed[1].y - i, 2) + pow(seed[1].z - k, 2)) <= 5) {
	//				resultedPath[k][xd] = 1.0;
	//			}*/
	//		}
	//	}
	//}

	//thresholding3dFunctionN(maskThreshold, Length, Width, Height, thresmin, thresmax, 0.0, 1.0);
	//std::string thresholded = outputPath + "threshold.raw";
	//store3dRawData<dataType>(maskThreshold, Length, Width, Height, thresholded.c_str());

	////Fast sweeping to find the point with the higest distance
	//Distance_Map_Params distParameters = { 0.4, 1.0, 0.0, 100000, 0.001 };
	//DistanceMapMethod D_method = FAST_SWEEP;
	//computeDistanceMap(distanceMap, maskThreshold, Length, Width, Height, distParameters, D_method);

	////finding of the point with the highest distance
	//dataType distanceMax = -1;
	//dataType i_max, j_max, k_max;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (distanceMap[k][x_new(i, j, Length)] >= distanceMax) {
	//				distanceMax = distanceMap[k][x_new(i, j, Length)];
	//				i_max = i; j_max = j; k_max = k;
	//			}
	//		}
	//	}
	//}
	//seed[0].x = seed[1].x = j_max;
	//seed[0].y = seed[1].y = i_max;
	//seed[1].z = seed[1].z = k_max;

	//fastMarching3D_N(imageData, distanceFunc, potentialFunc, Length, Width, Height, seed);
	//std::string distance = outputPath + "distanceMarch.raw";
	//store3dRawData<dataType>(distanceFunc, Length, Width, Height, distance.c_str());

	//distance = outputPath + "distanceSweep.raw";
	//store3dRawData<dataType>(distanceMap, Length, Width, Height, distance.c_str());

	//distance = outputPath + "potential.raw";
	//store3dRawData<dataType>(potentialFunc, Length, Width, Height, distance.c_str());

	//////shortestPath3d(distanceFunc, resultedPath, Length, Width, Height, 1.0, seed);
	//////std::string resultedImagePath = outputPath + "minimalPath.raw";
	//////store3dRawData<dataType>(resultedPath, Length, Width, Height, resultedImagePath.c_str());

	//delete[] seed;
	//for (k = 0; k < Height; k++) {
	//	delete[] imageData[k];
	//	delete[] distanceFunc[k];
	//	delete[] potentialFunc[k];
	//	delete[] resultedPath[k];
	//	delete[] distanceMap[k];
	//	delete[] maskThreshold[k];
	//}
	//delete[] imageData;
	//delete[] distanceFunc;
	//delete[] potentialFunc;
	//delete[] resultedPath;
	//delete[] distanceMap;
	//delete[] maskThreshold;

	//for (k = 0; k < Height; k++) {
	//	delete[] imageData[k];
	//}
	//delete[] imageData;
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

	//Point3D* seed = (Point3D*)malloc(sizeof(Point3D));
	//seed->x = i_max; seed->y = j_max; seed->z = k_max;
	//regionGrowing(resampledImage, segmented, status, Length, Width, HeightNew, thresmin, thresmax, seed);

	//std::string segmentedImagePath = outputPath + "segmentedThreeTimesErosion.raw";
	//store3dRawData<dataType>(segmented, Length, Width, HeightNew, segmentedImagePath.c_str());

	//------------------------------------------------------------------------------------------------------

	//Manual croppping

	//kMin depends on the interpolated image
	//size_t kMin = 332, iMin = 126, jMin = 180, kn, in, jn;
	//const size_t heightNew = 170, lengthNew = 180, widthNew = 180, dim2dNew = lengthNew * widthNew;

	////Downsampled x2
	//size_t kMin = 160, iMin = 60, jMin = 85, kn = 0, in = 0, jn = 0;
	//const size_t heightNew = 90, lengthNew = 100, widthNew = 100, dim2dNew = lengthNew * widthNew;

	//Downsampled x4
	size_t kMin = 83, iMin = 30, jMin = 40, kn, in, jn;
	const size_t heightNew = 40, lengthNew = 50, widthNew = 50, dim2dNew = lengthNew * widthNew;

	dataType** croppedImage = new dataType* [heightNew];
	dataType** maskThreshold = new dataType* [heightNew];
	dataType** distanceMap = new dataType* [heightNew];
	dataType** initialSegment = new dataType* [heightNew];
	for (k = 0; k < heightNew; k++) {
		croppedImage[k] = new dataType[dim2dNew];
		distanceMap[k] = new dataType[dim2dNew];
		maskThreshold[k] = new dataType[dim2dNew];
		initialSegment[k] = new dataType[dim2dNew];
	}
	if (croppedImage == NULL || distanceMap == NULL || maskThreshold == NULL || initialSegment == NULL) {
		return false;
	}

	for (k = 0, kn = kMin; k < heightNew; k++, kn++) {
		for (i = 0, in = iMin; i < lengthNew; i++, in++) {
			for (j = 0, jn = jMin; j < widthNew; j++, jn++) {
				x = x_new(i, j, lengthNew);
				croppedImage[k][x] = imageData[kn][x_new(in, jn, Length)];
				maskThreshold[k][x] = croppedImage[k][x];
				distanceMap[k][x] = 0.0;
				initialSegment[k][x] = 0.0;
			}
		}
	}

	std::string croppedVolumePath = outputPath + "croppedImage.raw";
	store3dRawData<dataType>(croppedImage, lengthNew, widthNew, heightNew, croppedVolumePath.c_str());

	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
	}
	delete[] imageData;

	//--------------------------Point with the highest distance-------------------------------------------------

	//dataType** distance = new dataType * [Height];
	//dataType** maskThreshold = new dataType* [Height];
	//for (k = 0; k < Height; k++) {
	//	distance[k] = new dataType[dim2D];
	//	maskThreshold[k] = new dataType[dim2D];
	//}
	//if (distance == NULL || maskThreshold == NULL)
	//	return false;

	//for (i = 0; i < Height; i++) {
	//	for (j = 0; j < Width * Length; j++) {
	//		maskThreshold[i][j] = imageData[i][j];
	//	}
	//}

	//std::string loadedImage = outputPath + "beforeThresh.raw";
	//store3dRawData<dataType>(maskThreshold, Length, Width, Height, loadedImage.c_str());

	thresholding3dFunctionN(maskThreshold, lengthNew, widthNew, heightNew, thresmin, thresmax, 0.0, 1.0);
	//thresholding3dFunctionN(imageData, Length, Width, Height, thresmin, thresmax, 0.0, 1.0);
	std::string thresholded = outputPath + "thresholded.raw";
	store3dRawData<dataType>(maskThreshold, lengthNew, widthNew, heightNew, thresholded.c_str());
	//store3dRawData<dataType>(imageData, Length, Width, Height, thresholded.c_str());
	 
	//Fast sweeping to find the point with the higest distance
	Distance_Map_Params distParameters = {0.4, 1.0, 0.0, 100000, 0.001};
	DistanceMapMethod D_method = FAST_SWEEP;
	computeDistanceMap(distanceMap, maskThreshold, lengthNew, widthNew, heightNew, distParameters, D_method);
	//computeDistanceMap(distanceMap, imageData, Length, Width, Height, distParameters, D_method);

	std::string distanceMapPath = outputPath + "distance.raw";
	store3dRawData<dataType>(distanceMap, lengthNew, widthNew, heightNew, distanceMapPath.c_str());

	//finding of the point with the highest distance
	dataType distanceMax = -1;
	dataType i_max, j_max, k_max;
	for (k = 0; k < heightNew; k++) {
		for (i = 0; i < lengthNew; i++) {
			for (j = 0; j < widthNew; j++) {
				if (distanceMap[k][x_new(i, j, lengthNew)] >= distanceMax) {
					distanceMax = distanceMap[k][x_new(i, j, lengthNew)];
					i_max = i; j_max = j; k_max = k;
				}
			}
		}
	}

	////finding of the point with the highest distance
	//dataType distanceMax = -1;
	//dataType i_max, j_max, k_max;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (distanceMap[k][x_new(i, j, Length)] >= distanceMax) {
	//				distanceMax = distanceMap[k][x_new(i, j, Length)];
	//				i_max = i; j_max = j; k_max = k;
	//			}
	//		}
	//	}
	//}
	
	//cout << "Maximal distance for the cropped volume : " << distanceMax << endl;
	//cout << "Coordinates of the highest distance x =  " << i_max << ",y =  " << j_max << ",z =  " << k_max << endl;

	////for (k = 0; k < Height; k++) {
	////	for (i = 0; i < Length; i++) {
	////		for (j = 0; j < Width; j++) {
	////			xd = x_new(i, j, Length);
	////			if (sqrt(pow(i - i_max, 2) + pow(j - j_max, 2) + pow(k - k_max, 2)) <= 5) {
	////				maskThreshold[k][xd] = 1.0;
	////			}
	////			else {
	////				maskThreshold[k][xd] = 0.0;
	////			}
	////		}
	////	}
	////}
	////std::string ballImagePath = outputPath + "ball.raw";
	////store3dRawData<dataType>(maskThreshold, lengthNew, widthNew, heightNew, ballImagePath.c_str());

	//for (k = 0; k < heightNew; k++) {
	//	//delete[] croppedImage[k]; delete[] initialSegment[k];
	//	//delete[] distanceMap[k];
	//	delete[] maskThreshold[k];
	//}
	////delete[] croppedImage; delete[] initialSegment;
	////delete[] distanceMap;
	//delete[] maskThreshold;

	//for (k = 0; k < Height; k++) {
	//	delete[] imageData[k];
	//	delete[] resampledImage[k];
	//	delete[] distanceMap[k];
	//}
	//delete[] imageData;
	//delete[] resampledImage;
	//delete[] distanceMap;

	//------------------------------Segmentation--------------------------------------------------------
	 
	//Segmentation parameters for real image
	size_t numb_centers = 1; 
	Point3D* centerSeg = new Point3D[numb_centers];
	centerSeg->x = i_max; centerSeg->y = j_max; centerSeg->z = k_max; //---> used for one center
	
	//////used for multiple centers
	////centerSeg[0].x = 40; centerSeg[0].y = 131; centerSeg[0].z = 99;
	////centerSeg[1].x = 46; centerSeg[1].y = 103; centerSeg[1].z = 99;
	////centerSeg[2].x = 50; centerSeg[2].y = 73; centerSeg[2].z = 99;
	////centerSeg[3].x = 63; centerSeg[3].y = 57; centerSeg[3].z = 99;
	////centerSeg[4].x = 113; centerSeg[4].y = 31; centerSeg[4].z = 99;
	//////If we want to start with the segmentatation function originally implemented in the library

	////generateInitialSegmentationFunctionForMultipleCentres(initialSegment, lengthNew, widthNew, heightNew, centerSeg, 0.5, 30, numb_centers);
	generateInitialSegmentationFunctionForMultipleCentres(initialSegment, lengthNew, widthNew, heightNew, centerSeg, 1.0, 10, numb_centers);
	std::string segmFolderPath = outputPath + "segmentation/";
	store3dRawData<dataType>(initialSegment, lengthNew, widthNew, heightNew, (segmFolderPath + std::string("_seg_func_000.raw")).c_str());

	Image_Data segment; segment.height = heightNew; segment.length = lengthNew; segment.width = widthNew; segment.imageDataPtr = croppedImage;
	rescaleNewRange(segment.imageDataPtr, lengthNew, widthNew, heightNew, 0.0, 1.0);

	Segmentation_Parameters segmentParameters; segmentParameters.coef = 100000; segmentParameters.eps2 = 1e-6; segmentParameters.gauss_seidelTolerance = 1e-6;
	segmentParameters.h = 1.0; segmentParameters.maxNoGSIteration = 100; segmentParameters.maxNoOfTimeSteps = 500; segmentParameters.mod = 10;
	segmentParameters.numberOfTimeStep = 500; segmentParameters.omega_c = 1.5; segmentParameters.segTolerance = 1e-10; segmentParameters.tau = 8.0;
	segmentParameters.coef_conv = 1.0; segmentParameters.coef_dif = 1.0;

	Filter_Parameters filter_Parameters; filter_Parameters.coef = 1e-2; filter_Parameters.edge_detector_coefficient = 1; filter_Parameters.eps2 = 1e-4;
	filter_Parameters.h = 1.0; filter_Parameters.maxNumberOfSolverIteration = 100; filter_Parameters.omega_c = 1.5; filter_Parameters.p = 1;
	filter_Parameters.sigma = 1e-3; filter_Parameters.timeStepSize = 0.25; filter_Parameters.timeStepsNum = 1; filter_Parameters.tolerance = 1e-3;

	unsigned char outputPathPtr[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/";
	//subsurfSegmentation(segment, initialSegment, segmentParameters, filter_Parameters, centerSeg, numb_centers, outputPathPtr);
	//generalizedSubsurfSegmentation(segment, initialSegment, segmentParameters, filter_Parameters, centerSeg, numb_centers, outputPathPtr);

	const SegmentationMethod model = SUBSURF_MODEL;
	segmentImage(segment, initialSegment, segmentParameters, filter_Parameters, centerSeg, numb_centers, outputPathPtr, model);

	//inputImagePath = outputPath + "segmentation/_seg_func_5000.raw";
	//if (load3dArrayRAW<dataType>(croppedImage, lengthNew, widthNew, heightNew, inputImagePath.c_str()) == false)
	//{
	//	printf("inputImagePath does not exist\n");
	//}

	//thresholding3dFunctionN(croppedImage, lengthNew, widthNew, heightNew, 0.3, 1.0, 1.0, 0.0);
	//computeDistanceMap(distanceMap, croppedImage, lengthNew, widthNew, heightNew, distParameters, D_method);

	//std::string thresholded = outputPath + "New_distance.raw";
	//store3dRawData<dataType>(distanceMap, lengthNew, widthNew, heightNew, thresholded.c_str());

	//New initial segmentation function
	
	//for (k = 0; k < heightNew; k++) {
	//	for (i = 0; i < lengthNew * widthNew; i++) {
	//		initialSegment[k][i] = (dataType)(1. / (distanceMap[k][i] + 1.0));
	//	}
	//}
	//rescaleNewRange(distanceMap, lengthNew, widthNew, heightNew, 0.0, 1.0);

	//for (k = 0, kn = kMin; k < heightNew; k++, kn++) {
	//	for (i = 0, in = iMin; i < lengthNew; i++, in++) {
	//		for (j = 0, jn = jMin; j < widthNew; j++, jn++) {
	//			x = x_new(i, j, lengthNew);
	//			croppedImage[k][x] = imageData[kn][x_new(in, jn, Length)];
	//		}
	//	}
	//}

	//segment.imageDataPtr = croppedImage;
	//rescaleNewRange(segment.imageDataPtr, lengthNew, widthNew, heightNew, 0, 1);
	//unsigned char seg2PathPtr[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/segmentationNew/";
	//generalizedSubsurfSegmentation(segment, initialSegment, segmentParameters, filter_Parameters, centerSeg, numb_centers, seg2PathPtr, 1.0, 1.0);

	delete[] centerSeg;
	for (k = 0; k < heightNew; k++) {
		delete[] maskThreshold[k];
		delete[] distanceMap[k];
		delete[] initialSegment[k];
		delete[] croppedImage[k];
	}
	delete[] maskThreshold;
	delete[] distanceMap;
	delete[] croppedImage;
	delete[] initialSegment;

	//for (k = 0; k < Height; k++) {
	//	delete[] imageData[k];
	//	//delete[] distanceMap[k];
	//}
	//delete[] imageData;
	////delete[] distanceMap;

	return EXIT_SUCCESS;
}
