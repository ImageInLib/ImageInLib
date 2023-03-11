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

#define thresmin 995
#define thresmax 1213

#define epsilonFast 0.001


int main() {

	size_t i, j, k, x;

	//image Dimensions
	const size_t Width = 512;
	const size_t Length = 512;
	const size_t Height = 406;
	const size_t dim2D = Width * Length;

	//-------------Real 3D image -------------------------
	
	//Preparation of image data pointers
	dataType** imageData = (dataType**)malloc(Height * sizeof(dataType*));
	short** image = (short**)malloc(Height * sizeof(short*));
	for (k = 0; k < Height; k++) {
		imageData[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		image[k] = (short*)malloc(dim2D * sizeof(short));
	}
	if (imageData == NULL || image == NULL) {
		return false;
	}
	std::string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/";
	std::string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";
	//std::string inputImagePath = inputPath + "patient2.raw";
	//std::string inputImagePath = inputPath + "Slices304.raw";
	////std::string inputImagePath = inputPath + "patient2_filtered.raw";
	////std::string inputImagePath = inputPath + "filteredK100.raw";
	////std::string inputImagePath = inputPath + "filteredK1.raw"; filteredHeatEQ
	//std::string inputImagePath = inputPath + "filteredHeatEQ.raw";
	std::string inputImagePath = inputPath + "filteredMC.raw";

	//if (load3dArrayRAW<short>(image, Length, Width, Height, inputImagePath.c_str()) == false)
	//{
	//	printf("inputImagePath does not exist\n");
	//}
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			x = x_new(i, j, Length);
	//			imageData[k][x] = (dataType)image[k][x];
	//		}
	//	}
	//}

	//Freeing image pointer 
	for(k = 0; k < Height; k++) {
		free(image[k]);
	}
	free(image);

	load3dArrayRAW<dataType>(imageData, Length, Width, Height, inputImagePath.c_str());

	//---------- 3D cropping -------------------------------------------
	//const size_t LengthNew = 150, WidthNew = 150, HeightNew = 100;
	const size_t LengthNew = 200, WidthNew = 200, HeightNew = 120;
	dataType** cropped = (dataType**)malloc(HeightNew * sizeof(dataType*));
	dataType** distanceMap3D = (dataType**)malloc(HeightNew * sizeof(dataType*));
	for (k = 0; k < HeightNew; k++) {
		cropped[k] = (dataType*)malloc(LengthNew * WidthNew * sizeof(dataType*));
		distanceMap3D[k] = (dataType*)malloc(LengthNew * WidthNew * sizeof(dataType*));
	}
	if (cropped == NULL || distanceMap3D == NULL)
		return false;

	for (k = 0; k < HeightNew; k++) {
		for (i = 0; i < LengthNew; i++) {
			for (j = 0; j < WidthNew; j++) {
				x = x_new(i, j, LengthNew);
				//cropped[k][x] = imageData[k + 190][x_new(i + 180, j + 230, Length)];
				cropped[k][x] = 0; //imageData[k + 144][x_new(i + 160, j + 200, Length)];
				//cropped[k][x] = imageData[k][x_new(i, j, Length)];
				distanceMap3D[k][x] = 0;
			}
		}
	}

	//std::string cropped3D = outputPath + "croppedImage.raw";
	//store3dRawData<dataType>(cropped, LengthNew, WidthNew, HeightNew, cropped3D.c_str());
	
	////Point3D centerBall; centerBall.x = 50; centerBall.y = 50; centerBall.z = 50;
	////unsigned char ball3D[] = "C:/Users/Konan Allaly/Documents/Tests/output/ball.vtk";
	////generateSphere(cropped, centerBall, Length, Width, Height, 10, 0.0, ball3D);

	//Freeing imageData pointer 
	for (k = 0; k < Height; k++) {
		free(imageData[k]);
	}
	free(imageData);

	//----- Fast Marching and Path finding for 3D image ----------------
	
	dataType** potential3D = (dataType**)malloc(HeightNew * sizeof(dataType*));
	dataType** path3D = (dataType**)malloc(HeightNew * sizeof(dataType*));
	dataType** label3D = (dataType**)malloc(HeightNew * sizeof(dataType*));
	for (k = 0; k < HeightNew; k++) {
		potential3D[k] = (dataType*)malloc(LengthNew * WidthNew * sizeof(dataType));
		path3D[k] = (dataType*)malloc(LengthNew * WidthNew * sizeof(dataType));
		label3D[k] = (dataType*)malloc(LengthNew * WidthNew * sizeof(dataType*));
	}
	if (potential3D == NULL || path3D == NULL || label3D == NULL) {
		return false;
	}

	//Initialization
	for (k = 0; k < HeightNew; k++) {
		for (i = 0; i < LengthNew; i++) {
			for (j = 0; j < WidthNew; j++) {
				x = x_new(i, j, LengthNew);
				potential3D[k][x] = 0;
				path3D[k][x] = 0;
				label3D[k][x] = 3;
			}
		}
	}

	//std::string distance3D = outputPath + "loaded.raw";
	//store3dRawData<dataType>(cropped, LengthNew, WidthNew, HeightNew, distance3D.c_str());

	Point3d* seed = (Point3d*)malloc(2 * sizeof(Point3d));
	seed[0].x = 102; seed[0].y = 61, seed[0].z = 0;
	seed[1].x = 132; seed[1].y = 113, seed[1].z = 119;

	//std::cout << "Start : " << cropped[0][x_new(103, 59, WidthNew)] << std::endl;
	//std::cout << "\nEnd : " << cropped[119][x_new(132, 117, WidthNew)] << std::endl;

	//dataType** gradientX = (dataType**)malloc(HeightNew * sizeof(dataType*));
	//dataType** gradientY = (dataType**)malloc(HeightNew * sizeof(dataType*));
	//dataType** gradientZ = (dataType**)malloc(HeightNew * sizeof(dataType*));
	//for (k = 0; k < HeightNew; k++) {
	//	gradientX[k] = (dataType*)malloc(LengthNew * WidthNew * sizeof(dataType));
	//	gradientY[k] = (dataType*)malloc(LengthNew * WidthNew * sizeof(dataType));
	//	gradientZ[k] = (dataType*)malloc(LengthNew * WidthNew * sizeof(dataType));
	//}

	//compute3dImageGradient(cropped, gradientX, gradientY, gradientZ, LengthNew, WidthNew, HeightNew, 1.0);
	//for (k = 0; k < HeightNew; k++) {
	//	for (i = 0; i < LengthNew; i++) {
	//		for (j = 0; j < WidthNew; j++) {
	//			x = x_new(i, j, LengthNew);
	//			path3D[k][x] = gradientX[k][x] * gradientX[k][x] + gradientY[k][x] * gradientY[k][x] + gradientZ[k][x] * gradientZ[k][x];
	//		}
	//	}
	//}
	//std::string gradientPath = outputPath + "normOfGradient.raw";
	//store3dRawData<dataType>(path3D, LengthNew, WidthNew, HeightNew, gradientPath.c_str());
	//store3dRawData<dataType>(gradientX, LengthNew, WidthNew, HeightNew, gradientPath.c_str());
	//store3dRawData<dataType>(gradientY, LengthNew, WidthNew, HeightNew, gradientPath.c_str());
	//store3dRawData<dataType>(gradientZ, LengthNew, WidthNew, HeightNew, gradientPath.c_str());

	//compute3dPotential(cropped, potential3D, LengthNew, WidthNew, HeightNew, seed);
	//std::string loadedImagePath = outputPath + "potential.raw";
	//store3dRawData<dataType>(potential3D, LengthNew, WidthNew, HeightNew, loadedImagePath.c_str());

	//fastMarching3d_N(cropped, distanceMap3D, potential3D, LengthNew, WidthNew, HeightNew, seed);
	//std::string distance3D = outputPath + "distanceMap.raw";
	//store3dRawData<dataType>(distanceMap3D, LengthNew, WidthNew, HeightNew, distance3D.c_str());
	//distance3D = outputPath + "potential.raw";
	//store3dRawData<dataType>(potential3D, LengthNew, WidthNew, HeightNew, distance3D.c_str());

	//shortestPath3d(distanceMap3D, path3D, LengthNew, WidthNew, HeightNew, 1.0, seed);
	//std::string resultPath3D = outputPath + "path.raw";
	//store3dRawData<dataType>(path3D, LengthNew, WidthNew, HeightNew, resultPath3D.c_str());

	free(seed);

	//free(gradientX);
	//free(gradientY);
	//free(gradientZ);


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
	//std::string filteredImagePath = outputPath + "filteredMC.raw";
	//store3dRawData<dataType>(imageData, Length, Width, Height, filteredImagePath.c_str());

	//---------------Fast marching already implemented in the library -------------

	cropped[50][1000] = 1;
	////thresholding3dFunctionN(cropped, Length, Width, Height, thresmin, thresmax, 0, 1);
	Distance_Map_Params distanceParamerters; distanceParamerters.h = 1.0; distanceParamerters.initValue = 1000000000;
	distanceParamerters.tau = 0.4; distanceParamerters.tolerance = 0.5; distanceParamerters.objectPixel = 1;
	DistanceMapMethod distMethod = FAST_MARCH; //FAST_SWEEP; //BRUTE_FORCE; //ROUY_TOURIN;       
	computeDistanceMap(distanceMap3D, cropped, LengthNew, WidthNew, HeightNew, distanceParamerters, distMethod);
	std::string ImagePath3D = outputPath + "distanceMap3DN.raw";
	store3dRawData<dataType>(distanceMap3D, LengthNew, WidthNew, HeightNew, ImagePath3D.c_str());

	for (k = 0; k < HeightNew; k++) {
		free(cropped[k]);
		free(distanceMap3D[k]);
		free(path3D[k]);
		free(label3D[k]);
		//free(gradientX[k]);
		//free(gradientY[k]);
		//free(gradientZ[k]);
	}
	free(cropped);
	free(distanceMap3D);
	free(path3D);
	free(label3D);

	////Freeing pointers 
	//for (k = 0; k < Height; k++) {
	//	free(cropped[k]);
	//	free(distanceMap3D[k]);
	//}
	//free(cropped);
	//free(distanceMap3D);

	//---------------Real 2D image -----------------------

	//Extract Slice of interest

	//dataType* real2dImageData = (dataType*)malloc(dim2D * sizeof(dataType));
	//dataType* loadingPtr = (dataType*)malloc(dim2D * sizeof(dataType));
	//dataType* maskThresh = (dataType*)malloc(dim2D * sizeof(dataType));

	//i = 0;
	//size_t cst = 238; // test 1 and 2
	//size_t cst = 234; // test 3
	//size_t cst = 285; // test 4
	////size_t cst = 280; // test 5
	//for (k = 0; k < Height; k++) {
	//	for (j = 0; j < Width; j++) {
	//		x = x_new(cst, j, Length);
	//		loadingPtr[i] = imageData[k][x];
	//		i++;
	//	}
	//}

	////Freeing imageData pointer 
	//for (k = 0; k < Height; k++) {
	//	free(imageData[k]);
	//}
	//free(imageData);

	//for (i = 0; i < Height; i++) {
	//	for (j = 0; j < Width; j++) {
	//		real2dImageData[x_new(i, j, Height)] = loadingPtr[x_new(Height - i - 1, Width - j - 1, Height)];
	//	}
	//}

	//std::string real2dImagePath = outputPath + "loaded2dImage.raw";
	//store2dRawData<dataType>(real2dImageData, Height, Width, real2dImagePath.c_str());

	//Point2D* startAndEnd = (Point2D*)malloc(2 * sizeof(Point2D));
	 
	////test 1
	//startAndEnd[0].x = 292; startAndEnd[0].y = 336;
	//startAndEnd[1].x = 333; startAndEnd[1].y = 192;

	////test 2
	//startAndEnd[0].x = 246; startAndEnd[0].y = 277;
	//startAndEnd[1].x = 250; startAndEnd[1].y = 135;

	////test 3
	//startAndEnd[0].x = 205; startAndEnd[0].y = 296;
	//startAndEnd[1].x = 240; startAndEnd[1].y = 188;

	////test 4
	//startAndEnd[1].x = 224; startAndEnd[1].y = 194;
	//startAndEnd[0].x = 225 /*213*/; startAndEnd[0].y = 130 /*133*/;
	////startAndEnd[0].x = 213; startAndEnd[0].y = 133;

	////test 5
	//startAndEnd[0].x = 231; startAndEnd[0].y = 203;
	//startAndEnd[1].x = 205; startAndEnd[1].y = 140;

	//printf("\nI(p0) = %f", real2dImageData[x_new(224, 194, Width)]);
	//printf("\nI(p1) = %f", real2dImageData[x_new(225, 130, Width)]);
	////printf("\nI(p') = %f", real2dImageData[x_new(213, 133, Width)]);

	////Manual 2D rescalling
	//dataType diffOld = 4000, diffNew = 1;
	//dataType scale_factor = (diffNew) / (diffOld);
	//for (i = 0; i < Height; i++) {
	//	for (j = 0; j < Width; j++) {
	//		x = x_new(j, i, Width);
	//		real2dImageData[x] = scale_factor * (real2dImageData[x] - 4000) + 1;
	//	}
	//}

	//---------------Threshold------------------------------

	////Manual thresholding
	//for (k = 0; k < Height; k++) {
	//	for (j = 0; j < Width; j++) {
	//		if (real2dImageData[x_new(k, j, Height)] >= thresmin && real2dImageData[x_new(k, j, Height)] <= thresmax) {
	//			maskThresh[x_new(k, j, Height)] = 1;
	//		}
	//		else {
	//			maskThresh[x_new(k, j, Height)] = 0;
	//		}
	//	}
	//}

	//-----------Create artificial image--------------------

	//const size_t Height = 200;
	//const size_t Width = 200;
	//const size_t dim2D = Width * Height;
	//size_t i, j, x;
	// 
	//dataType* artificial2dImage = (dataType*)malloc(dim2D * sizeof(dataType));
	//std::string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";

	//////Draw snake
	////for (i = 0; i < Height; i++) {
	////	for (j = 0; j < Width; j++) {
	////		artificial2dImage[x_new(j, i, Width)] = 0;
	////	}
	////}
	////size_t i1 = 50, j1 = 100, radius = 50;
	////for (i = 0; i < 100; i++) {
	////	for (j = j1; j < ( j1 + radius) ; j++) {
	////		if (sqrt((i1 - i) * (i1 - i) + (j1 - j) * (j1 - j)) <= radius) {
	////			artificial2dImage[x_new(j, i, Width)] = 1;
	////		}
	////		if (sqrt((i1 - i) * (i1 - i) + (j1 - j) * (j1 - j)) <= 20) {
	////			artificial2dImage[x_new(j, i, Width)] = 0;
	////		}
	////	}
	////}
	////size_t i2 = 120;
	////for (i = 70; i < 170; i++) {
	////	for (j = 50; j < 100; j++) {
	////		if (sqrt((i2 - i) * (i2 - i) + (j1 - j) * (j1 - j)) <= radius) {
	////			artificial2dImage[x_new(j, i, Width)] = 1;
	////		}
	////		if (sqrt((i2 - i) * (i2 - i) + (j1 - j) * (j1 - j)) <= 20) {
	////			artificial2dImage[x_new(j, i, Width)] = 0;
	////		}
	////	}
	////}
	////Point2D* startAndEnd = (Point2D*)malloc(2 * sizeof(Point2D));
	////startAndEnd[1].x = 113; startAndEnd[1].y = 18;
	////startAndEnd[0].x = 79; startAndEnd[0].y = 148;

	////Draw L
	//for (i = 0; i < Height; i++) {
	//	for (j = 0; j < Width; j++) {
	//		x = x_new(j, i, Width);
	//		artificial2dImage[x] = 0;
	//	}
	//}
	//for (i = 0; i < Height; i++) {
	//	for (j = 0; j < Width; j++) {
	//		x = x_new(j, i, Width);
	//		if ((j >= 10 && j <= 40) && (i >= 10 && i <= 190)) {
	//			artificial2dImage[x] = 1;
	//		}
	//		if ((j >= 10 && j <= 190) && (i >= 160 && i <= 190)) {
	//			artificial2dImage[x] = 1;
	//		}
	//	}
	//}
	//Point2D* startAndEnd = (Point2D*)malloc(2 * sizeof(Point2D));
	//startAndEnd[1].x = 25; startAndEnd[1].y = 15;
	//startAndEnd[0].x = 181; startAndEnd[0].y = 175;

	////Draw a cercle
	//for (i = 0; i < Height; i++) {
	//	for (j = 0; j < Width; j++) {
	//		artificial2dImage[x_new(j, i, Width)] = 0;
	//	}
	//}
	//for (i = 0; i < Height; i++) {
	//	for (j = 0; j < Width; j++) {
	//		x = x_new(j, i, Width);
	//		if (sqrt((i - 100) * (i - 100) + (j - 100) * (j - 100)) <= 70) {
	//			artificial2dImage[x] = 1;
	//		}
	//		if (sqrt((i - 100) * (i - 100) + (j - 100) * (j - 100)) <= 50) {
	//			artificial2dImage[x] = 0;
	//		}
	//	}
	//}
	//Point2D* startAndEnd = (Point2D*)malloc(2 * sizeof(Point2D));
	//startAndEnd[0].x = 40; startAndEnd[0].y = 88;
	//startAndEnd[1].x = 155 /*144*/; /*106*/  startAndEnd[1].y = 116 /*58*/; /*158*/
	
	//std::string pathArtificialImage = outputPath + "artificial.raw";
	//store2dRawData<dataType>(artificial2dImage, Height, Width, pathArtificialImage.c_str());

	////Draw maze
	//for (i = 0; i < Height; i++) {
	//	for (j = 0; j < Width; j++) {
	//		x = x_new(j, i, Width);
	//	}
	//}

	//------------------------------------------------------

	//dataType* distanceMap = (dataType*)malloc(dim2D * sizeof(dataType));
	//dataType* potentialPtr = (dataType*)malloc(dim2D * sizeof(dataType));
	//dataType* pathPtr = (dataType*)malloc(dim2D * sizeof(dataType));
	//dataType* gradientX = (dataType*)malloc(dim2D * sizeof(dataType));
	//dataType* gradientY = (dataType*)malloc(dim2D * sizeof(dataType));
	//dataType* gradientNorm = (dataType*)malloc(dim2D * sizeof(dataType));

	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			x = x_new(j, i, Width);
	//			real2dImageData[x] = imageData[k][x];
	//		}
	//	}
	//}

	//for (i = 0; i < Height; i++) {
	//	for (j = 0; j < Width; j++) {
	//		x = x_new(j, i, Width);
	//		potentialPtr[x] = 0;
	//		distanceMap[x] = 0;
	//		pathPtr[x] = 0;
	//		gradientX[x] = 0;
	//		gradientY[x] = 0;
	//	}
	//}

	//fastMarching2d(real2dImageData, distanceMap, potentialPtr, Length, Width, startAndEnd);

	//fastMarching2d(artificial2dImage, distanceMap, potentialPtr, Height, Width, startAndEnd);
	//fastMarching2d(maskThresh, distanceMap, potentialPtr, Height, Width, startAndEnd);

	////computeImageGradient(real2dImageData, gradientX, gradientY, Height, Width, 1.0);
	//for (i = 0; i < Height; i++) {
	//	for (j = 0; j < Width; j++) {
	//		x = x_new(j, i, Width);
	//		gradientNorm[x] = sqrt(gradientX[x] * gradientX[x] + gradientY[x] * gradientY[x]);
	//	}
	//}
	//std::string gradientOutPut = outputPath + "normOfGradient.raw";
	////store2dRawData<dataType>(gradientNorm, Height, Width, gradientOutPut.c_str());

	//std::string distanceOutPut = outputPath + "distance.raw";
	////store2dRawData<dataType>(distanceMap, Height, Width, distanceOutPut.c_str());
	//distanceOutPut = outputPath + "potential.raw";
	////store2dRawData<dataType>(potentialPtr, Height, Width, distanceOutPut.c_str());

	//shortestPath2d(distanceMap, pathPtr, Height, Width, 1.0, startAndEnd);

	//for (i = 0; i < Height; i++) {
	//	for (j = 0; j < Width; j++) {
	//		x = x_new(j, i, Width);
	//		if (pathPtr[x] == 1) {
	//			distanceMap[x] = 0;
	//			real2dImageData[x] = 0;
	//			//artificial2dImage[x] = 0;
	//		}
	//	}
	//}

	//std::string outPath = outputPath + "ImagePlusPath.raw";
	//store2dRawData<dataType>(real2dImageData, Height, Width, outPath.c_str());

	//store2dRawData<dataType>(artificial2dImage, Height, Width, outPath.c_str());

	//outPath = outputPath + "DistancePlusPath.raw";
	//store2dRawData<dataType>(distanceMap, Height, Width, outPath.c_str());

	////========================Add start and Point =================
	//distanceMap[x_new(startPoint->x, startPoint->y - 1, Width)] = 0;
	//distanceMap[x_new(startPoint->x, startPoint->y - 1, Width)] = 0;
	//distanceMap[x_new(startPoint->x, startPoint->y + 1, Width)] = 0;
	//distanceMap[x_new(startPoint->x - 1, startPoint->y, Width)] = 0;
	//distanceMap[x_new(startPoint->x - 1, startPoint->y - 1, Width)] = 0;
	//distanceMap[x_new(startPoint->x - 1, startPoint->y + 1, Width)] = 0;
	//distanceMap[x_new(startPoint->x + 1, startPoint->y, Width)] = 0;
	//distanceMap[x_new(startPoint->x + 1, startPoint->y - 1, Width)] = 0;
	//distanceMap[x_new(startPoint->x + 1, startPoint->y + 1, Width)] = 0;
	//distanceMap[x_new(endPoint->x, endPoint->y, Width)] = 0;
	//distanceMap[x_new(endPoint->x, endPoint->y - 1, Width)] = 0;
	//distanceMap[x_new(endPoint->x, endPoint->y + 1, Width)] = 0;
	//distanceMap[x_new(endPoint->x - 1, endPoint->y - 1, Width)] = 0;
	//distanceMap[x_new(endPoint->x - 1, endPoint->y + 1, Width)] = 0;
	//distanceMap[x_new(endPoint->x - 1, endPoint->y, Width)] = 0;
	//distanceMap[x_new(endPoint->x + 1, endPoint->y - 1, Width)] = 0;
	//distanceMap[x_new(endPoint->x + 1, endPoint->y + 1, Width)] = 0;
	//distanceMap[x_new(endPoint->x + 1, endPoint->y, Width)] = 0;
	////=============================================================

	////free memory
	//for (k = 0; k < Height; k++) {
	//	free(imageData[k]); free(image[k]);
	//}
	//free(imageData); free(image);

	//for (k = 0; k < Height; k++) {
	//	free(potential3D[k]);
	//	free(distanceMap3D[k]);
	//	free(path3D[k]);
	//	free(cropped[k]);
	//}
	//free(potential3D); free(distanceMap3D); free(path3D);
	//free(cropped);
	////free(seed);

	//free(real2dImageData); free(loadingPtr); free(maskThresh);

	////free(artificial2dImage);
	//free(distanceMap); free(potentialPtr); free(pathPtr);
	//free(gradientX); free(gradientY); free(gradientNorm);
	//free(startAndEnd);

	//free(seed);

	//for (k = 0; k < Height; k++) {
	//	free(gradientX[k]);
	//	free(gradientY[k]);
	//	free(gradientZ[k]);
	//}
	//free(gradientX); free(gradientY); free(gradientZ);

	return EXIT_SUCCESS;
}
