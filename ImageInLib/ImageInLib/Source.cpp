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

//#define thresmin 995
//#define thresmax 1213


int main() {

	size_t i, j, k, x;

	//image Dimensions
	const size_t Width = 512;
	const size_t Length = 512;
	const size_t Height = 406;
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
		delete[] image[k];
	}
	delete[] image;

	load3dArrayRAW<dataType>(imageData, Length, Width, Height, inputImagePath.c_str());
	
	//std::string loaded3D = outputPath + "loaded.raw";
	//store3dRawData<dataType>(imageData, Length, Width, Height, loaded3D.c_str());

	//---------- 3D cropping -------------------------------------------

	const size_t LengthNew = 130, WidthNew = 130, HeightNew = 220;

	dataType** cropped = new dataType* [HeightNew];
	dataType** distanceMap3D = new dataType* [HeightNew];
	dataType** potential = new dataType * [HeightNew];
	dataType** path3D = new dataType * [HeightNew];
	bool** status = new bool * [HeightNew];
	for (k = 0; k < HeightNew; k++) {
		cropped[k] = new dataType [LengthNew * WidthNew];
		distanceMap3D[k] = new dataType[LengthNew * WidthNew];
		potential[k] = new dataType[LengthNew * WidthNew];
		path3D[k] = new dataType[LengthNew * WidthNew];
		status[k] = new bool[LengthNew * WidthNew];
	}
	if (cropped == NULL || distanceMap3D == NULL || potential == NULL || path3D == NULL || status == NULL)
		return false;

	for (k = 0; k < HeightNew; k++) {
		for (i = 0; i < LengthNew; i++) {
			for (j = 0; j < WidthNew; j++) {
				x = x_new(i, j, LengthNew);
				cropped[k][x] = imageData[k + 70][x_new(i + 190, j + 220, Length)];
				//cropped[k][x] = imageData[k + 150][x_new(i + 160, j + 200, Length)];
				//potential[k][x] = 0;
				//cropped[k][x] = 0; //imageData[k + 144][x_new(i + 160, j + 200, Length)];
				//cropped[k][x] = imageData[k][x_new(i, j, Length)];
				distanceMap3D[k][x] = 0;
				path3D[k][x] = 0;
				status[k][x] = false;
			}
		}
	}

	std::string cropped3D = outputPath + "croppedImage.raw";
	store3dRawData<dataType>(cropped, LengthNew, WidthNew, HeightNew, cropped3D.c_str());

	//------------ //Draw 3D spiral //---------------------------

	//dataType p = 0.19, a = 0.029, b = 0.191;
	//dataType del = atan(p / (2 * pi * b));
	//size_t n = 5;
	//u in [0,2*pi], with 10 points
	//v in [0,2*pi*n] with 300 points
	//dataType x1 = b + a * cos(u);
	//dataType x2 = -a * sin(u) * sin(del)
	//x = (x1*sin(v)) + (x2*cos(v));
	//y = (-x1*cos(v)) + (x2*sin(v));
	//z = ((p * v) / (2 * pi)) + (a * sin(u) * cos(del));

	//Point3D* seedP = (Point3D*)malloc(2 * sizeof(Point3D));
	//seedP[0].x = 77; seedP[0].y = 55, seedP[0].z = 0;
	//seedP[1].x = 110; seedP[1].y = 116, seedP[1].z = 90;
	
	////Point3D centerBall; centerBall.x = 50; centerBall.y = 50; centerBall.z = 50;
	////unsigned char ball3D[] = "C:/Users/Konan Allaly/Documents/Tests/output/ball.vtk";
	//generateSphere(cropped, centerBall, Length, Width, Height, 10, 0.0, ball3D);

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

	//cropped[20][x_new(50, 50, LengthNew)] = 1;
	//////thresholding3dFunctionN(cropped, Length, Width, Height, thresmin, thresmax, 0, 1);

	Point3D* seedP = (Point3D*)malloc(sizeof(Point3D));
	seedP->x = 100; seedP->y = 87; seedP->z = 198;

	Distance_Map_Params distanceParamerters; distanceParamerters.h = 1.0; distanceParamerters.initValue = 1000000000;
	distanceParamerters.tau = 0.4; distanceParamerters.tolerance = 0.5; distanceParamerters.objectPixel = 1;
	DistanceMapMethod distMethod = FAST_MARCH; //FAST_SWEEP; //BRUTE_FORCE; //ROUY_TOURIN;       
	computeDistanceMap(distanceMap3D, cropped, LengthNew, WidthNew, HeightNew, distanceParamerters, distMethod, seedP);
	std::string ImagePath3D = outputPath + "distance.raw";
	store3dRawData<dataType>(distanceMap3D, LengthNew, WidthNew, HeightNew, ImagePath3D.c_str());

	//Freeing imageData pointer 
	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
	}
	delete[] imageData;
	delete[] cropped;
	delete[] distanceMap3D;
	delete[] potential;
	delete[] path3D;
	//free(seedP);
 
	//---------------- Fast Marching new implementation ---------------------------------

	//Point3d* seedP = (Point3d*)malloc(2 * sizeof(Point3d));
	//seedP[1].x = 121; seedP[1].y = 102, seedP[1].z = 169;
	//seedP[0].x = 92; seedP[0].y = 58, seedP[0].z = 45;

	//Point3D* seedPoint = (Point3D*)malloc(sizeof(Point3D));
	//seedPoint->x = 288; seedPoint->y = 158; seedPoint->z = 163;

	//compute3dPotential(cropped, potential, LengthNew, WidthNew, HeightNew, seedP);
	//std::string potential3D = outputPath + "potential.raw";
	//store3dRawData<dataType>(potential, LengthNew, WidthNew, HeightNew, potential3D.c_str());

	//cout << "Value fisrt point : " << cropped[seedP[0].z][x_new(seedP[0].x, seedP[0].y, WidthNew)] << endl;;
	//cout << "\nValue second point : " << cropped[seedP[1].z][x_new(seedP[1].x, seedP[1].y, WidthNew)] << endl;
	//cout << "Difference potential " << potential[seedP[0].z][x_new(seedP[0].x, seedP[0].y, WidthNew)] - potential[seedP[1].z][x_new(seedP[1].x, seedP[1].y, WidthNew)];

	//dataType thres_min = 0.02; //potential[seedP[0].z][x_new(seedP[0].x, seedP[0].y, WidthNew)];
	//dataType thres_max = 0.05; //2 * thres_min;

	//regionGrowing(potential, distanceMap3D, status, LengthNew, Width, HeightNew, thres_min, thres_max, seedPoint);
	//std::string segmented = outputPath + "segmentedImage.raw";
	//store3dRawData<dataType>(distanceMap3D, LengthNew, WidthNew, HeightNew, segmented.c_str());

	//fastMarching3d_N(cropped, distanceMap3D, potential, LengthNew, WidthNew, HeightNew, seedP);
	//std::string DistancePath = outputPath + "distanceRealSpeed.raw";
	//store3dRawData<dataType>(distanceMap3D, LengthNew, WidthNew, HeightNew, DistancePath.c_str());

	//shortestPath3d(distanceMap3D, path3D, LengthNew, WidthNew, HeightNew, 1.0, seedP);
	//std::string ImagePath3D = outputPath + "resultPath.raw";
	//store3dRawData<dataType>(path3D, LengthNew, WidthNew, HeightNew, ImagePath3D.c_str());

	////free(seedP);
	//for (k = 0; k < HeightNew; k++) {
	//	delete[] cropped[k];
	//	delete[] distanceMap3D[k];
	//	delete[] potential[k];
	//	delete[] path3D[k];
	//}
	//delete[] cropped;
	//delete[] distanceMap3D;
	//delete[] potential;
	//delete[] path3D;

	//---------------Real 2D image -----------------------

	////Extract Slice of interest

	//dataType* real2dImageData = new dataType [Height * Width];
	//dataType* loadingPtr = new dataType[Height * Width];

	//i = 0;
	////size_t cst = 238; // test 1 and 2
	//size_t cst = 234; // test 3
	////size_t cst = 285; // test 4
	////size_t cst = 280; // test 5
	//for (k = 0; k < Height; k++) {
	//	for (j = 0; j < Width; j++) {
	//		x = x_new(cst, j, Length);
	//		loadingPtr[i] = imageData[k][x];
	//		i++;
	//	}
	//}

	//for (i = 0; i < Height; i++) {
	//	for (j = 0; j < Width; j++) {
	//		real2dImageData[x_new(i, j, Height)] = loadingPtr[x_new(Height - i - 1, Width - j - 1, Height)];
	//	}
	//}

	//std::string real2dImagePath = outputPath + "loaded2dImage.raw";
	//store2dRawData<dataType>(real2dImageData, Height, Width, real2dImagePath.c_str());

	//delete[] loadingPtr;
	//for (k = 0; k < Height; k++) {
	//	delete[] imageData[k];
	//}
	//delete[] imageData;

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


	//-----------Create artificial image--------------------

	//const size_t Height = 200;
	//const size_t Width = 200;
	//const size_t dim2D = Width * Height;
	//size_t i, j, x;
	// 
	//dataType* artificial2dImage = new dataType[dim2D];
	//std::string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";

	////Draw snake
	//for (i = 0; i < Height; i++) {
	//	for (j = 0; j < Width; j++) {
	//		artificial2dImage[x_new(j, i, Width)] = 0;
	//	}
	//}
	//size_t i1 = 50, j1 = 100, radius = 50;
	//for (i = 0; i < 100; i++) {
	//	for (j = j1; j < ( j1 + radius) ; j++) {
	//		if (sqrt((i1 - i) * (i1 - i) + (j1 - j) * (j1 - j)) <= radius) {
	//			artificial2dImage[x_new(j, i, Width)] = 1;
	//		}
	//		if (sqrt((i1 - i) * (i1 - i) + (j1 - j) * (j1 - j)) <= 20) {
	//			artificial2dImage[x_new(j, i, Width)] = 0;
	//		}
	//	}
	//}	 
	//size_t i2 = 120;
	//for (i = 70; i < 170; i++) {
	//	for (j = 50; j < 100; j++) {
	//		if (sqrt((i2 - i) * (i2 - i) + (j1 - j) * (j1 - j)) <= radius) {
	//			artificial2dImage[x_new(j, i, Width)] = 1;
	//		}
	//		if (sqrt((i2 - i) * (i2 - i) + (j1 - j) * (j1 - j)) <= 20) {
	//			artificial2dImage[x_new(j, i, Width)] = 0;
	//		}
	//	}
	//}
	//Point2D* startAndEnd = (Point2D*)malloc(2 * sizeof(Point2D));
	//startAndEnd[1].x = 103; startAndEnd[1].y = 16;
	//startAndEnd[0].x = 95; startAndEnd[0].y = 157;

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
	//
	//std::string pathArtificialImage = outputPath + "artificial.raw";
	//store2dRawData<dataType>(artificial2dImage, Height, Width, pathArtificialImage.c_str());

	////Draw maze
	//for (i = 0; i < Height; i++) {
	//	for (j = 0; j < Width; j++) {
	//		x = x_new(j, i, Width);
	//	}
	//}

	//----------- 2D Fast Marching ------------------------

	//dataType* distanceMap = new dataType[Height * Width];
	//dataType* potentialPtr = new dataType[Height * Width];
	//dataType* pathPtr = new dataType[Height * Width];

	//////for (k = 0; k < Height; k++) {
	//////	for (i = 0; i < Length; i++) {
	//////		for (j = 0; j < Width; j++) {
	//////			x = x_new(j, i, Width);
	//////			real2dImageData[x] = imageData[k][x];
	//////		}
	//////	}
	//////}

	//for (i = 0; i < Height; i++) {
	//	for (j = 0; j < Width; j++) {
	//		x = x_new(j, i, Width);
	//		potentialPtr[x] = 0;
	//		distanceMap[x] = 0;
	//		pathPtr[x] = 0;
	//	}
	//}

	//fastMarching2d(real2dImageData, distanceMap, potentialPtr, Height, Width, startAndEnd);
	////fastMarching2d(artificial2dImage, distanceMap, potentialPtr, Height, Width, startAndEnd);
	//std::string distanceOutPut = outputPath + "distance.raw";
	//store2dRawData<dataType>(distanceMap, Height, Width, distanceOutPut.c_str());
	//distanceOutPut = outputPath + "potential.raw";
	//store2dRawData<dataType>(potentialPtr, Height, Width, distanceOutPut.c_str());

	//shortestPath2d(distanceMap, pathPtr, Height, Width, 1.0, startAndEnd);
	//distanceOutPut = outputPath + "path.raw";
	//store2dRawData<dataType>(pathPtr, Height, Width, distanceOutPut.c_str());

	//for (i = 0; i < Height; i++) {
	//	for (j = 0; j < Width; j++) {
	//		x = x_new(j, i, Width);
	//		if (pathPtr[x] == 1) {
	//			real2dImageData[x] = 0;
	//			//artificial2dImage[x] = 0;
	//		}
	//	}
	//}

	//std::string outPath = outputPath + "ImagePlusPath.raw";
	//store2dRawData<dataType>(real2dImageData, Height, Width, outPath.c_str());
	////store2dRawData<dataType>(artificial2dImage, Height, Width, outPath.c_str());

	////outPath = outputPath + "DistancePlusPath.raw";
	////store2dRawData<dataType>(distanceMap, Height, Width, outPath.c_str());

	//free(startAndEnd);
	//delete[] real2dImageData;
	////delete[] artificial2dImage;
	//delete[] distanceMap;
	//delete[] potentialPtr;
	//delete[] pathPtr;

	return EXIT_SUCCESS;
}
