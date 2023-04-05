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

	size_t i, j, k, x, xd, m, n;

	//image Dimensions
	const size_t Width = 512;
	const size_t Length = 512;
	const size_t Height = 406; /*607*/ /*508*/
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
	//std::string inputImagePath = inputPath + "patient1b.raw";
	//std::string inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/output/filteredP1bMC.raw";
	//std::string inputImagePath = inputPath + "patient7.raw";
	//std::string inputImagePath = inputPath + "Slices304.raw";

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

	load3dArrayRAW<dataType>(imageData, Length, Width, Height, inputImagePath.c_str());
	
	//std::string loaded3D = outputPath + "loaded.raw";
	//store3dRawData<dataType>(imageData, Length, Width, Height, loaded3D.c_str());

	//Freeing image pointer 
	for (k = 0; k < Height; k++) {
		delete[] image[k];
	}
	delete[] image;

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

	//---------------Fast marching already implemented in the library -------------

	const size_t LengthNew = 130, WidthNew = 130, HeightNew = 220;
	dataType** cropped = new dataType * [HeightNew];
	dataType** distanceMap3D = new dataType * [HeightNew];
	dataType** potential = new dataType * [HeightNew];
	dataType** path3D = new dataType * [HeightNew];
	bool** status = new bool* [HeightNew];
	for (k = 0; k < HeightNew; k++) {
		cropped[k] = new dataType[LengthNew * WidthNew];
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
				//cropped[k][x] = imageData[k + 150][x_new(i + 160, j + 200, Length)];
				//cropped[k][x] = 0; //imageData[k + 144][x_new(i + 160, j + 200, Length)];
				//cropped[k][x] = imageData[k][x_new(i, j, Length)];
				cropped[k][x] = 0; // imageData[k + 70][x_new(i + 190, j + 220, Length)];
				potential[k][x] = 0;
				distanceMap3D[k][x] = 0;
				distanceMap3D[k][x] = 0;
				path3D[k][x] = 0;
				status[k][x] = false;
			}
		}
	}

	cropped[0][0] = 1;
	//thresholding3dFunctionN(cropped, LengthNew, WidthNew, HeightNew, 1050, 1200, 0, 1);

	Point3D* seedPoint = (Point3D*)malloc(2 * sizeof(Point3D));
	seedPoint[0].x = 100; seedPoint[0].y = 87; seedPoint[0].z = 198;
	seedPoint[1].x = 70; seedPoint[1].y = 37; seedPoint[1].z = 77;

	Distance_Map_Params distanceParamerters; distanceParamerters.h = 1.0; distanceParamerters.initValue = 1000000000;
	distanceParamerters.tau = 0.4; distanceParamerters.tolerance = 0.5; distanceParamerters.objectPixel = 1;
	DistanceMapMethod distMethod = FAST_MARCH;  // ROUY_TOURIN; // FAST_SWEEP; //  BRUTE_FORCE;
	dataType firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
	computeDistanceMap(distanceMap3D, cropped, LengthNew, WidthNew, HeightNew, distanceParamerters, distMethod);
	dataType secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
	cout << "Execution time : " << secondCpuTime - firstCpuTime << endl;
	std::string distance3D = outputPath + "distance.raw";
	store3dRawData<dataType>(distanceMap3D, LengthNew, WidthNew, HeightNew, distance3D.c_str());

	//Freeing imageData pointer 
	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
	}
	delete[] imageData;
	delete[] cropped;
	delete[] distanceMap3D;
	delete[] potential;
	delete[] path3D;
	free(seedPoint);

	
	return EXIT_SUCCESS;
}
