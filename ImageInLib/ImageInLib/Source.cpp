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

//#define thresmin 995
//#define thresmax 1213


int main() {

	size_t i, j, k, xd, m, n;

	//image Dimensions
	const size_t Width = 64; /*512*/
	const size_t Length = 64; /*512*/
	const size_t Height = 64; /*406*/ /*607*/ /*508*/
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
	//std::string inputImagePath = inputPath + "filteredMC.raw";

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

	//load3dArrayRAW<dataType>(imageData, Length, Width, Height, inputImagePath.c_str());
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

	//--------------- Generate artificial Image -----

	//// Artificial Sphere
	size_t x = Length / 2, y = Width / 2, z = Height / 2;
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

	Point3D c = { x,y,z };
	dataType radius = 20;
	dataType small_radius = 6;
	//size_t pitch_ball = 30;
	dataType fill_value = 1.0;

	std::string artificialImage = outputPath + "ballsWithHoles.raw";

	//generateSphere(imageData, c, Length, Width, Height, radius, fill_value, artificialImage.c_str());
	generateSphereWithSixHoles(imageData, c, Length, Width, Height, radius, small_radius, fill_value, artificialImage.c_str());
	//ballsOnHelix(imageData, pitch_ball, Length, Width, Height, radius, small_radius, fill_value);

	store3dRawData<dataType>(imageData, Length, Width, Height, artificialImage.c_str());

	//--------------- Test segmentation -------------

	Point3D* centerSeg = new Point3D[1];
	size_t nb_centers = 1;
	centerSeg->x = x; centerSeg->y = y; centerSeg->z = z;

	dataType** initialSegment = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		initialSegment[k] = new dataType[Width * Height];
	}

	////If we want to start with the segmentatation function originally implemented in the library
	generateInitialSegmentationFunctionForMultipleCentres(initialSegment, Length, Width, Height, centerSeg, 0.5, 10, nb_centers);

	std::string segmFolderPath = outputPath + "segmentation/";
	//store3dRawData<dataType>(initialSegment, Length, Width, Height, (segmFolderPath + std::string("_seg_func_000.raw")).c_str());

	Image_Data segment; segment.height = Height; segment.length = Length; segment.width = Width; segment.imageDataPtr = imageData;
	//rescaleNewRange(segment.imageDataPtr, Length, Width, Height, 0, 1);
	Segmentation_Parameters segmentParameters; segmentParameters.coef = 10000; segmentParameters.eps2 = 1e-6; segmentParameters.gauss_seidelTolerance = 1e-3;
	segmentParameters.h = 1.0; segmentParameters.maxNoGSIteration = 100; segmentParameters.maxNoOfTimeSteps = 50; segmentParameters.mod = 1;
	segmentParameters.numberOfTimeStep = 50; segmentParameters.omega_c = 1.5; segmentParameters.segTolerance = 1e-4; segmentParameters.tau = 8;

	Filter_Parameters filterParameters; filterParameters.coef = 1e-2; filterParameters.edge_detector_coefficient = 1; filterParameters.eps2 = 1e-4;
	filterParameters.h = 1.0; filterParameters.maxNumberOfSolverIteration = 100; filterParameters.omega_c = 1.5; filterParameters.p = 1;
	filterParameters.sigma = 1e-3; filterParameters.timeStepSize = 1.2; filterParameters.timeStepsNum = 1; filterParameters.tolerance = 4 * 1e-4;

	unsigned char outputPathPtr[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/";
	//subsurfSegmentation(segment, initialSegment, segmentParameters, filterParameters, centerSeg, nb_centers, outputPathPtr);
	generalizedSubsurfSegmentation(segment, initialSegment, segmentParameters, filterParameters, centerSeg, nb_centers, outputPathPtr, 0.5, 0.2);

	delete[] centerSeg;
	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
		delete[] initialSegment[k];
	}
	delete[] imageData;
	delete[] initialSegment;

	return EXIT_SUCCESS;
}
