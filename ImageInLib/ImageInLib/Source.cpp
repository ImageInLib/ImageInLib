#pragma warning(disable : 6011)
#pragma warning(disable : 4244)
#pragma warning(disable : 4477)
#pragma warning(disable : 4267)
#pragma warning(disable : 4700)

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <time.h>

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
#include "../src/segmentation3D_subsurf.h"
#include "../src/segmentation3d_gsubsurf.h"
#include "../src/image_difference.h"
#include "../src/noise_generator.h"
#include "../src/imageInterpolation.h"
#include "distanceForPathFinding.h"
#include "distanceMaps.h"
#include "segmentation.h"
#include "segmentation2d.h"

#define thresmin 995
#define thresmax 1213


int main() {

	size_t i, j, k;

	//image Dimensions
	size_t Width = 512;
	size_t Length = 512;
	const size_t dim2D = Width * Length;

	std::string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/";
	std::string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";

	dataType* imageData = new dataType[dim2D];
	if (imageData == NULL)
		return false;
	
	//std::string inputImagePath = inputPath + "Slices304.raw";
	//std::string inputImagePath = "C:/Users/Konan Allaly/Documents/ProjectFor2dTest/output/filtered_gmc.raw";
	std::string inputImagePath = "C:/Users/Konan Allaly/Documents/MPI/input/Image.raw";
	Operation operation = LOAD_DATA;

	manageRAWFile2D<dataType>(imageData, Length, Width, inputImagePath.c_str(), operation, false);
	//rescaleToZeroOne2d(imageData, Length, Width);

	//=========== generate circle ========

	//for (i = 0; i < Length; i++) {
	//	for (j = 0; j < Width; j++) {
	//		if (sqrt((i - 32) * (i - 32) + (j - 32) * (j - 32)) < 25) {
	//			imageData[x_new(i, j, Length)] = 1;
	//		}
	//		else {
	//			imageData[x_new(i, j, Length)] = 0;
	//		}
	//	}
	//	for (j = 0; j < Width; j++) {
	//		if (sqrt((i - 32) * (i - 32) + (j - 32) * (j - 32)) < 23) {
	//			imageData[x_new(i, j, Length)] = 0;
	//		}
	//	}
	//}

	//=======================

	Image_Data2D imageToBeSegmented;
	imageToBeSegmented.imageDataPtr = imageData; imageToBeSegmented.height = Length; imageToBeSegmented.width = Width;

	dataType* initialSegment = new dataType[dim2D];
	if (initialSegment == NULL) return false;
	initialize2dArrayD(initialSegment, Length, Width, 0.0);

	point2D* center = new point2D[1];
	center->x = 187; center->y = 254;

	generateInitialSegmentationFunction(initialSegment, Length, Width, center, 0.5, 10.0);

	std::string segmentationPath = outputPath + "/segmentation/";

	std::string initialSegmPath = segmentationPath + "_seg_func_0000.raw";
	operation = STORE_DATA;
	//manageRAWFile2D<dataType>(initialSegment, Length, Width, initialSegmPath.c_str(), operation, false);

	Filter_Parameters filtering_parameters;
	filtering_parameters.timeStepSize = 1.0; filtering_parameters.h = 1.0;
	filtering_parameters.omega_c = 1.4; filtering_parameters.tolerance = 1e-3;
	filtering_parameters.maxNumberOfSolverIteration = 100; filtering_parameters.timeStepsNum = 1;
	filtering_parameters.eps2 = 0.000001; filtering_parameters.edge_detector_coefficient = 100;

	heatImplicit2dScheme(imageToBeSegmented, filtering_parameters);
	std::string filteredImagePath = outputPath + "filtered.raw";
	manageRAWFile2D<dataType>(imageToBeSegmented.imageDataPtr, Length, Width, filteredImagePath.c_str(), operation, false);

	Segmentation_Parameters segmentation_parameters;
	segmentation_parameters.h = 1.0; segmentation_parameters.coef = 100; segmentation_parameters.eps2 = 0.00001;
	segmentation_parameters.maxNoGSIteration = 100; segmentation_parameters.omega_c = 1.5; segmentation_parameters.segTolerance = 1e-3;
	segmentation_parameters.tau = 10.0; segmentation_parameters.numberOfTimeStep = 100; segmentation_parameters.maxNoOfTimeSteps = 100;
	segmentation_parameters.mod = 1; segmentation_parameters.coef_conv = 5.0, segmentation_parameters.coef_dif = 1.0;

	gsubsurf(imageToBeSegmented, initialSegment, segmentationPath, filtering_parameters, segmentation_parameters);
	//subsurf(imageToBeSegmented, initialSegment, segmentationPath, filtering_parameters, segmentation_parameters);

	//=======================

	//std::string outputImagePath = outputPath + "loaded.raw";
	//operation = STORE_DATA;
	//manageRAWFile2D<dataType>(imageData, Length, Width, outputImagePath.c_str(), operation, false);
	 
	delete[] imageData;
	delete[] initialSegment;
	delete[] center;

	return EXIT_SUCCESS;
}
