#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"
#include "../src/filter_params.h"
#include "../src/heat_equation.h"


	//for 2D image using 2D arrays
	bool labelling2D(int** imageDataPtr, int** segmentedImage, bool** statusArray, const size_t xDim, const size_t yDim, int object);

	//For 3D Images
	bool labelling3D(dataType** imageDataPtr, int** segmentedImage, bool** statusArray, const size_t xDim, const size_t yDim, const size_t zDim, dataType object);

	bool regionGrowing(dataType** imageDataPtr, int** segmentedImage, bool** statusArray, const size_t xDim, const size_t yDim, const size_t zDim, dataType thres_min, dataType thres_max, Filter_Parameters smoothParameters);
	

	//===========================
	//Statistics after labelling
	bool labellingStats(int** segmentedImage, int* CountingArray, const size_t xDim, const size_t yDim, const char* pathPtr);

	bool minorRegionRemoval(int** imageDataPtr, int** segmentedImage, int* CountingArray, const size_t xDim, const size_t yDim, int size);

	bool initialization2dArray(int** imageDataPtr, const size_t xDim, const size_t yDim, int value);

	bool rescalingTo2D(int** imageDataPtr, int xDim, int yDim, int minNew, int maxNew);

	bool sortArray(int* valuePtr, int sizeArray);
	//=============================

	////Compute image gradient
	//bool imageGradient(dataType imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType thresMin, dataType threshMax);

#ifdef __cplusplus
}
#endif
