#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"


	//for 2D image using 2D arrays
	bool labelling2D(int** imageDataPtr, int** segmentedImage, bool** statusArray, const size_t xDim, const size_t yDim, int object);

	//Statistics after labelling
	bool labellingStats(int** segmentedImage, int* CountingArray, const size_t xDim, const size_t yDim, const char* pathPtr);

	//Erosion ---> Shrink objects and remove boundaries pixels
	bool erosion2D(int** imageDataPtr, const size_t xDim, const size_t yDim, int background, int object);

	//===========================
	bool minorRegionRemoval(int** imageDataPtr, int** segmentedImage, int* CountingArray, const size_t xDim, const size_t yDim, int size);

	bool initialization2dArray(int** imageDataPtr, const size_t xDim, const size_t yDim, int value);

	bool rescalingTo2D(int** imageDataPtr, int xDim, int yDim, int minNew, int maxNew);

	bool sortArray(int* valuePtr, int sizeArray);
	//=============================

	//For 3D Images
	bool labelling3D(dataType** imageDataPtr, int** segmentedImage, bool** statusArray, const size_t xDim, const size_t yDim, const size_t zDim, dataType object);

	bool erosion3D(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background);

	bool dilatation3D(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background);


#ifdef __cplusplus
}
#endif
