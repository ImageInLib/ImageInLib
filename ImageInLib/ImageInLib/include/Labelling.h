#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"


	//for 2D image using 2D arrays
	bool regionLabelling2D(int** imageDataPtr, int** segmentedImage, const size_t xDim, const size_t yDim, int background, int object);

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

	//for 3D images
	bool fixEquivalence(int** segmentedImage, const size_t xDim, size_t x, size_t y, size_t z, int minV, int maxV, bool parallize, size_t nbtreads);

	bool regionLabelling3D(dataType** imageDataPtr, int** segmentedImage, const size_t xDim, const size_t yDim, const size_t zDim, dataType background, dataType object, bool parallize, size_t nbtreads);

	bool labelling3D(dataType** imageDataPtr, int** segmentedImage, bool** statusArray, const size_t xDim, const size_t yDim, const size_t zDim, dataType object);



#ifdef __cplusplus
}
#endif
