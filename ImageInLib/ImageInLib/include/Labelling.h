#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"


	//for 2D image using 2D arrays
	bool regionLabelling(int** imageDataPtr, int** segmentedImage, int xDim, int yDim, int background, int object);

	bool labelling(int* imageDataPtr, int* segmentedImage, bool* statusArray, int xDim, int yDim, int object);

	bool minorRegionRemoval(int** imageDataPtr, int** segmentedImage, int* CountingArray, int xDim, int yDim, int size);

	bool initialization2dArray(int** imageDataPtr, int xDim, int yDim, int value);

	bool rescalingTo2D(int** imageDataPtr, int xDim, int yDim, int minNew, int maxNew);

	bool sortArray(int* valuePtr, int sizeArray);


	//for 3D images
	bool fixEquivalence(dataType** segmentedImage, const size_t xDim, size_t x, size_t y, size_t z, dataType minV, dataType maxV, bool parallize, size_t nbtreads);

	bool regionLabelling3D(dataType** imageDataPtr, dataType** segmentedImage, const size_t xDim, const size_t yDim, const size_t zDim, dataType background, dataType object, bool parallize, size_t nbtreads);

	bool labelling3D(dataType* imageDataPtr, dataType* segmentedImage, bool* statusArray, int xDim, int yDim, int zDim, dataType object);



#ifdef __cplusplus
}
#endif
