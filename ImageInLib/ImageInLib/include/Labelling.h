#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"
#include "../src/filter_params.h"
#include "../src/heat_equation.h"

	typedef struct {
		size_t x, y;
		bool label;
	}point2dLabelling;

	//for 2D image using 2D arrays
	bool labelling2D(int** imageDataPtr, int** segmentedImage, bool** statusArray, const size_t xDim, const size_t yDim, int object);

	//For 3D Images
	bool labelling3D(dataType** imageDataPtr, int** segmentedImage, bool** statusArray, const size_t xDim, const size_t yDim, const size_t zDim, dataType object);

	bool regionGrowing(dataType** imageDataPtr, dataType** segmentedImage, bool** statusArray, const size_t xDim, const size_t yDim, const size_t zDim, dataType thres_min, dataType thres_max, Point3D* seedPoint);

#ifdef __cplusplus
}
#endif
