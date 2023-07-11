#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"
#include "../src/data_storage.h"
#include "../src/heat_equation.h"
#include "../src/filter_params.h"
#include "../src/segmentation3D_subsurf.h"

	typedef struct
	{
		dataType* East;
		dataType* West;
		dataType* North;
		dataType* South;
	} neighPtrs;

	typedef struct
	{
		size_t x, y;
	} point2d;

	//dataType min(dataType a, dataType b);

	dataType l2norm(dataType* arrayPtr1, dataType* arrayPtr2, const size_t height, const size_t width, dataType h);

	bool rescaleToZeroOne2d(dataType* imageDataPtr, const size_t height, const size_t width);

	bool generateInitialSegmentationFunction(dataType* imageDataPtr, const size_t height, const size_t width, point2d* center, dataType v, dataType R);

	bool set2dDirichletBoundaryCondition(dataType* imageDataPtr, const size_t height, const size_t width);

	bool computeNormOfGradientDiamondCells(dataType* arrayPtr, neighPtrs neigbours, const size_t height, const size_t width, dataType h);

	bool epsilonRegularization(neighPtrs neighbours, const size_t height, const size_t width, dataType epsilon);

	bool subsurf(Image_Data2D imageData, dataType* initialSegment, const char * segmentPath, const Filter_Parameters smooth_parms, Segmentation_Parameters seg_parms);

	bool gsubsurf(Image_Data2D imageData, dataType* initialSegment, const char * segmentPath, const Filter_Parameters smooth_parms, Segmentation_Parameters seg_parms);

#ifdef __cplusplus
}
#endif