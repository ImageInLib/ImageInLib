#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"
#include "../src/segmentation3D_subsurf.h"
#include "segmentation2d.h"
#include "segmentation2D_lagrangean.h"

	typedef enum
	{
		SUBSURF_MODEL = 1,
		GSUBSURF_MODEL,
		LABELING,
		GSUBSURF_ATLAS_MODEL,
		CURVE_2D_OPEN_EXPLCIT
	} SegmentationMethod;

	void segmentImage(Image_Data inputImageData, void* pSegParameters, void* pfilterParameters,
		const SegmentationMethod model, unsigned char * outputPathPtr, void * resultSegment);

	void segment2dImage(Image_Data2D inputImageData, dataType* initialSegment, Segmentation_Parameters segParameters, Filter_Parameters filteringParameters,
		point2d* centers, const char* outputPathPtr, const SegmentationMethod model);

#ifdef __cplusplus
}
#endif