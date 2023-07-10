#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"
#include "../src/segmentation3D_subsurf.h"
#include "segmentation2d.h"

	typedef enum
	{
		SUBSURF_MODEL = 1,
		GSUBSURF_MODEL,
		LABELING,
		GSUBSURF_ATLAS_MODEL
	} SegmentationMethod;

	void segmentImage(Image_Data inputImageData, dataType** initialSegment, Segmentation_Parameters segParameters, Filter_Parameters filteringParameters,
		Point3D * centers, size_t no_of_centers, unsigned char * outputPathPtr, const SegmentationMethod model);

	void segment2dImage(Image_Data2D inputImageData, dataType* initialSegment, Segmentation_Parameters segParameters, Filter_Parameters filteringParameters,
		point2d* centers, const char* outputPathPtr, const SegmentationMethod model);

#ifdef __cplusplus
}
#endif