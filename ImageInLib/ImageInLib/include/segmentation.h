#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"
#include "../src/segmentation3D_subsurf.h"

	typedef enum
	{
		SUBSURF_MODEL = 1,
		GSUBSURF_ATLAS_MODEL
	} SegmentationMethod;

	void segmentImage(Image_Data inputImageData, Segmentation_Parameters segParameters, Filter_Parameters explicit_lhe_Parameters,
		Point3D * centers, size_t no_of_centers, unsigned char * outputPathPtr, const SegmentationMethod model);

#ifdef __cplusplus
}
#endif