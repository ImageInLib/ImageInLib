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
		GSUBSURF_ATLAS_MODEL,
		CURVE_2D_OPEN_EXPLCIT
	} SegmentationMethod;

	// Structure that holds the parameters used during 2d Lagrangean segmentation process.
	typedef struct
	{
		size_t num_time_steps;	// Number of time steps
		Curve2D* pinitial_condition;	// Initial curve
		size_t num_points;// Number of initial segmentation curve points
	} Lagrangean2DSegmentationParameters;


	void segmentImage(Image_Data inputImageData, void* pSegParameters, void* pfilterParameters,
		const SegmentationMethod model, unsigned char * outputPathPtr, void * resultSegment);

	void segment2dImage(Image_Data2D inputImageData, dataType* initialSegment, Segmentation_Parameters segParameters, FilterParameters filteringParameters,
		point2d* centers, const char* outputPathPtr, const SegmentationMethod model);

#ifdef __cplusplus
}
#endif