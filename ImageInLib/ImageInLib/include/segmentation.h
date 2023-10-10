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
		dataType time_step_size;	//discrete time step size
	} Lagrangean2DSegmentationParameters;

	/// <summary>
	/// The method segments the image according to given inoput paraneters
	/// </summary>
	/// <param name="pInputImageData">Pointer to Image_Data for 3D cases or pointer to Image_Data2D respectively representing an image to be segmented</param>
	/// <param name="pSegParameters">Pointer to Segmentation_Parameters representing segmentation parameters</param>
	/// <param name="pFilterParameters">Pointer to FilterParameters representing pre-filtering parameters</param>
	/// <param name="model">SegmentationMethod representing available segmentation methods</param>
	/// <param name="outputPathPtr">A pointer to unsigned char representing a destination path for resulting segmentation </param>
	/// <param name="pResultSegment">A pointer to Curve2D where results are returned. It is used by CURVE_2D_OPEN_EXPLCIT only</param>
	//TODO: utilize pResultSegment for another segmentation methods as well
	void segmentImage(void * pInputImageData, void* pSegParameters, void* pFilterParameters,
		const SegmentationMethod model, unsigned char * outputPathPtr, void * pResultSegment);

	void segment2dImage(Image_Data2D inputImageData, dataType* initialSegment, Segmentation_Parameters segParameters, FilterParameters filteringParameters,
		point2d* centers, const char* outputPathPtr, const SegmentationMethod model);

#ifdef __cplusplus
}
#endif