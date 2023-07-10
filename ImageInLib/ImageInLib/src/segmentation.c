#include "segmentation.h"
#include "segmentation3D_subsurf.h"
#include "segmentation3d_gsubsurf.h"
#include "segmentation3D_atlas.h"
#include "segmentation3d_gsubsurf.h"
#include "../include/segmentation2d.h"

void segmentImage(Image_Data inputImageData, dataType** initialSegment, Segmentation_Parameters segParameters, Filter_Parameters explicit_lhe_Parameters,
	Point3D* centers, size_t no_of_centers, unsigned char* outputPathPtr, const SegmentationMethod model)
{
	switch (model)
	{
	case SUBSURF_MODEL:
		subsurfSegmentation(inputImageData, initialSegment, segParameters, explicit_lhe_Parameters, centers, no_of_centers, outputPathPtr);
		break;
	case GSUBSURF_MODEL:
		generalizedSubsurfSegmentation(inputImageData, initialSegment, segParameters, explicit_lhe_Parameters, centers, no_of_centers, outputPathPtr);
	default:
		break;
	}
}

void segment2dImage(Image_Data2D inputImageData, dataType* initialSegment, Segmentation_Parameters segParameters, Filter_Parameters filteringParameters,
	point2d* centers, const char* outputPathPtr, const SegmentationMethod model) 
{
	switch (model) 
	{
	case SUBSURF_MODEL:
		subsurf(inputImageData, initialSegment, outputPathPtr, filteringParameters, segParameters);
		break;
	case GSUBSURF_MODEL:
		gsubsurf(inputImageData, initialSegment, outputPathPtr, filteringParameters, segParameters);
	}
}