#include "segmentation.h"
#include "segmentation3D_subsurf.h"
#include "segmentation3d_gsubsurf.h"
#include "segmentation3D_atlas.h"

void segmentImage(Image_Data inputImageData, dataType** initialSegment, Segmentation_Parameters segParameters, Filter_Parameters explicit_lhe_Parameters,
	Point3D * centers, size_t no_of_centers, unsigned char * outputPathPtr, const SegmentationMethod model)
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