#include "segmentation.h"
#include "segmentation3D_subsurf.h"
#include "segmentation3D_atlas.h"

void segmentImage(Image_Data inputImageData, Segmentation_Parameters segParameters, Filter_Parameters explicit_lhe_Parameters,
	Point3D * centers, size_t no_of_centers, unsigned char * outputPathPtr, const SegmentationMethod model)
{
	switch (model)
	{
	case SUBSURF_MODEL:
		subsurfSegmentation(inputImageData, segParameters, explicit_lhe_Parameters, centers, no_of_centers, outputPathPtr);
		break;
	default:
		break;
	}
}