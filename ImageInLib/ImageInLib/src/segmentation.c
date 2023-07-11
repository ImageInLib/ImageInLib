#include "segmentation.h"
#include "segmentation3D_subsurf.h"
#include "segmentation3d_gsubsurf.h"
#include "segmentation3D_atlas.h"
#include "segmentation3d_gsubsurf.h"
#include "../include/segmentation2d.h"
#include "segmentation2D_lagrangean.h"

void segmentImage(Image_Data inputImageData, void * pSegParameters, void * pfilterParameters,
	const SegmentationMethod model, unsigned char* outputPathPtr, void* resultSegment)
{
	switch (model)
	{
        case SUBSURF_MODEL:
        {
            Segmentation_Parameters* pSegmentationParams = (Segmentation_Parameters*)(pSegParameters);
            Filter_Parameters* explicitLheParameters = (Filter_Parameters*)(pfilterParameters);
            subsurfSegmentation(inputImageData, pSegmentationParams->pInitialCondition, *pSegmentationParams, *explicitLheParameters, 
                pSegmentationParams->pCenters, pSegmentationParams->no_of_centers, outputPathPtr);
            break;
        }
        case GSUBSURF_MODEL:
        {
            Segmentation_Parameters* pSegmentationParams = (Segmentation_Parameters*)pSegParameters;
            Filter_Parameters* explicitLheParameters = (Filter_Parameters*)(pfilterParameters);
            generalizedSubsurfSegmentation(inputImageData, pSegmentationParams->pInitialCondition, *pSegmentationParams, *explicitLheParameters, 
                pSegmentationParams->pCenters, pSegmentationParams->no_of_centers, outputPathPtr);
        }
        case CURVE_2D_OPEN_EXPLCIT:
            Lagrangean2DSegmentationParameters* pSegmentationParams = (Lagrangean2DSegmentationParameters*)pSegParameters;
            Curve2D* resultSegmentationCurve = (Curve2D*)resultSegment;
            lagrangeanExplicitOpen2DCurveSegmentation(inputImageData, pSegmentationParams, outputPathPtr, resultSegmentationCurve);
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