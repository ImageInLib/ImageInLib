#include "segmentation.h"
#include "segmentation3D_subsurf.h"
#include "segmentation3d_gsubsurf.h"
#include "segmentation3D_atlas.h"
#include "segmentation3d_gsubsurf.h"
#include "../include/segmentation2d.h"
#include "segmentation2D_lagrangean.h"
#include "segmentation3D_atlas.h"

void segmentImage(void * pInputImageData, void * pSegParameters, void * aSegParameters, void * aPcaParams, void * aTmpDataHolders, void * pfilterParameters,
	const SegmentationMethod model, unsigned char* outputPathPtr, void* resultSegment)
{
	switch (model)
	{
        case SUBSURF_MODEL:
        {
            Image_Data inputImageData = *(Image_Data *)pInputImageData;
            Segmentation_Parameters* pSegmentationParams = (Segmentation_Parameters*)(pSegParameters);
            FilterParameters* explicitLheParameters = (FilterParameters*)(pfilterParameters);
            subsurfSegmentation(inputImageData, pSegmentationParams->pInitialCondition, *pSegmentationParams, *explicitLheParameters, 
                pSegmentationParams->pCenters, pSegmentationParams->no_of_centers, outputPathPtr);
            break;
        }
        case GSUBSURF_MODEL:
        {
            Image_Data inputImageData = *(Image_Data*)pInputImageData;
            Segmentation_Parameters* pSegmentationParams = (Segmentation_Parameters*)pSegParameters;
            FilterParameters* explicitLheParameters = (FilterParameters*)(pfilterParameters);
            generalizedSubsurfSegmentation(inputImageData, pSegmentationParams->pInitialCondition, *pSegmentationParams, *explicitLheParameters, 
                pSegmentationParams->pCenters, pSegmentationParams->no_of_centers, outputPathPtr);
        }
        case CURVE_2D_OPEN_EXPLCIT:
            Image_Data2D inputImageData = *(Image_Data2D*)pInputImageData;
            Lagrangean2DSegmentationParameters* pSegmentationParams = (Lagrangean2DSegmentationParameters*)pSegParameters;
            Curve2D* resultSegmentationCurve = (Curve2D*)resultSegment;
            lagrangeanExplicitOpen2DCurveSegmentation(inputImageData, pSegmentationParams, outputPathPtr, resultSegmentationCurve);
        case GSUBSURF_ATLAS_MODEL:
            Image_Data imageData = *(Image_Data*)pInputImageData;
            Atlas_Segmentation_Parameters* segParameters = (Atlas_Segmentation_Parameters*)(pSegParameters);
            PCAData* pcaParams = (PCAData*)(aPcaParams);
            tmpDataHolders* tmpDataStepHolder = (tmpDataHolders*)(aTmpDataHolders);
            atlasSegmentationModel(imageData, *segParameters, pcaParams, tmpDataStepHolder, outputPathPtr);
        default:
            break;
	}
}

void segment2dImage(Image_Data2D inputImageData, dataType* initialSegment, Segmentation_Parameters segParameters, FilterParameters filteringParameters,
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