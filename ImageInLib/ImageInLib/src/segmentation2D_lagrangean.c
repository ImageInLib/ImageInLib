#include "segmentation2D_lagrangean.h"

bool lagrangeanExplicitOpen2DCurveSegmentation(Image_Data inputImage2D, const Lagrangean2DSegmentationParameters* pSegmentationParams, 
    unsigned char* pOutputPathPtr, Curve2D* pResultSegmentation)
{

    if (pResultSegmentation)
    {
        for (size_t i = 0; i < pResultSegmentation->numPoints; i++)
        {
            pResultSegmentation->pPoints[i] = pSegmentationParams->pInitialCondition->pPoints[i];
        }
    }

    return true;
}
