#include "segmentation2D_lagrangean.h"

bool lagrangeanExplicitOpen2DCurveSegmentation(Image_Data inputImage2D, const Lagrangean2DSegmentationParameters* pSegmentationParams,
    unsigned char* pOutputPathPtr, Curve2D* pResultSegmentation)
{
    if (pSegmentationParams == NULL || pResultSegmentation == NULL) {
        return false;
    }

    for(size_t t = 0; t < pSegmentationParams->num_time_steps; t++)
    {
        for (size_t i = 0; i < pResultSegmentation->numPoints; i++)
        {
            pResultSegmentation->pPoints[i] = pSegmentationParams->pinitial_condition->pPoints[i];
        }
    }

    return true;
}
