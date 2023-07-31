#include "segmentation2D_lagrangean.h"
//#include "common_functions.h"
#include <stdlib.h>


bool lagrangeanExplicitOpen2DCurveSegmentation(Image_Data inputImage2D, const Lagrangean2DSegmentationParameters* pSegmentationParams,
    unsigned char* pOutputPathPtr, Curve2D* pResultSegmentation)
{
    if (pSegmentationParams == NULL || pResultSegmentation == NULL) {
        return false;
    }

    const size_t sizeHeight = sizeof(dataType *) * inputImage2D.height;

    dataType** pvelocity_x = (dataType**)malloc(sizeHeight); // velocity component x
    dataType** pvelocity_y = (dataType**)malloc(sizeHeight); // velocity component x

    const size_t sizeWidth = sizeof(dataType) * inputImage2D.width;

    for (size_t i = 0; i < inputImage2D.height; i++)
    {
        pvelocity_x[i] = (dataType*)malloc(sizeWidth);
        pvelocity_y[i] = (dataType*)malloc(sizeWidth);
    }

    //definition of velocity field
    for (size_t i = 0; i < inputImage2D.height; i++)
    {
        for (size_t j = 0; j < inputImage2D.width; j++)
        {
            pvelocity_x[i][j] = (dataType)0;
            pvelocity_y[i][j] = (dataType)1.0;
        }
    }

    size_t curveSize = sizeof(Curve2D);
    size_t curvePoint2DSize = sizeof(CurvePoint2D);
    size_t curvePoint2DPartsSize = sizeof(Point2D) + sizeof(bool) + 2 * sizeof(void*);

    size_t current_i = 0;
    size_t current_j = 0;
    for (size_t i = 0; i < (pResultSegmentation->numPoints); i++)
    {
        pResultSegmentation->pPoints[i].x = pSegmentationParams->pinitial_condition->pPoints[i].x;
        pResultSegmentation->pPoints[i].y = pSegmentationParams->pinitial_condition->pPoints[i].y;
    }


    for(size_t t = 0; t < pSegmentationParams->num_time_steps; t++)
    {
        for (size_t i = 0; i < pResultSegmentation->numPoints; i++)
        {
            //let us move by curve points just in vector field
            current_i = (size_t)(pResultSegmentation->pPoints[i].y + 0.5);
            current_j = (size_t)(pResultSegmentation->pPoints[i].x + 0.5);

            //let us keep points inside the image
            if (current_i < 0) {
                current_i = 0;
            } else if (current_i >= inputImage2D.height) {
                current_i = inputImage2D.height - 1;
            }

            if (current_j < 0) {
                current_j = 0;
            } else if (current_j < inputImage2D.width) {
                current_j = inputImage2D.width - 1;
            }

            //it is just simple motion in vector field
            pResultSegmentation->pPoints[i].x = pResultSegmentation->pPoints[i].x +
                pSegmentationParams->time_step_size * pvelocity_x[current_i][current_j];

            pResultSegmentation->pPoints[i].y = pResultSegmentation->pPoints[i].y +
                pSegmentationParams->time_step_size * pvelocity_y[current_i][current_j];
        }
    }

    for (size_t i = 0; i < inputImage2D.height; i++)
    {
        free(pvelocity_x[i]);
        free(pvelocity_y[i]);
    }

    free(pvelocity_x);
    free(pvelocity_y);

    return true;
}
