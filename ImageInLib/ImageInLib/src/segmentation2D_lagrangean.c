#include "segmentation2D_lagrangean.h"
#include "common_functions.h"
#include <stdlib.h>
#include <math.h>
#include "file.h"

bool lagrangeanExplicit2DCurveSegmentation(Image_Data2D inputImage2D, const Lagrangean2DSegmentationParameters* pSegmentationParams,
    unsigned char* pOutputPathPtr, Curve2D* pResultSegmentation)
{
    if (pSegmentationParams == NULL || pResultSegmentation == NULL) {
        return false;
    }

    const size_t sizeHeight = sizeof(dataType *) * inputImage2D.height;

    const size_t dataDimension = inputImage2D.width * inputImage2D.height;
    const size_t dataSize = dataDimension * sizeof(dataType);

    dataType* pvelocity_x = (dataType*)malloc(dataSize); // velocity component x
    dataType* pvelocity_y = (dataType*)malloc(dataSize); // velocity component y
    dataType* abs_val_grad = (dataType*)malloc(dataSize); // absolute value of gradient

    dataType* edge_detector = (dataType*)malloc(dataSize); // edge detctor
    const dataType edge_detector_coef = 10;// 00000;
    const dataType hx = 1, hy = 1;      //spatial discretization step
    const dataType hx_c = 1, hy_c = 1;  //h for central differences
    Point2D current_grad;
    FiniteVolumeSize2D finite_volume_sz = { 1.0, 1.0 };

    const size_t sizeWidth = sizeof(dataType) * inputImage2D.width;

    //get absolute value of gradient
    size_t xd = 0;
    for (size_t i = 0; i < inputImage2D.height; i++)
    {
        for (size_t j = 0; j < inputImage2D.width; j++)
        {
            getGradient2D(inputImage2D, j, i, finite_volume_sz, &current_grad);
            xd = x_new(j, i, inputImage2D.width);
            abs_val_grad[xd] = norm(current_grad);
        }
    }

    //get edge detector
    
    for (size_t i = 0; i < dataDimension; i++)
    {
        edge_detector[i] = edgeDetector(abs_val_grad[i], edge_detector_coef);
    }

 /*   dataType dx, dy, dist;
    for (size_t i = 0; i < inputImage2D.height; i++)
    {
        for (size_t j = 0; j < inputImage2D.width; j++)
        {
            xd = x_new(j, i, inputImage2D.width);
            dx = pow(j - inputImage2D.width / 2.0, 2);
            dy = pow(i - inputImage2D.height / 2.0, 2);
            dist = sqrt(dx + dy);
            edge_detector[xd] = dist;
        }
    }*/

    Image_Data2D edgeDetector = { inputImage2D.height, inputImage2D.width, edge_detector };

    //get velocity
    for (size_t i = 0; i < inputImage2D.height; i++)
    {
        for (size_t j = 0; j < inputImage2D.width; j++)
        {
            getGradient2D(edgeDetector, j, i, finite_volume_sz, &current_grad);
            xd = x_new(j, i, inputImage2D.width);
            pvelocity_x[xd] = current_grad.x;
            pvelocity_y[xd] = current_grad.y;
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

    size_t iterPt = 0;
    size_t maxIterPt = pResultSegmentation->numPoints;

    if (pSegmentationParams->open_curve) {
        iterPt = 1;
        maxIterPt = pResultSegmentation->numPoints - 1;
    }
 
    for(size_t t = 0; t < pSegmentationParams->num_time_steps; t++)
    {
        for (size_t i = iterPt; i < maxIterPt; i++)
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
            } else if (current_j >= inputImage2D.width) {
                current_j = inputImage2D.width - 1;
            }

            //it is just simple motion in vector field
            xd = x_new(current_j, current_i, inputImage2D.width);
            pResultSegmentation->pPoints[i].x = pResultSegmentation->pPoints[i].x -
                pSegmentationParams->time_step_size * pvelocity_x[xd];

            pResultSegmentation->pPoints[i].y = pResultSegmentation->pPoints[i].y -
                pSegmentationParams->time_step_size * pvelocity_y[xd];
        }
    }

    free(pvelocity_x);
    free(pvelocity_y);

    free(abs_val_grad);
    free(edge_detector);

    return true;
}
