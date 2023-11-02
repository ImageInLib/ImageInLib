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

    CurvePoint2D* pOldCurve = (CurvePoint2D*)malloc(sizeof(CurvePoint2D) * pResultSegmentation->numPoints);
    Curve2D oldSegmentation = { pOldCurve, pResultSegmentation->numPoints };

    const size_t sizeHeight = sizeof(dataType *) * inputImage2D.height;

    const size_t dataDimension = inputImage2D.width * inputImage2D.height;
    const size_t dataSize = dataDimension * sizeof(dataType);

    dataType* pgrad_x = (dataType*)malloc(dataSize); // velocity component x
    dataType* pgrad_y = (dataType*)malloc(dataSize); // velocity component y
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

    dataType dx, dy, dist, max_dist = 0;

    Image_Data2D edgeDetector = { inputImage2D.height, inputImage2D.width, edge_detector };

    //get velocity
    for (size_t i = 0; i < inputImage2D.height; i++)
    {
        for (size_t j = 0; j < inputImage2D.width; j++)
        {
            getGradient2D(edgeDetector, j, i, finite_volume_sz, &current_grad);
            xd = x_new(j, i, inputImage2D.width);
            pgrad_x[xd] = current_grad.x;
            pgrad_y[xd] = current_grad.y;
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
    size_t maxIterPt = oldSegmentation.numPoints;

    if (pSegmentationParams->open_curve) {
        iterPt = 1;
        maxIterPt = oldSegmentation.numPoints - 1;
    }
 
    dataType vx, vy;
    dataType rx_l, rx_g, ry_l, ry_g;
    dataType nx, ny, dot, norm, tx, ty;

    for(size_t t = 0; t < pSegmentationParams->num_time_steps; t++)
    {
        for (size_t i = 0; i < (oldSegmentation.numPoints); i++)
        {
            oldSegmentation.pPoints[i].x = pResultSegmentation->pPoints[i].x;
            oldSegmentation.pPoints[i].y = pResultSegmentation->pPoints[i].y;
        }


        for (size_t i = iterPt; i < maxIterPt; i++)
        {
            //let us move by curve points just in vector field
            current_i = (size_t)(oldSegmentation.pPoints[i].y + 0.5);
            current_j = (size_t)(oldSegmentation.pPoints[i].x + 0.5);

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

            xd = x_new(current_j, current_i, inputImage2D.width);

            if (i == 0) {
                rx_l = oldSegmentation.pPoints[pResultSegmentation->numPoints - 1].x;
                ry_l = oldSegmentation.pPoints[pResultSegmentation->numPoints - 1].y;
                rx_g = oldSegmentation.pPoints[1].x;
                ry_g = oldSegmentation.pPoints[1].y;
            }
            else if (i == oldSegmentation.numPoints - 1) {
                rx_l = oldSegmentation.pPoints[i - 1].x;
                ry_l = oldSegmentation.pPoints[i - 1].y;
                rx_g = oldSegmentation.pPoints[0].x;
                ry_g = oldSegmentation.pPoints[0].y;
            }
            else {
                rx_l = oldSegmentation.pPoints[i - 1].x;
                ry_l = oldSegmentation.pPoints[i - 1].y;
                rx_g = oldSegmentation.pPoints[i + 1].x;
                ry_g = oldSegmentation.pPoints[i + 1].y;
            }

            norm = sqrt(pow(rx_l - rx_g, 2) + pow(ry_l - ry_g, 2));
            tx = (rx_g - rx_l) / norm;
            ty = (ry_g - ry_l) / norm;

            nx = ty;
            ny = -tx;

            dot = pgrad_x[xd] * nx + pgrad_y[xd] * ny;
            
            vx = -pgrad_x[xd];
            vy = -pgrad_y[xd];

            dot = vx* nx + vy * ny;

            vx = dot * nx;
            vy = dot * ny;

            //it is just simple motion in vector field

            pResultSegmentation->pPoints[i].x = oldSegmentation.pPoints[i].x +
                pSegmentationParams->time_step_size * vx;

            pResultSegmentation->pPoints[i].y = oldSegmentation.pPoints[i].y +
                pSegmentationParams->time_step_size * vy;
        }
    }

    free(pOldCurve);

    free(pgrad_x);
    free(pgrad_y);

    free(abs_val_grad);
    free(edge_detector);

    return true;
}
