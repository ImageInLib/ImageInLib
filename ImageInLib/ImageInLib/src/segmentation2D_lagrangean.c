#include "segmentation2D_lagrangean.h"
#include "common_functions.h"
#include <stdlib.h>
#include <math.h>
#include "file.h"
#include "solvers.h"

bool evolveBySingleStep(Image_Data2D * pimage,Image_Data2D* pedge, LinkedCurve* plinked_curve, SchemeData* pscheme_data, const Lagrangean2DSegmentationParameters* pparams);
void calculateCurvature(LinkedCurve* plinked_curve, SchemeData* pscheme_data);
void normal_velocity(Image_Data2D * pimage, Image_Data2D* pedge, LinkedCurve* plinked_curve, SchemeData* pscheme_data,
    void(*pget_velocity)(Image_Data2D*, double, double, double*, double*),
    void(*pget_g2)(Image_Data2D*, double, double, double*),
    const double eps, const double lambda);
void tang_velocity(LinkedCurve* plinked_curve, SchemeData* pscheme_data, const double omega);
bool semiCoefficients(LinkedCurve* plinked_curve, SchemeData* pscheme_data, const double eps, const double dt);

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

    dataType* edge_detector = (dataType*)malloc(dataSize); // edge detector
    const dataType edge_detector_coef = 10;// 00000;
    dataType* similar_intensity_detector = NULL;
    if (!pSegmentationParams->open_curve)
    {
        similar_intensity_detector = (dataType*)malloc(dataSize); // similar intensity detector
    }

    const dataType similar_intensity_detector_coef = 10;// 00000;
    const dataType hx = 1, hy = 1;      //spatial discretization step
    const dataType hx_c = 1, hy_c = 1;  //h for central differences
    Point2D current_grad;
    FiniteVolumeSize2D finite_volume_sz = { 1.0, 1.0 };

    dataType lambda = 1.0;
    
    if (!pSegmentationParams->open_curve) {
        lambda = pSegmentationParams->lambda;
    }

    const size_t sizeWidth = sizeof(dataType) * inputImage2D.width;

    //get absolute value of gradient
    size_t xd = 0;
    for (size_t i = 0; i < inputImage2D.height; i++)
    {
        for (size_t j = 0; j < inputImage2D.width; j++)
        {
            getGradient2D(inputImage2D.imageDataPtr, inputImage2D.width, inputImage2D.height, j, i, finite_volume_sz, &current_grad);
            xd = x_new(j, i, inputImage2D.width);
            abs_val_grad[xd] = norm(current_grad);
        }
    }

    //get edge detector
    
    for (size_t i = 0; i < dataDimension; i++)
    {
        edge_detector[i] = edgeDetector(abs_val_grad[i], edge_detector_coef);
    }

    Point2D centroid = getCurveCentroid(pSegmentationParams->pinitial_condition);
    size_t centroid_i = (size_t)(centroid.y + 0.5);
    size_t centroid_j = (size_t)(centroid.x + 0.5);
    dataType ref_intensity = inputImage2D.imageDataPtr[x_new(centroid_j, centroid_i, inputImage2D.width)];

    if (!pSegmentationParams->open_curve)
    {
        for (size_t i = 0; i < dataDimension; i++)
        {
            similar_intensity_detector[i] = similarIntensityDetector(inputImage2D.imageDataPtr[i], ref_intensity, similar_intensity_detector_coef);
        }
    }

    Image_Data2D edgeDetector = { inputImage2D.height, inputImage2D.width, edge_detector };

    //get velocity
    for (size_t i = 0; i < inputImage2D.height; i++)
    {
        for (size_t j = 0; j < inputImage2D.width; j++)
        {
            getGradient2D(edgeDetector.imageDataPtr, edgeDetector.width, edgeDetector.height, j, i, finite_volume_sz, &current_grad);
            xd = x_new(j, i, inputImage2D.width);
            pgrad_x[xd] = current_grad.x;
            pgrad_y[xd] = current_grad.y;
        }
    }

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
 
    //_l - lower, _g - greater, _c - current
    dataType vx, vy;
    dataType rx_l, rx_g, ry_l, ry_g, rx_c, ry_c;
    dataType nx, ny, dot, norm, tx, ty;
    dataType h_g, h_c;

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

            rx_c = oldSegmentation.pPoints[i].x;
            ry_c = oldSegmentation.pPoints[i].y;

            if (i == 0) {
                rx_l = oldSegmentation.pPoints[oldSegmentation.numPoints - 1].x;
                ry_l = oldSegmentation.pPoints[oldSegmentation.numPoints - 1].y;
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

            norm = (dataType)sqrt(pow(rx_l - rx_g, 2) + pow(ry_l - ry_g, 2));
            tx = (rx_g - rx_l) / norm;
            ty = (ry_g - ry_l) / norm;

            nx = ty;
            ny = -tx;

            dot = pgrad_x[xd] * nx + pgrad_y[xd] * ny;
            
            //basic velocity - gradients of edge detector
            vx = -pgrad_x[xd];
            vy = -pgrad_y[xd];

            //projection of velocity field in direction of normal vector
            dot = vx* nx + vy * ny;

            h_c = (dataType)sqrt(pow(rx_c - rx_l,2) + pow(ry_c - ry_l,2));
            h_g = (dataType)sqrt(pow(rx_g - rx_c, 2) + pow(ry_g - ry_c, 2));

            //mu * normal vector * projection + eps * normal * curvature
            vx = pSegmentationParams->mu * (lambda * dot * nx + 
                ((dataType)1.0 - lambda) * similar_intensity_detector[xd] * nx) +
                pSegmentationParams->eps * (dataType)(2.0 / (h_g + h_c)) * ((rx_g - rx_c) / h_g - ((rx_c - rx_l) / h_c));
            vy = pSegmentationParams->mu * dot * ny + (lambda * dot * ny +
                ((dataType)1.0 - lambda) * similar_intensity_detector[xd] * ny) +
                pSegmentationParams->eps * (dataType)(2.0 / (h_g + h_c)) * ((ry_g - ry_c) / h_g - ((ry_c - ry_l) / h_c));

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

    if (similar_intensity_detector)
    {
        free(similar_intensity_detector);
    }

    return true;
}

bool lagrangeanSemiImplicit2DCurveSegmentation(Image_Data2D inputImage2D, const Lagrangean2DSegmentationParameters* pSegmentationParams,
    unsigned char* pOutputPathPtr, Curve2D* pResultSegmentation)
{
    if (pSegmentationParams == NULL || pSegmentationParams->open_curve || pResultSegmentation == NULL) {
        return false;
    }

    if (pSegmentationParams->num_points < 3) {
        return false;
    }

    const dataType similar_intensity_detector_coef = 10;// 00000;
    const dataType edge_detector_coef = 10;
    const dataType hx = 1, hy = 1;      //spatial discretization step
    const dataType hx_c = 1, hy_c = 1;  //h for central differences
    Point2D current_grad;
    FiniteVolumeSize2D finite_volume_sz = { 1.0, 1.0 };
    const size_t dataDimension = inputImage2D.width * inputImage2D.height;
    const size_t dataSize = dataDimension * sizeof(dataType);
    dataType* ptmp = (dataType*)malloc(dataSize); // absolute value of gradient

    const size_t sizeWidth = sizeof(dataType) * inputImage2D.width;

    //get absolute value of gradient
    size_t xd = 0;
    for (size_t i = 0; i < inputImage2D.height; i++)
    {
        for (size_t j = 0; j < inputImage2D.width; j++)
        {
            getGradient2D(inputImage2D.imageDataPtr, inputImage2D.width, inputImage2D.height, j, i, finite_volume_sz, &current_grad);
            xd = x_new(j, i, inputImage2D.width);
            ptmp[xd] = norm(current_grad);
        }
    }

    //get edge detector

    for (size_t i = 0; i < dataDimension; i++)
    {
        ptmp[i] = edgeDetector(ptmp[i], edge_detector_coef);
    }

    Image_Data2D edge = { inputImage2D.height, inputImage2D.width, ptmp };

    Point2D centroid = getCurveCentroid(pSegmentationParams->pinitial_condition);
    size_t centroid_i = (size_t)(centroid.y + 0.5);
    size_t centroid_j = (size_t)(centroid.x + 0.5);
    dataType ref_intensity = inputImage2D.imageDataPtr[x_new(centroid_j, centroid_i, inputImage2D.width)];

    resetIDGenerator();
    //let us consider single curve without topological changes
 
    bool isOrientedPositively = true; //does not metter for open curves

    if (!pSegmentationParams->open_curve)
    {
        isOrientedPositively = isCurveOrientedPositively(pSegmentationParams->pinitial_condition);
    }
    
    LinkedCurve linked_curve = createLinkedCurve();
    initializeLinkedCurve(pSegmentationParams->pinitial_condition, &linked_curve, !isOrientedPositively, !pSegmentationParams->open_curve);
    
    if (!pSegmentationParams->open_curve)
    {
        //este treba premysliet, ktora cast tychto dat bude v zretazenom zozname
        size_t length_of_data = linked_curve.number_of_points + 2;
        SchemeData* pscheme_data = (SchemeData*)calloc(length_of_data, sizeof(SchemeData));
        //current_point = linked_curve.first_point;

        for (size_t it = 1, res_it = 0; it <= pSegmentationParams->num_time_steps; it++)
        {
            if (length_of_data < linked_curve.number_of_points + 2)
            {
                free(pscheme_data);
                length_of_data = linked_curve.number_of_points + 2;
                pscheme_data = (SchemeData*)calloc(length_of_data, sizeof(SchemeData));
            }
            //evolve curve
            evolveBySingleStep(&inputImage2D, &edge, &linked_curve, pscheme_data, pSegmentationParams);
        }

        free(pscheme_data);
        free(ptmp);

        LinkedPoint* pt = linked_curve.first_point;

        for (size_t i = 0; i < linked_curve.number_of_points; i++)
        {
            pResultSegmentation->pPoints[i].x = (dataType)pt->x;
            pResultSegmentation->pPoints[i].y = (dataType)pt->y;
            pt = pt->next;
        }

        releaseLinkedCurve(&linked_curve);

        return true;

    }

    return false;
}

bool evolveBySingleStep(Image_Data2D * pimage, Image_Data2D* pedge, LinkedCurve* plinked_curve, SchemeData* pscheme_data, const Lagrangean2DSegmentationParameters* pparams)
{
    if (plinked_curve == NULL ||
        pscheme_data == NULL ||
        pparams == NULL ||
        pimage->imageDataPtr == NULL || 
        pedge->imageDataPtr == NULL) {
        return false;
    }

    const double lambda = pparams->lambda;
    const double eps = pparams->eps;
    const double omega = pparams->omega;
    const double dt = pparams->time_step_size;

    calculateCurvature(plinked_curve, pscheme_data);
    normal_velocity(pimage, pedge, plinked_curve, pscheme_data, pparams->get_velocity, pparams->get_g2, eps, lambda);
    tang_velocity(plinked_curve, pscheme_data, omega);

    if (!semiCoefficients(plinked_curve, pscheme_data, eps, dt))
    {
        return false;
    }

    ////////////////////    X component ///////////////////////////////////////////////////////////
    LinkedPoint* current_point = plinked_curve->first_point;

    for (size_t i = 1; i <= plinked_curve->number_of_points; i++)
    {
        pscheme_data[i].ps = pscheme_data[i].m * current_point->x +
            0.5 * (pscheme_data[i].beta_ps_expl * (current_point->next->y - current_point->previous->y)) +
            0.25 * (fmin(-pscheme_data[i].alfa, 0.0) * (current_point->previous->x - current_point->next->x) +
                fmin(pscheme_data[i].alfa, 0.0) * (current_point->next->x - current_point->previous->x));

        current_point = current_point->next;
    }

    sherman_morris(pscheme_data, plinked_curve->number_of_points);

    ///////////////////    Y component   ///////////////////////////////////////////////////////////

    current_point = plinked_curve->first_point;
    for (size_t i = 1; i <= plinked_curve->number_of_points; i++)
    {
        pscheme_data[i].ps = pscheme_data[i].m * current_point->y -
            0.5 * (pscheme_data[i].beta_ps_expl * (current_point->next->x - current_point->previous->x)) +
            0.25 * (fmin(-pscheme_data[i].alfa, 0.0) * (current_point->previous->y - current_point->next->y) +
                fmin(pscheme_data[i].alfa, 0.0) * (current_point->next->y - current_point->previous->y));

        current_point = current_point->next;
    }

    current_point = plinked_curve->first_point;
    for (size_t i = 1; i <= plinked_curve->number_of_points; i++)
    {
        updatePoint(plinked_curve, current_point, pscheme_data[i].sol, current_point->y);
        current_point = current_point->next;
    }

    sherman_morris(pscheme_data, plinked_curve->number_of_points);

    current_point = plinked_curve->first_point;
    for (size_t i = 1; i <= plinked_curve->number_of_points; i++)
    {
        updatePoint(plinked_curve, current_point, current_point->x, pscheme_data[i].sol);
        current_point = current_point->next;
    }

    return true;
}


//beta preparation
void normal_velocity(Image_Data2D* pimage, Image_Data2D* pedge, LinkedCurve* plinked_curve, SchemeData* pscheme_data,
    void(*pget_velocity)(Image_Data2D*, double, double, double*, double*),
    void(*pget_g2)(Image_Data2D*, double, double, double*),
    const double eps, const double lambda)
{
    if (plinked_curve == NULL ||
        pscheme_data == NULL ||
        pget_velocity == NULL ||
        pget_g2 == NULL ||
        pimage == NULL ||
        pedge == NULL)
    {
        return;
    }

    const size_t number_of_points = plinked_curve->number_of_points;

    double  m_pdvx_ij = 0.0;
    double  m_pdvy_ij = 0.0;
    double m_pdg2_ij = 1.0;
    LinkedPoint* current_point = plinked_curve->first_point;

    //h_i: |i, i-1|
    double h_i = -1;
    double h_i_plus = -1;

    for (size_t i = 1; i <= number_of_points; i++)
    {
        h_i = current_point->previous->distance_to_next;
        h_i_plus = current_point->distance_to_next;

        (*pget_velocity)(pedge, current_point->x, current_point->y, &m_pdvx_ij, &m_pdvy_ij);

        (*pget_g2)(pimage, current_point->x, current_point->y, &m_pdg2_ij);

        pscheme_data[i].f = m_pdvx_ij * (current_point->next->y - current_point->previous->y) / (h_i_plus + h_i) -
            m_pdvy_ij * (current_point->next->x - current_point->previous->x) / (h_i_plus + h_i);

        pscheme_data[i].f = (1.0 - lambda) * m_pdg2_ij + lambda * pscheme_data[i].f;

        current_point = current_point->next;
    }

    for (size_t i = 1; i <= number_of_points; i++)
    {
        pscheme_data[i].beta_ps_expl = pscheme_data[i].f; //beta v pravej strane pre expl.cast//lam2 == 0
        pscheme_data[i].beta = pscheme_data[i].curvature * eps - pscheme_data[i].beta_ps_expl;//celkova beta
    }

    pscheme_data[0].f = pscheme_data[number_of_points].f;
    pscheme_data[number_of_points + 1].f = pscheme_data[1].f;
    pscheme_data[0].beta = pscheme_data[number_of_points].beta;
    pscheme_data[number_of_points + 1].beta = pscheme_data[1].beta;
    pscheme_data[0].beta_ps_expl = pscheme_data[number_of_points].beta_ps_expl;
    pscheme_data[number_of_points + 1].beta_ps_expl = pscheme_data[1].beta_ps_expl;
}

// alpha
void tang_velocity(LinkedCurve* plinked_curve, SchemeData* pscheme_data,
    const double omega)
{
    int i;
    double mean = 0;
    const size_t number_of_points = plinked_curve->number_of_points;
    const double curve_length = plinked_curve->length;
    LinkedPoint* current_point = plinked_curve->first_point;
    const double avg_length = curve_length / number_of_points;
    double h_i = -1;

    for (i = 1; i <= number_of_points; i++)
    {
        h_i = current_point->previous->distance_to_next;

        mean += pscheme_data[i].beta * pscheme_data[i].curvature * h_i;

        current_point = current_point->next;
    }
    mean /= curve_length;

    pscheme_data[1].alfa = 0.0;	// alfa prveho bodu v poradi - ten sa teda nebude hybat v tangencialnom smere

    double alpha_sum = 0;

    current_point = plinked_curve->first_point->next;
    for (i = 2; i <= number_of_points; i++)
    {
        h_i = current_point->previous->distance_to_next;

        pscheme_data[i].alfa = pscheme_data[i - 1].alfa +
            pscheme_data[i].beta * pscheme_data[i].curvature * h_i -
            mean * h_i + omega * (avg_length - h_i);

        alpha_sum += pscheme_data[i].alfa;

        current_point = current_point->next;
    }
    pscheme_data[0].alfa = pscheme_data[number_of_points].alfa;
    pscheme_data[number_of_points + 1].alfa = pscheme_data[1].alfa;
}

//coeffs
bool semiCoefficients(LinkedCurve* plinked_curve, SchemeData* pscheme_data,
    const double eps, const double dt)
{
    if (plinked_curve == NULL ||
        pscheme_data == NULL) {
        return false;
    }

    double h_i = -1;
    double h_i_plus = -1;
    LinkedPoint* current_point = plinked_curve->first_point;

    for (size_t i = 1; i <= plinked_curve->number_of_points; i++)
    {
        h_i = current_point->previous->distance_to_next;
        h_i_plus = current_point->distance_to_next;

        pscheme_data[i].b = -0.5 * fmax(-pscheme_data[i].alfa, 0.0) - 1.0 / h_i * eps;//lower diagonal
        pscheme_data[i].c = -0.5 * fmax(pscheme_data[i].alfa, 0.0) - 1.0 / h_i_plus * eps;//upper diagonal

        pscheme_data[i].m = (h_i_plus + h_i) / (2.0 * dt);
        pscheme_data[i].a = pscheme_data[i].m - (pscheme_data[i].b + pscheme_data[i].c);//stiffness matrix

        if (fabs(pscheme_data[i].b) + fabs(pscheme_data[i].c) > fabs(pscheme_data[i].a) || pscheme_data[i].a < 0) {
            //the matrix is not positive dominant
            return false;
        }

        current_point = current_point->next;
    }

    return true;
}

//vypocita sa mean curvature z troch susednych elementov - 4 bodov
void calculateCurvature(LinkedCurve* plinked_curve, SchemeData* pscheme_data)
{
    if (plinked_curve == NULL ||
        pscheme_data == NULL) {
        return;
    }

    const size_t curve_length = plinked_curve->number_of_points;
    double phi;
    LinkedPoint* current_point = plinked_curve->first_point;

    double h_i_minus = -1;
    double h_i_plus = -1;
    double h_i = -1;

    double x_i_minus_2, x_i_minus_1, x_i, x_i_plus_1;
    double y_i_minus_2, y_i_minus_1, y_i, y_i_plus_1;

    for (size_t i = 1; i <= curve_length; i++)
    {
        h_i_minus = current_point->previous->previous->distance_to_next;
        h_i_plus = current_point->distance_to_next;
        h_i = current_point->previous->distance_to_next;

        x_i_minus_2 = current_point->previous->previous->x;
        x_i_minus_1 = current_point->previous->x;
        x_i = current_point->x;
        x_i_plus_1 = current_point->next->x;

        y_i_minus_2 = current_point->previous->previous->y;
        y_i_minus_1 = current_point->previous->y;
        y_i = current_point->y;
        y_i_plus_1 = current_point->next->y;

        phi = (
            (x_i_minus_1 - x_i_minus_2) * (x_i_plus_1 - x_i) +
            (y_i_minus_1 - y_i_minus_2) * (y_i_plus_1 - y_i)
            ) / (h_i_plus * h_i_minus);

        if (phi > 1)
            phi = 1.;
        else if (phi < -1)
            phi = -1.;

        pscheme_data[i].curvature = 0.5 * signum(
            (x_i_minus_1 - x_i_minus_2) * (y_i_plus_1 - y_i) -
            (x_i_plus_1 - x_i) * (y_i_minus_1 - y_i_minus_2)
        ) * acos(phi) / h_i;

        current_point = current_point->next;
    }

    pscheme_data[0].curvature = pscheme_data[curve_length].curvature;
    pscheme_data[curve_length + 1].curvature = pscheme_data[curve_length + 1].curvature;
}
