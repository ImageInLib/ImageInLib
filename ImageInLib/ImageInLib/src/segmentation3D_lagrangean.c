#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "common_functions.h"
#include "file.h"
#include "segmentation2D_lagrangean.h"
#include "segmentation3D_lagrangean.h"
#include "solvers.h"

/// <summary>
/// This function evolve by one time step a given curve in vector field
/// </summary>
/// <param name="pDistanceMap">distance map used to construct the vector field</param>
/// <param name="plinked_curve">linked curve (curve points)</param>
/// <param name="pscheme_data">coefficients of the linear system</param>
/// <param name="pparams">function parameters</param>
/// <returns></returns>
bool evolveBySingleStep3D(Image_Data* pDistanceMap, LinkedCurve3D* plinked_curve, SchemeData3D* pscheme_data, const Lagrangean3DSegmentationParameters* pparams);

/// <summary>
/// Function to get the velocity field from the distance map
/// It computes the 3D gradient in given point
/// </summary>
/// <param name="pDistanceMap">distance map</param>
/// <param name="x">x component of the point</param>
/// <param name="y">y component of the point</param>
/// <param name="z">z component of the point</param>
/// <param name="vx">gradient x component of the point</param>
/// <param name="vy">gradient y component of the point</param>
/// <param name="vz">gradient z component of the point</param>
void getVelocity3D(Image_Data* pDistanceMap, const double x, const double y, const double z, double* vx, double* vy, double* vz);

/// <summary>
/// This function compute the normal velocity from the distance map
/// </summary>
/// <param name="pDistanceMap">distance map</param>
/// <param name="plinked_curve">linked curve, each point will be allocated with the corresponding velocity</param>
/// <param name="pscheme_data">scheme data</param>
/// <param name="eps">parameter epsilon</param>
/// <param name="mu">parameter mu</param>
void normal_velocity3D(Image_Data* pDistanceMap, LinkedCurve3D* plinked_curve, SchemeData3D* pscheme_data,
    const double eps, const double mu);

/// <summary>
/// This function compute the tangential velocity
/// </summary>
/// <param name="plinked_curve">linked curve, each point will be allocated with the corresponding velocity</param>
/// <param name="pscheme_data">scheme data</param>
/// <param name="omega">parameter omega</param>
void tang_velocity3D(LinkedCurve3D* plinked_curve, SchemeData3D* pscheme_data, const double omega);

/// <summary>
/// This function set the coefficients needed to solve the linear system
/// </summary>
/// <param name="plinked_curve">linked curve</param>
/// <param name="pscheme_data">scheme data</param>
/// <param name="eps">parameter epsilon</param>
/// <param name="dt">parameter of time step</param>
/// <returns></returns>
bool semiCoefficients3D(LinkedCurve3D* plinked_curve, SchemeData3D* pscheme_data, const double eps, const double dt);

//===========================================================

bool lagrangeanExplicit3DCurveSegmentation(Image_Data inputImage3D, const Lagrangean3DSegmentationParameters* pSegmentationParams,
    unsigned char* pOutputPathPtr, Curve3D* pResultSegmentation)
{

    if (pSegmentationParams == NULL || pResultSegmentation == NULL) {
        return false;
    }

    CurvePoint3D* pOldCurve = (CurvePoint3D*)malloc(sizeof(CurvePoint3D) * pResultSegmentation->numPoints);
    Curve3D oldSegmentation = { pOldCurve, pResultSegmentation->numPoints };
    
    const size_t sizeHeight = sizeof(dataType*) * inputImage3D.height;
    const size_t sliceDimension = inputImage3D.length * inputImage3D.width;
    const size_t sliceSize = sliceDimension * sizeof(dataType);

    dataType** pgrad_x = (dataType**)malloc(sizeHeight);       // velocity component x
    dataType** pgrad_y = (dataType**)malloc(sizeHeight);       // velocity component y
    dataType** pgrad_z = (dataType**)malloc(sizeHeight);       // velocity component z
    dataType** abs_val_grad = (dataType**)malloc(sizeHeight);  // absolute value of gradient
    dataType** edge_detector = (dataType**)malloc(sizeHeight); // edge detector
    for (size_t k = 0; k < inputImage3D.height; k++) {
        pgrad_x[k] = (dataType*)malloc(sliceSize);
        pgrad_y[k] = (dataType*)malloc(sliceSize);
        pgrad_z[k] = (dataType*)malloc(sliceSize);
        abs_val_grad[k] = (dataType*)malloc(sliceSize);
        edge_detector[k] = (dataType*)malloc(sliceSize);
        if (pgrad_x[k] == NULL || pgrad_y[k] == NULL || pgrad_z[k] == NULL ||
            abs_val_grad[k] == NULL || edge_detector[k] == NULL)
            return false;
    }
    if (pgrad_x == NULL || pgrad_y == NULL || pgrad_z == NULL ||
        abs_val_grad == NULL || edge_detector == NULL)
        return false;

    dataType** similar_intensity_detector = NULL;
    if (!pSegmentationParams->open_curve)
    {
        similar_intensity_detector = (dataType**)malloc(sizeHeight); // similar intensity detector
        for (size_t k = 0; k < inputImage3D.height; k++) {
            similar_intensity_detector[k] = (dataType*)malloc(sliceSize);
        }
    }

    const dataType hx = 1.0, hy = 1.0, hz = 1.0;                 //spatial discretization step
    const dataType hx_c = 1.0, hy_c = 1.0, hz_c = 1.0;           //h for central differences
    Point3D current_grad;
    FiniteVolumeSize3D finite_volume = { 1.0, 1.0, 1.0 };

    dataType lambda;
    if (!pSegmentationParams->open_curve) {
        lambda = pSegmentationParams->lambda;
    }
    else {
        lambda = 1.0;
    }

    for (size_t k = 0; k < inputImage3D.height; k++)
    {
        for (size_t i = 0; i < inputImage3D.length; i++)
        {
            for (size_t j = 0; j < inputImage3D.width; j++)
            {
                getGradient3D(inputImage3D.imageDataPtr, inputImage3D.width, inputImage3D.length, inputImage3D.height, j, i, k, finite_volume, &current_grad);
                size_t xd = x_new(j, i, inputImage3D.width);

                //Absolute gradient
                abs_val_grad[k][xd] = norm3D(current_grad);

                //Edge detector
                edge_detector[k][xd] = edgeDetector(abs_val_grad[k][xd], pSegmentationParams->edge_detector_coef);

                //Similar intensity detector
                if (!pSegmentationParams->open_curve) {
                    similar_intensity_detector[k][xd] = similarIntensityDetector(inputImage3D.imageDataPtr[k][xd], pSegmentationParams->reference_intensity, pSegmentationParams->intensity_coef);
                }
            }
        }
    }

    //get the velocity field
    for (size_t k = 0; k < inputImage3D.height; k++)
    {
        for (size_t i = 0; i < inputImage3D.length; i++)
        {
            for (size_t j = 0; j < inputImage3D.width; j++)
            {
                getGradient3D(edge_detector, inputImage3D.width, inputImage3D.length, inputImage3D.height, j, i, k, finite_volume, &current_grad);
                size_t xd = x_new(j, i, inputImage3D.width);
                pgrad_x[k][xd] = current_grad.x;
                pgrad_y[k][xd] = current_grad.y;
                pgrad_z[k][xd] = current_grad.z;
            }
        }
    }

    size_t current_i = 0;
    size_t current_j = 0;
    size_t current_k = 0;

    for (size_t i = 0; i < (pResultSegmentation->numPoints); i++)
    {
        pResultSegmentation->pPoints[i].x = pSegmentationParams->pinitial_condition->pPoints[i].x;
        pResultSegmentation->pPoints[i].y = pSegmentationParams->pinitial_condition->pPoints[i].y;
        pResultSegmentation->pPoints[i].z = pSegmentationParams->pinitial_condition->pPoints[i].z;
    }

    size_t iterPt = 0;
    size_t maxIterPt = oldSegmentation.numPoints;

    if (pSegmentationParams->open_curve) {
        iterPt = 1;
        maxIterPt = oldSegmentation.numPoints - 1;
    }

    //_l - lower, _g - greater, _c - current
    dataType vx, vy, vz;
    dataType rx_l, rx_g;
    dataType ry_l, ry_g;
    dataType rz_l, rz_g;
    dataType rx_c, ry_c, rz_c;
    dataType nx, ny, nz, dot, norm, tx, ty, tz;
    dataType h_g, h_c;
    dataType curv_x, curv_y, curv_z, mean_dist;

    dataType tau = pSegmentationParams->time_step_size;
    dataType mu = pSegmentationParams->mu;
    dataType eps = pSegmentationParams->eps;

    if (!pSegmentationParams->open_curve) 
    {
        for (size_t t = 0; t < pSegmentationParams->num_time_steps; t++)
        {
            for (size_t i = 0; i < (oldSegmentation.numPoints); i++)
            {
                oldSegmentation.pPoints[i].x = pResultSegmentation->pPoints[i].x;
                oldSegmentation.pPoints[i].y = pResultSegmentation->pPoints[i].y;
                oldSegmentation.pPoints[i].z = pResultSegmentation->pPoints[i].z;
            }

            for (size_t i = iterPt; i < maxIterPt; i++)
            {
                //let us move by curve points just in vector field
                current_i = (size_t)(oldSegmentation.pPoints[i].x + 0.5);
                current_j = (size_t)(oldSegmentation.pPoints[i].y + 0.5);
                current_k = (size_t)(oldSegmentation.pPoints[i].z + 0.5);

                //let us keep points inside the image
                if (current_i < 0) {
                    current_i = 0;
                }
                else if (current_i >= inputImage3D.length) {
                    current_i = inputImage3D.length - 1;
                }

                if (current_j < 0) {
                    current_j = 0;
                }
                else if (current_j >= inputImage3D.width) {
                    current_j = inputImage3D.width - 1;
                }

                if (current_k < 0) {
                    current_k = 0;
                }
                else if (current_k >= inputImage3D.height) {
                    current_k = inputImage3D.height - 1;
                }

                size_t xd = x_new(current_j, current_i, inputImage3D.width);

                rx_c = oldSegmentation.pPoints[i].x;
                ry_c = oldSegmentation.pPoints[i].y;
                rz_c = oldSegmentation.pPoints[i].z;

                if (i == 0) {
                    rx_l = oldSegmentation.pPoints[oldSegmentation.numPoints - 1].x;
                    ry_l = oldSegmentation.pPoints[oldSegmentation.numPoints - 1].y;
                    rz_l = oldSegmentation.pPoints[oldSegmentation.numPoints - 1].z;

                    rx_g = oldSegmentation.pPoints[1].x;
                    ry_g = oldSegmentation.pPoints[1].y;
                    rz_g = oldSegmentation.pPoints[1].z;
                }
                else if (i == oldSegmentation.numPoints - 1) {
                    rx_l = oldSegmentation.pPoints[i - 1].x;
                    ry_l = oldSegmentation.pPoints[i - 1].y;
                    rz_l = oldSegmentation.pPoints[i - 1].z;

                    rx_g = oldSegmentation.pPoints[0].x;
                    ry_g = oldSegmentation.pPoints[0].y;
                    rz_g = oldSegmentation.pPoints[0].z;
                }
                else {
                    rx_l = oldSegmentation.pPoints[i - 1].x;
                    ry_l = oldSegmentation.pPoints[i - 1].y;
                    rz_l = oldSegmentation.pPoints[i - 1].z;

                    rx_g = oldSegmentation.pPoints[i + 1].x;
                    ry_g = oldSegmentation.pPoints[i + 1].y;
                    rz_g = oldSegmentation.pPoints[i + 1].z;
                }

                h_c = (dataType)sqrt(pow(rx_c - rx_l, 2) + pow(ry_c - ry_l, 2) + pow(rz_c - rz_l, 2));
                h_g = (dataType)sqrt(pow(rx_g - rx_c, 2) + pow(ry_g - ry_c, 2) + pow(rz_g - ry_c, 2));

                norm = (dataType)sqrt(pow(rx_g - rx_l, 2) + pow(ry_g - ry_l, 2) + pow(rz_g - rz_l, 2));
                tx = (rx_g - rx_l) / norm;
                ty = (ry_g - ry_l) / norm;
                tz = (rz_g - rz_l) / norm;

                nx = ty * tz;
                ny = -2 * tx * tz;
                nz = tx * ty;

                dot = pgrad_x[current_k][xd] * nx + pgrad_y[current_k][xd] * ny + pgrad_z[current_k][xd] * nz;

                mean_dist = 2.0 / (h_g + h_c);
                curv_x = mean_dist * ((rx_g - rx_c) / h_g - ((rx_c - rx_l) / h_c));
                curv_y = mean_dist * ((ry_g - ry_c) / h_g - ((ry_c - ry_l) / h_c));
                curv_z = mean_dist * ((rz_g - rz_c) / h_g - ((rz_c - rz_l) / h_c));

                vx = mu * ((1 - lambda) * similar_intensity_detector[current_k][xd]
                    - lambda * dot) + eps * curv_x;

                vy = mu * ((1 - lambda) * similar_intensity_detector[current_k][xd]
                    - lambda * dot) + eps * curv_y;

                vz = mu * ((1 - lambda) * similar_intensity_detector[current_k][xd]
                    - lambda * dot) + eps * curv_z;

                pResultSegmentation->pPoints[i].x = oldSegmentation.pPoints[i].x + tau * vx * nx;
                pResultSegmentation->pPoints[i].y = oldSegmentation.pPoints[i].y + tau * vy * ny;
                pResultSegmentation->pPoints[i].z = oldSegmentation.pPoints[i].z + tau * vz * nz;

            }

        }
    }

    //Free memory
    free(pOldCurve);
    for (size_t k = 0; k < inputImage3D.height; k++) {
        free(pgrad_x[k]);
        free(pgrad_y[k]);
        free(pgrad_z[k]);
        free(abs_val_grad[k]);
        free(edge_detector[k]);
    }
    free(pgrad_x);
    free(pgrad_y);
    free(pgrad_z);
    free(abs_val_grad);
    free(edge_detector);

    if (!pSegmentationParams->open_curve)
    {
        for (size_t k = 0; k < inputImage3D.height; k++) {
            free(similar_intensity_detector[k]);
        }
        free(similar_intensity_detector);
    }

    return true;
}

bool lagrangeanSemiImplicit3DCurveSegmentation(Image_Data inputImage3D, const Lagrangean3DSegmentationParameters* pSegmentationParams,
    unsigned char* pOutputPathPtr, Curve3D* pResultSegmentation)
{
    //check if the pointers a well allocated
    if (pSegmentationParams == NULL || pResultSegmentation == NULL) {
        return false;
    }

    //check that we have at least three points at the biginning
    if (pSegmentationParams->num_points < 3) {
        return false;
    }

    //pOutputPathPtr : can be used to save the path during the motion
    //In the current implemention, just the final curve can be exported

    //let us consider single curve without topological changes
    resetIDGenerator();

    bool isOrientedPositively = true; //does not matter for open curves

    if (!pSegmentationParams->open_curve)
    {
        isOrientedPositively = is3dCurveOrientedPositively(pSegmentationParams->pinitial_condition);
    }

    //create initial linked curve
    LinkedCurve3D linked_curve = create3dLinkedCurve();

    //Initialize the linked curve
    initialize3dLinkedCurve(pSegmentationParams->pinitial_condition, &linked_curve, !isOrientedPositively, !pSegmentationParams->open_curve);

    size_t length_of_data = linked_curve.number_of_points + 2;
    SchemeData3D* pscheme_data = (SchemeData3D*)calloc(length_of_data, sizeof(SchemeData3D));

    size_t it = 1;
    double motion = 0.0;
    double distance_to_next = 0.0;
    do {

        if (length_of_data < linked_curve.number_of_points + 2)
        {
            free(pscheme_data);
            length_of_data = linked_curve.number_of_points + 2;
            pscheme_data = (SchemeData3D*)calloc(length_of_data, sizeof(SchemeData3D));
        }
        
        //evolve curve
        evolveBySingleStep3D(&inputImage3D, &linked_curve, pscheme_data, pSegmentationParams);

        it++;

    } while (it < pSegmentationParams->num_time_steps);

    LinkedPoint3D* pPoint = linked_curve.first_point;

    for (size_t i = 0; i < linked_curve.number_of_points; i++)
    {
        pResultSegmentation->pPoints[i].x = (dataType)pPoint->x;
        pResultSegmentation->pPoints[i].y = (dataType)pPoint->y;
        pResultSegmentation->pPoints[i].z = (dataType)pPoint->z;
        pPoint = pPoint->next;
    }

    free(pscheme_data);
    release3dLinkedCurve(&linked_curve);

    return true;
}

bool evolveBySingleStep3D(Image_Data* pDistanceMap, LinkedCurve3D* plinked_curve, SchemeData3D* pscheme_data, const Lagrangean3DSegmentationParameters* pparams)
{
    //check if the pointers a well allocated
    if (plinked_curve == NULL || pscheme_data == NULL || pparams == NULL ||
        pDistanceMap->imageDataPtr == NULL)
    {
        return false;
    }

    const double eps = pparams->eps;
    const double omega = pparams->omega;
    const double dt = pparams->time_step_size;
    const double mu = pparams->mu;

    //function to compute the normal velocity
    normal_velocity3D(pDistanceMap, plinked_curve, pscheme_data, eps, mu);

    //function to compute the tangential velocity
    tang_velocity3D(plinked_curve, pscheme_data, omega);

    if (!semiCoefficients3D(plinked_curve, pscheme_data, eps, dt))
    {
        return false;
    }

    LinkedPoint3D* current_point = plinked_curve->first_point;

    //////////////////////    X component ///////////////////////////////////////////////////////////
  
    for (size_t i = 1; i <= plinked_curve->number_of_points; i++)
    {

        if (pparams->open_curve && (i == 1 || i == plinked_curve->number_of_points))
        {
            pscheme_data[i].ps = current_point->x;
        }
        else
        {
            pscheme_data[i].ps = pscheme_data[i].m * current_point->x + mu * current_point->nvx * dt * pscheme_data[i].m;
        }
        current_point = current_point->next;
    }

    if (pparams->open_curve)
    {
        calculate_by_thomas3D(pscheme_data, plinked_curve->number_of_points);
    }
    else
    {
        sherman_morris3D(pscheme_data, plinked_curve->number_of_points);
    }

    /////////////////////    Y component   ///////////////////////////////////////////////////////////

    current_point = plinked_curve->first_point;
    for (size_t i = 1; i <= plinked_curve->number_of_points; i++)
    {

        if (pparams->open_curve && (i == 1 || i == plinked_curve->number_of_points))
        {
            pscheme_data[i].ps = current_point->y;
        }
        else
        {
            pscheme_data[i].ps = pscheme_data[i].m * current_point->y + mu * current_point->nvy * dt * pscheme_data[i].m;
        }
        current_point = current_point->next;
    }

    current_point = plinked_curve->first_point;
    for (size_t i = 1; i <= plinked_curve->number_of_points; i++)
    {
        update3dPoint(plinked_curve, current_point, pscheme_data[i].sol, current_point->y, current_point->z);
        current_point = current_point->next;
    }

    if (pparams->open_curve)
    {
        calculate_by_thomas3D(pscheme_data, plinked_curve->number_of_points);
    }
    else
    {
        sherman_morris3D(pscheme_data, plinked_curve->number_of_points);
    }

    /////////////////////    Z component   ///////////////////////////////////////////////////////////

    current_point = plinked_curve->first_point;
    for (size_t i = 1; i <= plinked_curve->number_of_points; i++)
    {

        if (pparams->open_curve && (i == 1 || i == plinked_curve->number_of_points))
        {
            pscheme_data[i].ps = current_point->z;
        }
        else
        {
            pscheme_data[i].ps = pscheme_data[i].m * current_point->z + mu * current_point->nvz * dt * pscheme_data[i].m;
        }
        current_point = current_point->next;
    }

    current_point = plinked_curve->first_point;
    for (size_t i = 1; i <= plinked_curve->number_of_points; i++)
    {
        update3dPoint(plinked_curve, current_point, current_point->x, pscheme_data[i].sol, current_point->z);
        current_point = current_point->next;
    }

    if (pparams->open_curve)
    {
        calculate_by_thomas3D(pscheme_data, plinked_curve->number_of_points);
    }
    else
    {
        sherman_morris3D(pscheme_data, plinked_curve->number_of_points);
    }

    current_point = plinked_curve->first_point;
    for (size_t i = 1; i <= plinked_curve->number_of_points; i++)
    {
        update3dPoint(plinked_curve, current_point, current_point->x, current_point->y, pscheme_data[i].sol);
        current_point = current_point->next;
    }

    return true;
}

void getVelocity3D(Image_Data* pDistanceMap, const double x, const double y, const double z, double* vx, double* vy, double* vz)
{

    //get the gradient components/coordinates
    const size_t x_dis = (size_t)x;
    const size_t y_dis = (size_t)y;
    const size_t z_dis = (size_t)z;
    size_t xd = x_new(x_dis, y_dis, pDistanceMap->width);
    Point3D current_grad;
    const FiniteVolumeSize3D finite_volume = { pDistanceMap->spacing.sx, pDistanceMap->spacing.sy, pDistanceMap->spacing.sz };

    getGradient3D(pDistanceMap->imageDataPtr, pDistanceMap->length, pDistanceMap->width, pDistanceMap->height, x_dis, y_dis, z_dis, finite_volume, &current_grad);

    *vx = current_grad.x;
    *vy = current_grad.y;
    *vz = current_grad.z;

}

void normal_velocity3D(Image_Data* pDistanceMap, LinkedCurve3D* plinked_curve, SchemeData3D* pscheme_data,
    const double eps, const double mu)
{

    //check if the pointers a well allocated
    if (plinked_curve == NULL || pscheme_data == NULL)
    {
        return;
    }

    const size_t number_of_points = plinked_curve->number_of_points;

    double h_i = -1, h_i_plus = -1;            // distance between two neighboring points
    double vx = 0, vy = 0, vz = 0;             // external velocity field components 
    double tx = 0, ty = 0, tz = 0;             //tangential vector components
    double som_dist = 0, dot = 0, norm_nv = 1.0; 
    double nvx = 0, nvy = 0, nvz = 0;          //normal velocity vector components
    double n1_x = 0, n1_y = 0, n1_z = 0;       //normal plan vector components
    double n2_x = 0, n2_y = 0, n2_z = 0;       //normal plan vector components
    double curv_x = 0, curv_y = 0, curv_z = 0; //discrete curvature vector components

    bool is_closed_curve = plinked_curve->first_point->previous != NULL;

    LinkedPoint3D* current_point = plinked_curve->first_point;

    for (size_t i = 1; i <= number_of_points; i++)
    {
        if (is_closed_curve || (i > 1 && i < plinked_curve->number_of_points))
        {
            h_i = current_point->previous->distance_to_next;
            h_i_plus = current_point->distance_to_next;
            som_dist = h_i_plus + h_i;

            //get the external velocity field
            getVelocity3D(pDistanceMap, current_point->x, current_point->y, current_point->z, &vx, &vy, &vz);

            tx = (current_point->next->x - current_point->previous->x) / som_dist;
            ty = (current_point->next->y - current_point->previous->y) / som_dist;
            tz = (current_point->next->z - current_point->previous->z) / som_dist;

            dot = tx * vx + ty * vy + tz * vz;

            current_point->nvx = vx - dot * tx;
            current_point->nvy = vy - dot * ty;
            current_point->nvz = vz - dot * tz;

            Point3D pnorm = { current_point->nvx, current_point->nvy, current_point->nvz };
            norm_nv = norm3D(pnorm);
            if (norm_nv == 0) {
                norm_nv = 1.0;
            }

            n1_x = nvx / norm_nv;
            n1_y = nvy / norm_nv;
            n1_z = nvz / norm_nv;

            n2_x = n1_y * tz - n1_z * ty;
            n2_y = n1_z * tx - n1_x * tz;
            n2_z = n1_x * ty - n1_y * tx;

            curv_x = (2.0 / som_dist) * (((current_point->next->x - current_point->x) / h_i_plus) - ((current_point->x - current_point->previous->x) / h_i));
            curv_y = (2.0 / som_dist) * (((current_point->next->y - current_point->y) / h_i_plus) - ((current_point->y - current_point->previous->y) / h_i));
            curv_z = (2.0 / som_dist) * (((current_point->next->z - current_point->z) / h_i_plus) - ((current_point->z - current_point->previous->z) / h_i));

            pscheme_data[i].k1 = curv_x * n1_x + curv_y * n1_y + curv_z * n1_z;
            pscheme_data[i].k2 = curv_x * n2_x + curv_y * n2_y + curv_z * n2_z;
            pscheme_data[i].u = eps * pscheme_data[i].k1 + mu * norm_nv;
            pscheme_data[i].v = eps * pscheme_data[i].k2;

        }
        else
        {
            pscheme_data[i].k1 = 0.0;
            pscheme_data[i].k2 = 0.0;
            pscheme_data[i].u = 0.0;
            pscheme_data[i].v = 0.0;
            current_point->nvx = 0.0;
            current_point->nvy = 0.0;
            current_point->nvz = 0.0;
        }

        current_point = current_point->next;
    }

    if (is_closed_curve)
    {
        pscheme_data[0].k1 = pscheme_data[number_of_points].k1;
        pscheme_data[number_of_points + 1].k1 = pscheme_data[1].k1;
        pscheme_data[0].k2 = pscheme_data[number_of_points].k2;
        pscheme_data[number_of_points + 1].k2 = pscheme_data[1].k2;
        pscheme_data[0].u = pscheme_data[number_of_points].u;
        pscheme_data[number_of_points + 1].u = pscheme_data[1].u;
        pscheme_data[0].v = pscheme_data[number_of_points].v;
        pscheme_data[number_of_points + 1].v = pscheme_data[1].v;
    }

}

void tang_velocity3D(LinkedCurve3D* plinked_curve, SchemeData3D* pscheme_data, const double omega)
{
    
    double mean = 0.0;
    const size_t number_of_points = plinked_curve->number_of_points;
    const double curve_length = plinked_curve->length;
    LinkedPoint3D* current_point = plinked_curve->first_point;
    double h_i = -1;

    bool is_curve_closed = plinked_curve->first_point->previous != NULL;

    double avg_length = 1.0;
    if (is_curve_closed) {
        avg_length = curve_length / (double)number_of_points;
    }
    else {
        avg_length = curve_length / (double)(number_of_points - 1);
    }
    

    for (size_t i = 1; i <= number_of_points; i++)
    {
        if (is_curve_closed || i > 1)
        {
            h_i = current_point->previous->distance_to_next;
        }
        else
        {
            h_i = current_point->distance_to_next;
        }

        mean += h_i * (pscheme_data[i].u * pscheme_data[i].k1 + pscheme_data[i].v * pscheme_data[i].k2);

        current_point = current_point->next;
    }

    mean /= curve_length;

    // the alpha of the first point in the sequence - it will therefore not move in the tangential direction
    pscheme_data[0].alfa = 0.0;
    pscheme_data[1].alfa = 0.0;

    current_point = plinked_curve->first_point;
    for (size_t i = 1; i <= number_of_points; i++)
    {
        if (i == 1) {
            h_i = current_point->distance_to_next;
        }
        else
        {
            h_i = current_point->previous->distance_to_next;
        }
        
        pscheme_data[i].alfa = pscheme_data[i - 1].alfa + h_i * (pscheme_data[i].u * pscheme_data[i].k1 + pscheme_data[i].v * pscheme_data[i].k2) - h_i * mean + omega * (avg_length - h_i);
        current_point = current_point->next;
    }

    pscheme_data[number_of_points].alfa = 0.0;
    pscheme_data[number_of_points + 1].alfa = 0.0;

    if (is_curve_closed)
    {
        pscheme_data[0].alfa = pscheme_data[number_of_points].alfa;
        pscheme_data[number_of_points + 1].alfa = pscheme_data[1].alfa;
    }

}

bool semiCoefficients3D(LinkedCurve3D* plinked_curve, SchemeData3D* pscheme_data, const double eps, const double dt)
{
    if (plinked_curve == NULL || pscheme_data == NULL)
    {
        return false;
    }

    double h_i = -1;
    double h_i_plus = -1;
    LinkedPoint3D* current_point = plinked_curve->first_point;
    LinkedPoint3D* previous_point;

    bool is_curve_closed = plinked_curve->first_point->previous != NULL;

    for (size_t iter = 1; iter <= plinked_curve->number_of_points; iter++)
    {

        if (is_curve_closed || (iter > 1 && iter < plinked_curve->number_of_points))
        {
            previous_point = current_point->previous;
            h_i = previous_point->distance_to_next;

            if (is_curve_closed || iter < plinked_curve->number_of_points)
            {
                h_i_plus = current_point->distance_to_next;
            }
            else
            {
                h_i_plus = h_i;
            }

            pscheme_data[iter].a = -eps / h_i + 0.5 * pscheme_data[iter].alfa;      //lower diagonal
            pscheme_data[iter].c = -eps / h_i_plus - 0.5 * pscheme_data[iter].alfa; //upper diagonal
            pscheme_data[iter].m = (h_i_plus + h_i) / (2.0 * dt);
            pscheme_data[iter].b = pscheme_data[iter].m - (pscheme_data[iter].a + pscheme_data[iter].c);//diagonal

            if ((fabs(pscheme_data[iter].a) + fabs(pscheme_data[iter].c)) > fabs(pscheme_data[iter].b) || pscheme_data[iter].b < 0) {
                return false;
            }

        }
        else
        {
            if (iter == 1)
            {
                h_i = current_point->distance_to_next;
                h_i_plus = h_i;
            }
            else
            {
                h_i = current_point->previous->distance_to_next;
                h_i_plus = h_i;
            }
          
            pscheme_data[iter].a = 0.0;
            pscheme_data[iter].c = 0.0;
            pscheme_data[iter].m = (h_i_plus + h_i) / (2.0 * dt);
            pscheme_data[iter].b = 1.0;
        }
      
        current_point = current_point->next;
    }

    return true;
}