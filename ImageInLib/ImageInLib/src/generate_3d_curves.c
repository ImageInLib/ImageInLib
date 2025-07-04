#include <stdlib.h>
#include "common_math.h"
#include "common_functions.h"
#include "generate_2d_curves.h"
#include "generate_3d_curves.h"

double getSphereSurfaceArea(const double radius)
{
    if (radius < 0)
    {
        return 0;
    }
    return 4 * M_PI * radius * radius;
}

size_t howManyPointsFor3DStraightLineCurve(const Point3D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance)
{
    if (initialPointsCount != 2)
    {
        return 0;
    }
    const double curveLength = getPoint3DDistance(pInitialPoints[0], pInitialPoints[1]);
    return (size_t)((curveLength / pointsDistance) + 1.5);
}

size_t howManyPointsForSphereCurve(const Point3D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance)
{
    if (initialPointsCount != 2)
    {
        return 0;
    }
    double radius = getPoint3DDistance(pInitialPoints[0], pInitialPoints[1]);
    double sphereSurface = getSphereSurfaceArea(radius);
    double singlePointSurface = pointsDistance * pointsDistance;
    return (size_t)((sphereSurface / singlePointSurface) + 0.5);
}

bool generate3DLineCurve(Curve3D* pCurve, const Point3D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance)
{
    //TODO
    return true;
}

bool generateSphereCurve(Curve3D* pCurve, const Point3D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance)
{
    if (pCurve == NULL || pInitialPoints == NULL)
    {
        return false;
    }

    size_t numberOfPoints = howManyPointsForSphereCurve(pInitialPoints, initialPointsCount, pointsDistance);
    CurvePoint3D* pcircle_points = malloc(sizeof(CurvePoint3D) * numberOfPoints);
    pCurve->pPoints = pcircle_points;
    pCurve->numPoints = numberOfPoints;
    const double radius = getPoint3DDistance(pInitialPoints[0], pInitialPoints[1]);

    Point3D center = pInitialPoints[0];
    
    double phi = M_PI * (3 - sqrt(5)); //golden angle
    for (size_t i = 0; i < numberOfPoints; i++) {
        dataType z = 1 - ((double)i / (double)(numberOfPoints - 1)) * 2;
        double altitude = sqrt(1 - z * z);
        double theta = phi * (double)i;
        pCurve->pPoints[i].x = (dataType)(center.x + radius * altitude * cos(theta));
        pCurve->pPoints[i].y = (dataType)(center.y + radius * altitude * sin(theta));
        pCurve->pPoints[i].z = (dataType)(center.z + radius * z);
        pCurve->pPoints[i].isEndPoint = false;
    }
    return true;
}
