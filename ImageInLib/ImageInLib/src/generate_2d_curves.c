#include "generate_2d_curves.h"
#include "common_math.h"
#include "common_functions.h"
#include <stdlib.h>

bool generateCircleCurve(Curve2D* pcurve, const Point2D* pinitial_points, const size_t initial_points_count, const double points_distance)
{
    if (pcurve == NULL || pinitial_points == NULL)
    {
        return false;
    }

    const size_t circlePointsCount = howManyPointsForCircleCurve(pinitial_points, initial_points_count, points_distance);


    CurvePoint2D* pcircle_points = malloc(sizeof(CurvePoint2D) * circlePointsCount);
    pcurve->pPoints = pcircle_points;
    pcurve->numPoints = circlePointsCount;

    const double radius = getPoint2DDistance(pinitial_points[0], pinitial_points[1]);

    double anlgleStep = 2 * M_PI / (double)circlePointsCount;
    Point2D center = pinitial_points[0];

    for (size_t i = 0; i < circlePointsCount; i++)
    {
        double phi = anlgleStep * (double)i;
        pcurve->pPoints[i].x = (dataType)(center.x + radius * cos(phi));
        pcurve->pPoints[i].y = (dataType)(center.y + radius * sin(phi));
        pcurve->pPoints[i].isEndPoint = false;
    }

    return true;
}

size_t howManyPointsForCircleCurve(const Point2D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance)
{
    if (initialPointsCount != 2)
    {
        return 0;
    }

    const double radius = getPoint2DDistance(pInitialPoints[0], pInitialPoints[1]);
    const double perimeter = getCirclePerimeter(radius);
    return (size_t)((perimeter / pointsDistance) + 0.5);
}

double getCirclePerimeter(const double radius)
{
    if (radius < 0)
    {
        return 0;
    }

    return 2 * M_PI * radius;
}

bool generateStraightLineCurve(Curve2D* pCurve, const Point2D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance)
{
    if(pInitialPoints == NULL || pCurve  == NULL)
    {
        return false;
    }

    const size_t line_points_count = howManyPointsForStraightLineCurve(pInitialPoints, initialPointsCount, pointsDistance);

    CurvePoint2D* p_line_points = (CurvePoint2D*)malloc(sizeof(CurvePoint2D) * line_points_count);

    if (p_line_points == NULL) {
        return false;
    }

    pCurve->pPoints = p_line_points;
    pCurve->numPoints = line_points_count;

    const double length = getPoint2DDistance(pInitialPoints[0], pInitialPoints[1]);

    double line_step = length / (double)(line_points_count - 1);
    Point2D first_point = pInitialPoints[0];
    Point2D last_point = pInitialPoints[1];

    Point2D tangentVector = { (dataType)((last_point.x - first_point.x)/ length), (dataType)((last_point.y - first_point.y)/ length) };

    for (size_t i = 0; i < line_points_count; i++)
    {
        double parameter = line_step * (double)i;
        pCurve->pPoints[i].x = (dataType)(first_point.x + parameter * tangentVector.x);
        pCurve->pPoints[i].y = (dataType)(first_point.y + parameter * tangentVector.y);
        
        if (i == 0 || i == line_points_count - 1)
        {
            pCurve->pPoints[i].isEndPoint = true;
        }
        else
        {
            pCurve->pPoints[i].isEndPoint = false;
        }
    }

    return true;
}

size_t howManyPointsForStraightLineCurve(const Point2D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance)
{
    if (initialPointsCount != 2)
    {
        return 0;
    }

    const double length = getPoint2DDistance(pInitialPoints[0], pInitialPoints[1]);
    return (size_t)((length / pointsDistance) + 1.5);
}
