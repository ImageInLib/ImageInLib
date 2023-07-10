#include "generate_2d_curves.h"
#include "common_math.h"

bool generateCircleCurve(Point2D* pCurve, const size_t circlePointsCount, const Point2D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance)
{
    if (howManyPointsForCircleCurve(pInitialPoints, initialPointsCount, pointsDistance) != circlePointsCount)
    {
        return false;
    }

    const double radius = getPoint2DDistance(pInitialPoints[0], pInitialPoints[1]);

    double anlgleStep = 2 * M_PI / (double)circlePointsCount;
    Point2D center = pInitialPoints[0];

    for (size_t i = 0; i < circlePointsCount; i++)
    {
        double phi = anlgleStep * (double)i;
        pCurve[i].x = center.x + radius * cos(phi);
        pCurve[i].y = center.y + radius * sin(phi);
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
