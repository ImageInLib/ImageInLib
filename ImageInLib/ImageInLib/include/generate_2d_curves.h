#pragma once
#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef GENERATE2DCURVES_H
#define GENERATE2DCURVES_H

    // INCLUDEs
#include <math.h>
#include <stdbool.h> // Boolean function bool
#include "common_functions.h"
// Structures

// Macros

// Function Prototypes
/*
* generateCircleCurve(Point2D* , const Point2D , const double, const double)
* Function to create 2d circle - set of 2 points (coordinate pairs)
* pCurve - pointer to resulting ser of 2d points - 2d curve
* pInitialPoints - initial points - expected 2 points, 1st is center, 2nd is point laying on the circle
* initialPointsCount - count of initial points - expected count is 2 points
* pointsDistance - scalar value representing expceted distance of neighbouring points (real value can differ, because we want to get equidistant discrete curve )
* the function returns boolean value depending on result - true for success / false for failure
*/

    bool generateCircleCurve(Point2D* pCurve, const size_t circlePointsCount, const Point2D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance);
    /*
    * howManyPointsForCircleCurve(const Point2D, const size_t circlePointsCount, const double, const double)
    * Function returns expected number of points needed for discrete circle secified by give parameters
    * circlePointsCount - expected count of discrete circle points
    * pInitialPoints - initial points - expected 2 points, 1st is center, 2nd is point laying on the circle
    * initialPointsCount - count of initial points - expected count is 2 points
    * pointsDistance - scalar value representing expceted distance of neighbouring points (real value can differ, because we want to get equidistant discrete curve )
    * the function returns sclar value representing count of points needed for the circle
    */
    size_t howManyPointsForCircleCurve(const Point2D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance);

    /*
    * double getCirclePerimeter(const double radius)
    * The function calculates a perimeter of a circle given by radius.
    * radius - the radius of the circle
    * the function returns scalar value representing the perimeter of the circle
    */
    double getCirclePerimeter(const double radius);

#endif // !GENERATE2DCURVES_H

#ifdef __cplusplus
}
#endif