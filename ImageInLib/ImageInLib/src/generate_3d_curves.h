#pragma once
#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef GENERATE3DCURVES_H
#define GENERATE3DCURVES_H

    // INCLUDEs
#include <math.h>
#include <stdbool.h>
#include "common_functions.h"

    // Function Prototypes

    /// <summary>
    /// The function calculates the surface of a sphere given by radius.
    /// </summary>
    /// <param name="radius">radius of the sphere</param>
    /// <returns>returns scalar value corresponding to the surface of the sphere</returns>
    double getSphereSurfaceArea(const double radius);

    /// <summary>
    /// Function returns expected number of points needed for discrete circle secified by give parameters
    /// </summary>
    /// <param name="pInitialPoints">initial points - expected 2 points, the first and the last line point</param>
    /// <param name="initialPointsCount">count of initial points - expected count is 2 points</param>
    /// <param name="pointsDistance">scalar value representing expceted distance of neighbouring points (real value can differ, because we want to get equidistant discrete curve )</param>
    /// <returns>the function returns sclar value representing count of points needed for the straight line</returns>
    size_t howManyPointsFor3DStraightLineCurve(const Point3D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance);

    /// <summary>
    /// Function returns expected number of points needed for discrete sphere specified by given parameters
    /// </summary>
    /// <param name="pInitialPoints">initial points - expected 2 points, 1st is center, 2nd is point laying on the circle</param>
    /// <param name="initialPointsCount">count of initial points - expected count is 2 points</param>
    /// <param name="pointsDistance">scalar value representing expected distance of neighbouring points (real value can differ, because we want to get equidistant discrete curve )</param>
    /// <returns></returns>
    size_t howManyPointsForSphereCurve(const Point3D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance);

    /// <summary>
    /// Function to create 3d line - set of 2 points (coordinate pairs)
    /// </summary>
    /// <param name="pCurve">pointer for resulting 3d curve, Be carefull, the function allocates memory for curve points
    /// hold in pPoints. It needs to be released out of the function</param>
    /// <param name="pInitialPoints">initial points - expected 2 points, the first and the last line point</param>
    /// <param name="initialPointsCount">count of initial points - expected count is 2 points</param>
    /// <param name="pointsDistance">scalar value representing expceted distance of neighbouring points (real value can differ, because we want to get equidistant discrete curve )</param>
    /// <returns>the function returns boolean value depending on result - true for success / false for failure</returns>
    bool generate3DLineCurve(Curve3D* pCurve, const Point3D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance);

    /// <summary>
    /// The function to create a sphere - set of 2 points (coordinate pairs)
    /// </summary>
    /// <param name="pCurve">pointer to resulting set of 3d points - 3d curve</param>
    /// <param name="pInitialPoints">initial points - expected 2 points, 1st is center, 2nd is point laying on the circle</param>
    /// <param name="initialPointsCount">count of initial points - expected count is 2 points</param>
    /// <param name="pointsDistance">scalar value representing expceted distance of neighbouring points (real value can differ, because we want to get equidistant discrete curve )</param>
    /// <returns>the function returns boolean value depending on result - true for success / false for failure</returns>
    bool generateSphereCurve(Curve3D* pCurve, const Point3D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance);

#endif // !GENERATE3DCURVES_H

#ifdef __cplusplus
}
#endif
