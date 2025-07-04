#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"

	typedef enum
	{
		SPHERE = 1,
		SOLID_SPHERE = 2,
		SPHERE_WITH_HOLES = 3,
		CUBOID = 4,
		CIRCLE_2D_CURVE = 5,
		LINE_2D_CURVE = 6,
		SPHERE_3D_CURVE = 7,
		LINE_3D_CURVE = 8
	} ShapeType;

	void generateShape(dataType **inputDataPtr, unsigned char *outputDataPtr, Point3D center, Point3D blockCorner, dataType *fillBlockDimension,
		size_t length, size_t width, size_t height, dataType sphereRadius, dataType smallRadius, dataType fillValue, ShapeType method);

	/// <summary>
	/// The function generates given 2D curve based on initial points and chosen shape
	/// </summary>
	/// <param name="pcurve">Destination curve structure object. The curve will be filled by generated one.</param>
	/// <param name="pInitial_points">Initial points needed for a curve generation (straight line needs 2 points - two ending points, 
	/// circle needs 2 points - center and curve point)</param>
	/// <param name="initial_points_count">Count of the initial points (straight line needs 2 points,, 
	/// circle needs 2 points)</param>
	/// <param name="points_distance">Wanted distance between neighboring points - defines discrete curve points density</param>
	/// <param name="shape">Wanted shape to be generated</param>
	/// <returns>Returns true, if the function succeeded, false otherwise</returns>
	bool generate2DCurve(Curve2D* pcurve, const Point2D* pinitial_points, const size_t initial_points_count, const double points_distance, ShapeType shape);

	/// <summary>
	/// The function calculating numebr of points needed for a 2D curve generation
	/// </summary>
	/// <param name="pinitial_points">Initial points neede for a curve generation (straight line needs 2 points - two ending points, 
	/// circle needs 2 points - center and curve point)</param>
	/// <param name="initial_points_count">Count of the initial points (straight line needs 2 points,, 
	/// circle needs 2 points)</param>
	/// <param name="points_distance">Wanted distance between neighboring points - defines discrete curve points density</param>
	/// <param name="method"></param>
	/// <returns>Returns the number of resulting curve points</returns>
	size_t getNumberOfExpected2DCurvePoints(const Point2D* pinitial_points, const size_t initial_points_count, const double points_distance, ShapeType shape);

	/// <summary>
	///  The function generates given 3D curve based on initial points and chosen shape
	/// </summary>
	/// <param name="pcurve">Destination curve structure object. The curve will be filled by generated one.</param>
	/// <param name="pinitial_points">Initial points needed for a curve generation (3D line needs two ending points, 
	/// sphere (envelop) needs 2 points - center and curve point)</param>
	/// <param name="initial_points_count">Count of the initial points (3D line needs 2 points,, 
	/// sphere needs 2 points)</param>
	/// <param name="points_distance">Wanted distance between neighboring points - defines discrete curve points density</param>
	/// <param name="shape">Wanted shape to be generated</param>
	/// <returns>Returns true, if the function succeeded, false otherwise</returns>
	bool generate3DCurve(Curve3D* pCurve, const Point3D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance, ShapeType shape);

	/// <summary>
	/// The function calculating number of points needed for a 3D curve generation
	/// </summary>
	/// <param name="pinitial_points">Initial points needed for a curve generation (3D line needs two ending points, 
	/// sphere needs 2 points - center and curve point)</param>
	/// <param name="initial_points_count">Count of the initial points (3D line needs 2 points, 
	/// sphere needs 2 points)</param>
	/// <param name="points_distance">Wanted distance between neighboring points - defines discrete curve points density</param>
	/// <param name="method"></param>
	/// <returns>Returns the number of resulting curve points</returns>
	size_t getNumberOfExpected3DCurvePoints(const Point3D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance, ShapeType shape);

#ifdef __cplusplus
}
#endif