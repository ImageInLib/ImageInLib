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
		CIRCLE_CURVE = 5
	} ShapeType;

	void generateShape(dataType **inputDataPtr, unsigned char *outputDataPtr, Point3D center, Point3D blockCorner, dataType *fillBlockDimension,
		size_t length, size_t width, size_t height, dataType sphereRadius, dataType smallRadius, dataType fillValue, ShapeType method);

	bool generate2DCurve(Point2D* pCurve, const size_t curvePointsCount, const Point2D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance, ShapeType method);

	size_t getNumberOfExpected2DCurvePoints(const Point2D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance, ShapeType method);

#ifdef __cplusplus
}
#endif