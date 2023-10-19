#include "common_functions.h"
#include "trajectories.h"
#include "shapeGeneration.h"
#include "generate_3d_shapes.h"
#include "generate_2d_curves.h"

void generateShape(dataType **inputDataPtr, unsigned char *outputDataPtr, Point3D center, Point3D blockCorner, dataType *fillBlockDimension,
	size_t length, size_t width, size_t height, dataType sphereRadius, dataType smallRadius, dataType fillValue, ShapeType method)
{
	switch (method)
	{
	case SPHERE:
		generateSphere(inputDataPtr, center, length, width, height, sphereRadius, fillValue, outputDataPtr);
		break;
	case SOLID_SPHERE:
		fillBall3D(inputDataPtr, height, length, width, sphereRadius, center, fillValue);
		break;
	case SPHERE_WITH_HOLES:
		generateSphereWithSixHoles(inputDataPtr, center, length, width, height,
			sphereRadius, smallRadius, fillValue, outputDataPtr);
		break;
	case CUBOID:
		fillBlock3D(inputDataPtr, height, length, width, blockCorner, fillBlockDimension, fillValue);
		break;
	default:
		break;
	}
}

bool generate2DCurve(Curve2D* pCurve, const Point2D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance, ShapeType method)
{
	switch (method)
	{
	case CIRCLE_2D_CURVE:
		return generateCircleCurve(pCurve, pInitialPoints, initialPointsCount, pointsDistance);
	case LINE_2D_CURVE:
		return generateStraightLineCurve(pCurve, pInitialPoints, initialPointsCount, pointsDistance);
	default:
		return false;
	}
}

size_t getNumberOfExpected2DCurvePoints(const Point2D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance, ShapeType method)
{
	switch (method)
	{
	case CIRCLE_2D_CURVE:
		return howManyPointsForCircleCurve(pInitialPoints, initialPointsCount, pointsDistance);
	case LINE_2D_CURVE:
		return howManyPointsForStraightLineCurve(pInitialPoints, initialPointsCount, pointsDistance);
	default:
		return 0;
	}
}