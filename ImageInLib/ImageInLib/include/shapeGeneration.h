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
		CUBOID = 4
	} ShapeType;

	void generateShape(dataType **inputDataPtr, unsigned char *outputDataPtr, Point3D center, Point3D blockCorner, dataType *fillBlockDimension,
		size_t length, size_t width, size_t height, dataType sphereRadius, dataType smallRadius, dataType fillValue, ShapeType method);

#ifdef __cplusplus
}
#endif