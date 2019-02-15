#pragma once
#include "common_functions.h"

typedef struct
{
	Point3D translation;
	Point3D scaling;
	Point3D rotation;
	dataType centroid[3];
	dataType imageBackground;
} affineTransformParameters;