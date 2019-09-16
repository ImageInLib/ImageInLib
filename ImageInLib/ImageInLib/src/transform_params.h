#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"
#include <stdbool.h>

	typedef struct
	{
		Point3D translation;
		Point3D scaling;
		Point3D rotation;
		dataType centroid[3];
		dataType imageBackground;
		dataType imageForeground;
		bool parallelize;
	} Affine_Transform_Parameters;

#ifdef __cplusplus
}
#endif