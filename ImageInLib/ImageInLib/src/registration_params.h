#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"
#include <stdbool.h>

	//==============================================================================
	// STRUCTS
	typedef struct {
		Point3D translation;
		Point3D rotation;
		Point3D scaling;
	} Affine_Parameter;

	typedef struct
	{
		dataType tolerance;
		dataType step_size;
		dataType rotation_weight;
		dataType scaling_weight;
		dataType translation_weight;
		dataType h;
		size_t max_iterations;
		size_t rand_points;
		Affine_Parameter affineResults;
		dataType imageBackground;
		dataType imageForeground;
		bool displayRegistrationOutputs;
		char fPathname[100];
	} Registration_Params;

	typedef enum
	{
		GRADIENT_DESCENT = 1,
		STOCHASTIC_DESCENT
	} Optimization_Method;

#ifdef __cplusplus
}
#endif