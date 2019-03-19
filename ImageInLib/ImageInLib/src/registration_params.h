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
	} AffineParameter;

	typedef struct
	{
		dataType tolerance;
		dataType step_size;
		dataType rotation_weight;
		dataType scaling_weight;
		dataType translation_weight;
		dataType h;
		size_t max_iterations;
		AffineParameter affineResults;
		dataType imageBackground;
		bool displayRegistrationOutputs;
	} registrationParams;

	typedef enum
	{
		GRADIENT_DESCENT = 1,
		STOCHASTIC_DESCENT
	} optimizationMethod;

#ifdef __cplusplus
}
#endif