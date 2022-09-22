#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"
	//==============================================================================
	typedef enum
	{
		SALT_AND_PEPPER = 1,
		ADDITIVE_NOISE = 2,
		MULTIPLICATIVE_NOISE = 3,
		STRUCTURAL_NOISE = 4
	} NoiseType;
	//==============================================================================
	typedef struct {
		float salt_pepper_density, multiplicative_variance;
		int additive_value;
		dataType fgMin, bgMax;
	} NoiseParameters;
	//==============================================================================
	void addNoiseToImage(dataType ** array3DPtr, const size_t xDim, const size_t yDim, const size_t zDim, const NoiseParameters parameters, const NoiseType method);
	//==============================================================================
#ifdef __cplusplus
}
#endif