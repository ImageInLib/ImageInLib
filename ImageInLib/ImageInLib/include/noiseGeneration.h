#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"

	typedef enum
	{
		SALT_AND_PEPPER = 1,
		ADDITIVE_NOISE = 2,
		MULTIPLICATIVE_NOISE = 3
	} NoiseType;

	void addNoiseToImage(dataType ** array3DPtr, const size_t xDim, const size_t yDim, const size_t zDim,
		float density, const NoiseType method, dataType pepper); //adding of pepper

#ifdef __cplusplus
}
#endif