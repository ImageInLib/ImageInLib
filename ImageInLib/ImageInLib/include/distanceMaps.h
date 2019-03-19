#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "..\src\common_functions.h"

	typedef enum
	{
		FAST_SWEEP = 1,
		FAST_MARCH,
		ROUY_TOURIN,
		BRUTE_FORCE
	} distanceMapMethod;

	typedef struct
	{
		dataType tau;
		dataType h;
		dataType objectPixel;
		dataType initValue;
		dataType tolerance;
	}distanceMapParams;

	void computeDistanceMap(void **, void **, const size_t imageLength, const size_t imageWidth, const size_t imageHeight, distanceMapParams distParams, distanceMapMethod method);

#ifdef __cplusplus
}
#endif