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
	} DistanceMapMethod;

	typedef struct
	{
		dataType tau;
		dataType h;
		dataType objectPixel;
		dataType initValue;
		dataType tolerance;
	} Distance_Map_Params;

	void computeDistanceMap(void **, void **, const size_t imageLength, const size_t imageWidth, const size_t imageHeight, Distance_Map_Params distParams, DistanceMapMethod method);

#ifdef __cplusplus
}
#endif