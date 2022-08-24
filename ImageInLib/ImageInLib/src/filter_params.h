#ifdef __cplusplus
extern "C" {
#endif

#include "common_functions.h"
#pragma once
	typedef struct
	{
		float timeStepSize;
		float h;
		float sigma;
		float edge_detector_coefficient;
		float omega_c;
		float tolerance;
		float eps2;
		size_t p;
		size_t timeStepsNum;
		size_t maxNumberOfSolverIteration;
		size_t maxNumberOftimeSteps;
	} Filter_Parameters;
	// Parameters used in filtering

#ifdef __cplusplus
}
#endif