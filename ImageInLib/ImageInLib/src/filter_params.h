#ifdef __cplusplus
extern "C" {
#endif

#include "common_functions.h"
#pragma once
	typedef struct
	{
		double timeStepSize;
		double h;
		double sigma;
		double edge_detector_coefficient;
		double omega_c;
		double tolerance;
		double eps2;
		double coef;
		size_t p;
		size_t timeStepsNum;
		size_t maxNumberOfSolverIteration;
	} Filter_Parameters;
	// Parameters used in filtering

#ifdef __cplusplus
}
#endif