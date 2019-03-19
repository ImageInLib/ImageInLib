#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "registration_params.h"

	typedef struct
	{
		registrationParams regParams;
		optimizationMethod gdescentMethod;
		size_t initEstimateShape;
		size_t meanCalcSteps;
	} shapeAnalysisParameters;

	typedef struct
	{
		dataType ** eigenvectors;
		dataType * eigenvalues;
		size_t princomp;
	} pcaParams;

#ifdef __cplusplus
}
#endif
