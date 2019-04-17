#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "registration_params.h"

	typedef struct
	{
		Registration_Params regParams;
		Optimization_Method gdescentMethod;
		size_t initEstimateShape;
		size_t meanCalcSteps;
	} shape_Analysis_Parameters;

	typedef struct
	{
		dataType ** eigenvectors;
		dataType * eigenvalues;
		size_t princomp;
	} PCA_Params;

#ifdef __cplusplus
}
#endif
