#ifdef __cplusplus
extern "C" {
#endif

	/*
	* Author: Markjoe Olunna UBA
	* Purpose: ImageInLife project - 4D Image Segmentation Methods
	* Language:  C
	*/
#pragma once
#ifndef HEAT_EQUATION_H
#define HEAT_EQUATION_H

	// INCLUDEs
#include <stdbool.h>
#include <stddef.h>
#include "common_functions.h" // ImageData structs, Reflection functions
#include "filter_params.h"
// MACROs

// STRUCTs
// FUNCTION PROTOTYPES
/*explicitParameters
* Function To Perform Heat Explicit Scheme
*/
	bool nonLinearHeatExplicitScheme(Image_Data inputImageData, Filter_Parameters explicitParameters);
	/*
	* Function To Perform Heat Gauss-Seidel Method Implicit Scheme
	*/
	bool nonLinearHeatImplicitScheme(Image_Data inputImageData, Filter_Parameters implicitParameters, size_t numberOfTimeStep);

	bool geodesicMeanCurvatureTimeStep(Image_Data inputImageData, Filter_Parameters filterParameters,
		const size_t maxNumberOfSolverIteration, dataType  coef, dataType  eps2, size_t numberOfTimeStep);

	bool meanCurvatureTimeStep(Image_Data inputImageData, Filter_Parameters filterParameters,
		const size_t maxNumberOfSolverIteration, dataType  eps2, size_t numberOfTimeStep);

#endif // !HEAT_QUATION_H

#ifdef __cplusplus
}
#endif
