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
	bool nonLinearHeatExplicitScheme(Image_Data inputImageData, FilterParameters explicitParameters);
	/*
	* Function To Perform Heat Gauss-Seidel Method Implicit Scheme
	*/
	bool nonLinearHeatImplicitScheme(Image_Data inputImageData, FilterParameters implicitParameters);

	bool geodesicMeanCurvatureTimeStep(Image_Data inputImageData, FilterParameters filterParameters);

	bool meanCurvatureTimeStep(Image_Data inputImageData, FilterParameters filterParameters);

#endif // !HEAT_QUATION_H

#ifdef __cplusplus
}
#endif
