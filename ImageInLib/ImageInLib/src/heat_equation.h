#ifdef __cplusplus
extern "C" {
#endif

#pragma once
	//#ifndef HEAT_EQUATION_H
	//#define HEAT_EQUATION_H
	//==============================================================================
	/*
	* Header file to contain general functions, variables, structs
	* These can be used in different heat filtering algorithms
	* Avoids redefinitions elsewhere
	*/
	//==============================================================================
	// Includes
#include "common_functions.h"
#include "filter_params.h"
//==============================================================================
// Macro's
//==============================================================================
// STRUCTs
//==============================================================================
// ENUMs
//==============================================================================
// Global Variables
//==============================================================================
// Function Prototypes
/*
* Function To Perform Heat Explicit Scheme
*/
	void heatExplicitScheme(ImageData toExplicitImage, const FilterParameters explicitParameters);
	/*
	* Function To Perform Heat Gauss-Seidel Method Implicit Scheme
	*/
	void heatImplicitScheme(ImageData toImplicitImage, const FilterParameters implicitParameters);
	//==============================================================================
	//#endif // !HEAT_QUATION_H

#ifdef __cplusplus
}
#endif
