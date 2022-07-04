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
	void heatExplicitScheme(Image_Data toExplicitImage, const Filter_Parameters explicitParameters);
	/*
	* Function To Perform Heat Gauss-Seidel Method Implicit Scheme
	*/
	void heatImplicitScheme(Image_Data toImplicitImage, const Filter_Parameters implicitParameters);
	//==============================================================================
	//#endif // !HEAT_QUATION_H

#ifdef __cplusplus
}
#endif
