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
	//3D functions
	void heatExplicitScheme(Image_Data toExplicitImage, const FilterParameters explicitParameters);
	/*
	* Function To Perform Heat Gauss-Seidel Method Implicit Scheme
	*/
	void heatImplicitScheme(Image_Data toImplicitImage, const FilterParameters implicitParameters);
	//==============================================================================
	
	//2D Functions
	/*
	* void heat2dExplicitScheme(Image_Data2D imageData, const FilterParameters explicitParameters)
	* imageData:	Image_Data2D structure to handle 2D images, 
	*					height, width image dimensions
	*					imageDataPtr pointer for pixels value
	* explicitParameters: FilterParameters structure for heat explicit scheme. we need:
						timeStepSize
	                    h
						timeStepsNum
	*/
	void heat2dExplicitScheme(Image_Data2D imageData, const FilterParameters explicitParameters);

	/*
	* void heatImplicit2dScheme(Image_Data2D imageData, const FilterParameters explicitParameters)
	* imageData:	Image_Data2D structure to handle 2D images,
	*					height, width image dimensions
	*					imageDataPtr pointer for pixels value
	* implicitParameters: FilterParameters structure for heat implicit scheme. we need:
						timeStepSize
						h
						timeStepsNum
						omega_c
						tolerance

	*/
	void heatImplicit2dScheme(Image_Data2D imageData, const FilterParameters implicitParameters);

	//#endif // !HEAT_QUATION_H

#ifdef __cplusplus
}
#endif
