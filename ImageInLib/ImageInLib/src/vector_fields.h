#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef VECTOR_FIELDS
#define VECTOR_FIELDS
	//==============================================================================
	/*
	* Header file to contain general functions, variables, structs
	* These can be used in different vector field
	* Avoids redefinitions elsewhere
	*/
	//==============================================================================
	// INCLUDEs
#include "common_functions.h" // ImageData Structure
//==============================================================================
// MACROs
//==============================================================================
// STRUCTs
// Needed in calculation of vector field points
	typedef struct
	{
		dataType h;
		size_t p;
	} Vector_Parameters; // Needed in calculation of vector field points
	// Structure for the Vector Field Direction
	//==============================================================================
	typedef struct
	{
		// Container to store the vector field calculated values.
		dataType ** fieldPtr;
	} Vector_Direction;
	//==============================================================================
	// Function Prototypes
	/*
	* Function generating 3D vector field
	* Field Direction according to the 6 faces - Top, Bottom, East, West, North, South
	* vectorVariables - Variables used in calculation of field values
	* vectorDataPtr - The calculated vector field values are stored in this struct pointer
	* inputDataPtr Input data for which we calculate it's vector fields
	* fieldDirection String for the finite difference method we want to calculate/evaluate
	* Coeff is the k used in calculation for gradient
	*/
	void generate3DVector(Point3D ** vectorPtr, Image_Data inputDataPtr, Vector_Parameters vectorVariables, enum FiniteDifference fieldDirection, dataType coef);
	//==============================================================================
#endif // !VECTOR_FIELDS

#ifdef __cplusplus
}
#endif
