#ifdef __cplusplus
extern "C" {
#endif

	/*
	* Author: Markjoe Olunna UBA
	* Purpose: ImageInLife project - 4D Image Segmentation Methods
	* Language:  C
	*/
#pragma once
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h> // Standard lib for input and output functions

	// INCLUDEs
#include "common_functions.h" // ImageData structs, Reflection functions
#include "heat_equation.h"
#include "filter_params.h"

// MACROs

// STRUCTs

//Structure that holds the parameters used during SUBSURF segmentation process.
	typedef struct
	{
		size_t maxNoGSIteration;// Maximum number of Gauss-Seidel iterations
		dataType coef; // constant K in the Perona-Malik function G for the image
		dataType eps2; // epsilon is the regularization factor (Evans-Spruck)
		size_t numberOfTimeStep;// Number of current time step
		size_t maxNoOfTimeSteps;// Maximum number of time step
		size_t mod; // Kind of writing density
		dataType segTolerance; // Tolerance for stopping of the segmentation process
		dataType tau, h, omega_c, gauss_seidelTolerance;
		/* h is the Grid size, tau is time step for the segmentation process,
		   omega_c is the relaxation parameter in SOR implementation using Gauss-Seidel, gauss_seidelTolerance is the acceptable
		   tolerance for Gauss-Seidel iterations*/
	} Segmentation_Parameters;

	// Structure that holds the coefficients for PM function G calculated from the image
	typedef struct
	{
		dataType ** GePtr; // East coefficient pointer for PM function G
		dataType ** GwPtr; // West coefficient pointer for PM function G
		dataType ** GnPtr; // North coefficient pointer for PM function G
		dataType ** GsPtr; // South coefficient pointer for PM function G
		dataType ** GtPtr; // Top coefficient pointer for PM function G
		dataType ** GbPtr; // Bottom coefficient pointer for PM function G
	} Gradient_Pointers;

	// Structure that holds the coefficients for calculations
	typedef struct
	{
		dataType ** e_Ptr; // East coefficient pointer
		dataType ** w_Ptr; // West coefficient pointer
		dataType ** n_Ptr; // North coefficient pointer
		dataType ** s_Ptr; // South coefficient pointer
		dataType ** t_Ptr; // Top coefficient pointer
		dataType ** b_Ptr; // Bottom coefficient pointer
	} Coefficient_Pointers;

	//Structure that holds the image and its dimensions
	typedef struct {
		// Image Dimensions
		size_t height, length, width; // Absolute Dimension

		dataType **segmentationFuntionPtr; // Segmentation funtion
		dataType **inputImageToBeSegmented; // input image to be segmented
	} Segment_Image_Data;

	//Structure that holds the extended image and its dimensions
	typedef struct {
		dataType **prevSol_extPtr; // Pointer to the Previous solution
		dataType **gauss_seidelPtr; // Pointer to array used during Gauss-Seidel iterations
		dataType **prevSol_p_extPtr; //
	} Extended_ImageData;

	// FUNCTION PROTOTYPES
	/*generateInitialSegmentationFunctionForMultipleCentres generates initial segmentation function for multiple centers
	-center_s is the center which is given. This center has x, y and z coordinates
	-R usually corresponds to the halved inner diameter of the image domain
	-1/v gives a maximum of the segmentation function
	-counter counts the numbers of centers available*/
	bool generateInitialSegmentationFunctionForMultipleCentres(dataType **inputDataArrayPtr, size_t length, size_t width, size_t height,
		Point3D *centers, dataType v, dataType R, size_t no_of_centers);

	/*subsurfSegmentationTimeStep performs one time step of SUBSURF segmentation
	- prevSol_extPtr is pointer the
	- gauss_seidelPtr
	- inputImageData is structure that holds the image and its dimensions
	- segParameters is structure that holds the parameters used during SUBSURF segmentation process.
	- GPtrs is structure that holds the coefficients for PM function G
	- CoefPtrs is structure that holds the coefficients for calculations
	- no_of_centers no of centers (usually more than one during segmentation of multiple cells)
	- center_x, center_y, center_z are pointers to the center coordinates*/
	bool subsurfSegmentationTimeStep(dataType **prevSol_extPtr, dataType **gauss_seidelPtr, Segment_Image_Data inputImageData, Gradient_Pointers GPtrs,
		Segmentation_Parameters segParameters, Coefficient_Pointers CoefPtrs, Point3D * centers, size_t no_of_centers);

	/*subsurfSegmentation performs the whole segmentation process
	- Allocation and deallocation of memory are done by the function
	- Fuction call to subsurfSegmentationTimeStep is made by this function also
	- Storage of segmentation result to file is also done by this fuction
	- Generation or loading of object to be segmented is done this function
	- Finally, it also generates initial segmentation function
	- no_of_centers is no of centers (usually more than one during segmentation of multiple cells)
	- center_x, center_y, center_z are pointers to the center coordinates*/
	bool subsurfSegmentation(Image_Data inputImageData, Segmentation_Parameters segParameters, Filter_Parameters explicit_lhe_Parameters,
		Point3D * centers, size_t no_of_centers, unsigned char * outputPathPtr);

	/* gFunctionForImageToBeSegmented manages computation of norm of presmoothed image to be segmented
	-It first makes a copy of the image to be segmented, presmooth the copy and evaluate diffusion coefficients on the presmoothed image.
	-It save the result in the GPtrs
	-imageToBeSegPtr is pointer to the image to be segmented
	-extendedCoefPtr is pointer to the copy of the image to be segmented
	-GPtrs is structure that holds the coefficients for PM function G
	-segParameters is structure that holds the parameters used during SUBSURF segmentation process
	-explicit_lhe_Parameters is structure that holds the parameters used for presmoothing or solving LHE*/
	bool gFunctionForImageToBeSegmented(Image_Data inputImageData, dataType **extendedCoefPtr, Gradient_Pointers GPtrs,
		Segmentation_Parameters segParameters, Filter_Parameters explicit_lhe_Parameters);

	/* gaussSeidelCoefficients calculates coefficient used during Gauss-Seidel iterations
	- inputImageData is structure that holds the image and its dimensions
	- segParameters is structure that holds the parameters used during SUBSURF segmentation process.
	- GPtrs is structure that holds the coefficients for PM function G
	- CoefPtrs is structure that holds the coefficients for calculations
	*/
	bool gaussSeidelCoefficients(dataType **extendedCoefPtr, Segment_Image_Data inputImageData, Gradient_Pointers GPtrs,
		Coefficient_Pointers CoefPtrs, Segmentation_Parameters segParameters);

	bool rescaleToIntervalZeroOne(dataType **imagePtr, size_t length, size_t width, size_t height);

	/*rescaleLocallyToIntervalZeroOne rescales each local neighborhood in an image of given dimensions to interval [0,1]
	-R is the radius of the neighborhood
	-counter keeps track of the count
	-center_x, center_y, center_z are the center coordinates*/
	bool rescaleLocallyToIntervalZeroOne(dataType **imagePtr, size_t length, size_t width, size_t height,
		dataType center_x, dataType center_y, dataType center_z, dataType R, size_t counter);

#ifdef __cplusplus
}
#endif