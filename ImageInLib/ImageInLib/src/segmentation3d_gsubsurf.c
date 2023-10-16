#include <stdio.h> // Standard lib for input and output functions
#include <stdlib.h>
#include <time.h>
#include <math.h> // Maths functions i.e. pow, sin, cos
#include <stdbool.h> // Boolean function bool
#include <string.h>
#include <common_vtk.h>
#include "file.h"
#include "heat_equation.h"
#include "non_linear_heat_equation.h"
#include "segmentation3D_subsurf.h"
#include "segmentation3d_gsubsurf.h"
#include "image_norm.h"
#include "data_initialization.h"
#include "edgedetection.h"
#include "data_storage.h"
#include "generate_3D_shapes.h"
#include "common_functions.h"
#include "Common_Math.h"
#include "setting_boundary_values.h"
#include "ctype.h"
#include "filter_params.h"
#include "vtk_params.h"

bool generalizedSubsurfSegmentation(Image_Data inputImageData, dataType** initialSegment, Segmentation_Parameters segParameters, FilterParameters implicit_lhe_Parameters,
	Point3D * centers, size_t no_of_centers, unsigned char* outputPathPtr) {

	size_t i, j, k;
	size_t height = inputImageData.height, length = inputImageData.length, width = inputImageData.width;
	size_t dim2D = length * width;
	size_t height_ext = height + 2;
	size_t length_ext = length + 2;
	size_t width_ext = width + 2;
	size_t dim2D_ext = length_ext * width_ext;
	
	dataType h = segParameters.h;
	dataType coef_conv = segParameters.coef_conv;
	dataType coef_dif = segParameters.coef_dif;

	dataType difference_btw_current_and_previous_sol;

	dataType** gauss_seidelPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	dataType** prevSol_extPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);

	dataType** imageToBeSegPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** segmFuntionPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** edgeGradientPtr = (dataType**)malloc(sizeof(dataType*) * height);

	dataType** e_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** w_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** n_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** s_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** t_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** b_Ptr = (dataType**)malloc(sizeof(dataType*) * height);

	dataType** VePtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** VwPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** VnPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** VsPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** VtPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** VbPtr = (dataType**)malloc(sizeof(dataType*) * height);

	//checks if the memory was allocated
	if (gauss_seidelPtr == NULL || prevSol_extPtr == NULL || imageToBeSegPtr == NULL || segmFuntionPtr == NULL || edgeGradientPtr == NULL ||
		e_Ptr == NULL || w_Ptr == NULL || n_Ptr == NULL || s_Ptr == NULL || t_Ptr == NULL || b_Ptr == NULL ||
		VePtr == NULL || VwPtr == NULL || VnPtr == NULL || VsPtr == NULL || VtPtr == NULL || VbPtr == NULL)
		return false;
	for (i = 0; i < height; i++)
	{
		imageToBeSegPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segmFuntionPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		edgeGradientPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);

		e_Ptr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		w_Ptr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		n_Ptr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		s_Ptr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		t_Ptr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		b_Ptr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);

		VePtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		VwPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		VnPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		VsPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		VtPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		VbPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);

		//checks if the memory was allocated
		if (imageToBeSegPtr[i] == NULL || segmFuntionPtr[i] == NULL || edgeGradientPtr[i] == NULL || e_Ptr[i] == NULL
			|| w_Ptr[i] == NULL || n_Ptr[i] == NULL || s_Ptr[i] == NULL || t_Ptr[i] == NULL || b_Ptr[i] == NULL
			|| VePtr[i] == NULL || VwPtr[i] == NULL || VnPtr == NULL || VsPtr == NULL || VtPtr == NULL || VbPtr == NULL)
			return false;
	}

	for (i = 0; i < height_ext; i++)
	{
		gauss_seidelPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
		prevSol_extPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
		//checks if the memory was allocated
		if (gauss_seidelPtr[i] == NULL || prevSol_extPtr[i] == NULL)
			return false;
	}

	//Initialization of structures
	Segment_Image_Data imageData;
	imageData.height = height;
	imageData.length = length;
	imageData.width = width;
	imageData.segmentationFuntionPtr = segmFuntionPtr;

	Coefficient_Pointers CoefPtrs;
	CoefPtrs.e_Ptr = e_Ptr;
	CoefPtrs.w_Ptr = w_Ptr;
	CoefPtrs.n_Ptr = n_Ptr;
	CoefPtrs.s_Ptr = s_Ptr;
	CoefPtrs.t_Ptr = t_Ptr;
	CoefPtrs.b_Ptr = b_Ptr;

	Gradient_Pointers VPtrs;
	VPtrs.GePtr = VePtr;
	VPtrs.GwPtr = VwPtr;
	VPtrs.GnPtr = VnPtr;
	VPtrs.GsPtr = VsPtr;
	VPtrs.GtPtr = VtPtr;
	VPtrs.GbPtr = VbPtr;

	//copy initial segmentation function
	copyDataToAnotherArray(initialSegment, segmFuntionPtr, height, length, width);

	copyDataToExtendedArea(segmFuntionPtr, gauss_seidelPtr, height, length, width);
	copyDataToExtendedArea(segmFuntionPtr, prevSol_extPtr, height, length, width);
	setBoundaryToZeroDirichletBC(gauss_seidelPtr, length_ext, width_ext, height_ext);
	setBoundaryToZeroDirichletBC(prevSol_extPtr, length_ext, width_ext, height_ext);

	//compute coefficients from presmoothed image
	generalizedGFunctionForImageToBeSegmented(inputImageData, edgeGradientPtr, VPtrs, segParameters, implicit_lhe_Parameters, coef_conv);

	//Array for name construction
	unsigned char  name[500];
	unsigned char  name_ending[200];
	Storage_Flags flags = { false,false };

	strcpy_s(name, sizeof name, outputPathPtr);
	sprintf_s(name_ending, sizeof(name_ending), "_edgeDetector.raw");
	strcat_s(name, sizeof(name), name_ending);
	store3dDataArrayD(edgeGradientPtr, length, width, height, name, flags);

	Storage_Flags storageFlags = { false, false };
	manageFile(edgeGradientPtr, length, width, height, name, STORE_DATA_RAW, BINARY_DATA, storageFlags);

	//loop for segmentation time steps
	size_t z = 0; // time steps counter
	do
	{
		z = z + 1;

		setBoundaryToZeroDirichletBC(gauss_seidelPtr, length_ext, width_ext, height_ext);
		setBoundaryToZeroDirichletBC(prevSol_extPtr, length_ext, width_ext, height_ext);

		//calcution of coefficients
		generalizedGaussSeidelCoefficients(imageData, edgeGradientPtr, CoefPtrs, VPtrs, segParameters, coef_conv);

		// Call to function that will evolve segmentation function in each discrete time step
		generalizedSubsurfSegmentationTimeStep(prevSol_extPtr, gauss_seidelPtr, imageData, segParameters, CoefPtrs, centers, no_of_centers);

		//Compute the L2 norm of the difference between the current and previous solutions
		difference_btw_current_and_previous_sol = l2normD(prevSol_extPtr, gauss_seidelPtr, length_ext, width_ext, height_ext, h);

		copyDataToAnotherArray(gauss_seidelPtr, prevSol_extPtr, height_ext, length_ext, width_ext);

		//writing density.
		if ((z % segParameters.savingFrequency) == 0)
		{
			strcpy_s(name, sizeof name, outputPathPtr);
			sprintf_s(name_ending, sizeof(name_ending), "_seg_func_%03zd.raw", z);
			strcat_s(name, sizeof(name), name_ending);

			Storage_Flags storageFlags = { false, false };
			manageFile(imageData.segmentationFuntionPtr, length, width, height, name, STORE_DATA_RAW, BINARY_DATA, storageFlags);
			printf("Step is %zd\n", z);
			printf("Error = %lf\n", difference_btw_current_and_previous_sol);
		}

	} while ((z <= segParameters.maxNoOfTimeSteps) && (difference_btw_current_and_previous_sol > segParameters.segTolerance));

	for (i = 0; i < height; i++)
	{
		free(imageToBeSegPtr[i]);
		free(segmFuntionPtr[i]);
		free(edgeGradientPtr[i]);

		free(e_Ptr[i]);
		free(w_Ptr[i]);
		free(n_Ptr[i]);
		free(s_Ptr[i]);
		free(t_Ptr[i]);
		free(b_Ptr[i]);

		free(VePtr[i]);
		free(VwPtr[i]);
		free(VnPtr[i]);
		free(VsPtr[i]);
		free(VtPtr[i]);
		free(VbPtr[i]);
	}
	free(imageToBeSegPtr);
	free(segmFuntionPtr);
	free(edgeGradientPtr);

	free(e_Ptr);
	free(w_Ptr);
	free(n_Ptr);
	free(s_Ptr);
	free(t_Ptr);
	free(b_Ptr);

	free(VePtr);
	free(VwPtr);
	free(VnPtr);
	free(VsPtr);
	free(VtPtr);
	free(VbPtr);

	for (i = 0; i < height_ext; i++)
	{
		free(prevSol_extPtr[i]);
		free(gauss_seidelPtr[i]);
	}
	free(prevSol_extPtr);
	free(gauss_seidelPtr);

	return true;
}

bool generalizedGFunctionForImageToBeSegmented(Image_Data inputImageData, dataType** edgeGradientPtr, Gradient_Pointers VPtrs,
	Segmentation_Parameters segParameters, FilterParameters implicit_lhe_Parameters, dataType coef_conv)
{
	//checks if the memory was allocated
	if (inputImageData.imageDataPtr == NULL || edgeGradientPtr == NULL || VPtrs.GePtr == NULL || VPtrs.GwPtr == NULL
		|| VPtrs.GnPtr == NULL || VPtrs.GsPtr == NULL || VPtrs.GtPtr == NULL || VPtrs.GbPtr == NULL)
		return false;

	size_t i, j, k, x;
	size_t kplus1, kminus1, iminus1, iplus1, jminus1, jplus1;
	size_t height = inputImageData.height;
	size_t length = inputImageData.length;
	size_t width = inputImageData.width;
	size_t dim2D = length * width;
	size_t k_ext, j_ext, i_ext, x_ext;
	size_t height_ext = height + 2;
	size_t length_ext = length + 2;
	size_t width_ext = width + 2;
	dataType h = segParameters.h, quotient = (dataType)(4.0 * h);
	dataType ux, uy, uz; //change in x, y and z respectively
	dataType u, uN, uS, uE, uW, uNW, uNE, uSE, uSW, Tu, TuN, TuS, TuE, TuW, TuNW, TuNE, TuSE, TuSW, //current and surrounding voxel values
		Bu, BuN, BuS, BuE, BuW, BuNW, BuNE, BuSE, BuSW;
	dataType norm_image_smoothed_e, norm_image_smoothed_w, norm_image_smoothed_n, norm_image_smoothed_s, norm_image_smoothed_t, norm_image_smoothed_b;
	dataType norm_image_smoothed_average;

	dataType** gradient_coef_ext = (dataType**)malloc(sizeof(dataType*) * height_ext);
	dataType** extendedCoefPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	for (k = 0; k < height_ext; k++) {
		gradient_coef_ext[k] = (dataType*)malloc(sizeof(dataType) * length_ext * width_ext);
		extendedCoefPtr[k] = (dataType*)malloc(sizeof(dataType) * length_ext * width_ext);
	}
	if (gradient_coef_ext == NULL || extendedCoefPtr == NULL) 
		return false;

	Image_Data presmoothingData;
	presmoothingData.height = height_ext;
	presmoothingData.length = length_ext;
	presmoothingData.width = width_ext;
	presmoothingData.imageDataPtr = extendedCoefPtr;

	//copy data to extended area which will be used for calculation of diffusion coefficients
	copyDataToExtendedArea(inputImageData.imageDataPtr, extendedCoefPtr, height, length, width);

	//perform reflection of the extended area to ensure zero Neumann boundary condition (for LHE)
	reflection3D(extendedCoefPtr, height_ext, length_ext, width_ext);

	//perfom presmoothing
	heatImplicitScheme(presmoothingData, implicit_lhe_Parameters);

	copyDataToReducedArea(inputImageData.imageDataPtr, extendedCoefPtr, height, length, width);

	//calculation of coefficients
	for (k = 0, k_ext = 1; k < height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				// 2D to 1D representation for i, j
				x_ext = x_new(i_ext, j_ext, length_ext);
				x = x_new(i, j, length);
				iminus1 = i_ext - 1;
				iplus1 = i_ext + 1;
				jplus1 = j_ext + 1;
				jminus1 = j_ext - 1;
				kplus1 = k_ext + 1;
				kminus1 = k_ext - 1;

				//values of voxels in the extended data container for presmoothed image
				u = extendedCoefPtr[k_ext][x_ext];
				uN = extendedCoefPtr[k_ext][x_new(i_ext, jminus1, length_ext)];
				uS = extendedCoefPtr[k_ext][x_new(i_ext, jplus1, length_ext)];
				uE = extendedCoefPtr[k_ext][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				uW = extendedCoefPtr[k_ext][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				uNW = extendedCoefPtr[k_ext][x_new(iminus1, jminus1, length_ext)];
				uNE = extendedCoefPtr[k_ext][x_new(iplus1, jminus1, length_ext)];
				uSE = extendedCoefPtr[k_ext][x_new(iplus1, jplus1, length_ext)];
				uSW = extendedCoefPtr[k_ext][x_new(iminus1, jplus1, length_ext)];
				Tu = extendedCoefPtr[kminus1][x_ext];
				TuN = extendedCoefPtr[kminus1][x_new(i_ext, jminus1, length_ext)];
				TuS = extendedCoefPtr[kminus1][x_new(i_ext, jplus1, length_ext)];
				TuE = extendedCoefPtr[kminus1][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				TuW = extendedCoefPtr[kminus1][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				TuNW = extendedCoefPtr[kminus1][x_new(iminus1, jminus1, length_ext)];
				TuNE = extendedCoefPtr[kminus1][x_new(iplus1, jminus1, length_ext)];
				TuSE = extendedCoefPtr[kminus1][x_new(iplus1, jplus1, length_ext)];
				TuSW = extendedCoefPtr[kminus1][x_new(iminus1, jplus1, length_ext)];
				Bu = extendedCoefPtr[kplus1][x_ext];
				BuN = extendedCoefPtr[kplus1][x_new(i_ext, jminus1, length_ext)];
				BuS = extendedCoefPtr[kplus1][x_new(i_ext, jplus1, length_ext)];
				BuE = extendedCoefPtr[kplus1][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				BuW = extendedCoefPtr[kplus1][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				BuNW = extendedCoefPtr[kplus1][x_new(iminus1, jminus1, length_ext)];
				BuNE = extendedCoefPtr[kplus1][x_new(iplus1, jminus1, length_ext)];
				BuSE = extendedCoefPtr[kplus1][x_new(iplus1, jplus1, length_ext)];
				BuSW = extendedCoefPtr[kplus1][x_new(iminus1, jplus1, length_ext)];

				//calculation of coefficients in the presmooted image data

				// Calculation of coefficients in East direction
				ux = (uE - u) / h;
				uy = ((uN + uNE) - (uS + uSE)) / quotient;
				uz = ((Tu + TuE) - (Bu + BuE)) / quotient;
				norm_image_smoothed_e = sqrt((ux * ux) + (uy * uy) + (uz * uz));

				// Calculation of coefficients in West direction
				ux = (uW - u) / h;
				uy = ((uNW + uN) - (uSW + uS)) / quotient;
				uz = ((TuW + Tu) - (BuW + Bu)) / quotient;
				norm_image_smoothed_w = sqrt((ux * ux) + (uy * uy) + (uz * uz));

				// Calculation of coefficients in North direction
				ux = ((uNE + uE) - (uNW + uW)) / quotient;
				uy = (uN - u) / h;
				uz = ((TuN + Tu) - (BuN + Bu)) / quotient;
				norm_image_smoothed_n = sqrt((ux * ux) + (uy * uy) + (uz * uz));

				// Calculation of coefficients in South direction
				ux = ((uE + uSE) - (uW + uSW)) / quotient;
				uy = (uS - u) / h;
				uz = ((TuS + Tu) - (BuS + Bu)) / quotient;
				norm_image_smoothed_s = sqrt((ux * ux) + (uy * uy) + (uz * uz));

				// Calculation of coefficients in Top direction
				ux = ((TuE + uE) - (TuW + uW)) / quotient;
				uy = ((TuN + uN) - (TuS + uS)) / quotient;
				uz = (Tu - u) / h;
				norm_image_smoothed_t = sqrt((ux * ux) + (uy * uy) + (uz * uz));

				// Calculation of coefficients in Bottom direction
				ux = ((BuW + uW) - (BuE + uE)) / quotient;
				uy = ((BuN + uN) - (BuS + uS)) / quotient;
				uz = (Bu - u) / h;
				norm_image_smoothed_b = sqrt((ux * ux) + (uy * uy) + (uz * uz));

				norm_image_smoothed_average = (dataType)((norm_image_smoothed_e + norm_image_smoothed_w + norm_image_smoothed_n + norm_image_smoothed_s +
					norm_image_smoothed_t + norm_image_smoothed_b) / 6.0);

				edgeGradientPtr[k][x] = edgeDetector(norm_image_smoothed_average * norm_image_smoothed_average, segParameters.coef);

			}
		}
	}

	copyDataToExtendedArea(edgeGradientPtr, gradient_coef_ext, height, length, width);
	reflection3D(gradient_coef_ext, height_ext, length_ext, width_ext);

	for (k = 0, k_ext = 1; k < height; k++, k_ext++) {
		for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
			for (j = 0, j_ext = 1; j < width; j++, j_ext++) {

				x = x_new(i, j, length);
				x_ext = x_new(i_ext, j_ext, length_ext);
				iminus1 = i_ext - 1;
				iplus1 = i_ext + 1;
				jplus1 = j_ext + 1;
				jminus1 = j_ext - 1;
				kplus1 = k_ext + 1;
				kminus1 = k_ext - 1;

				VPtrs.GePtr[k][x] = -coef_conv * (gradient_coef_ext[k_ext][x_new(iplus1, j_ext, length_ext)] - gradient_coef_ext[k_ext][x_new(iminus1, j_ext, length_ext)]) / (2 * h);
				VPtrs.GwPtr[k][x] = -coef_conv * (gradient_coef_ext[k_ext][x_new(iminus1, j_ext, length_ext)] - gradient_coef_ext[k_ext][x_new(iplus1, j_ext, length_ext)]) / (2 * h);
				VPtrs.GnPtr[k][x] = -coef_conv * (gradient_coef_ext[k_ext][x_new(i_ext, jminus1, length_ext)] - gradient_coef_ext[k_ext][x_new(i_ext, jplus1, length_ext)]) / (2 * h);
				VPtrs.GsPtr[k][x] = -coef_conv * (gradient_coef_ext[k_ext][x_new(i_ext, jplus1, length_ext)] - gradient_coef_ext[k_ext][x_new(i_ext, jminus1, length_ext)]) / (2 * h);
				VPtrs.GtPtr[k][x] = -coef_conv * (gradient_coef_ext[kminus1][x_ext] - gradient_coef_ext[kplus1][x_ext]) / (2 * h);
				VPtrs.GbPtr[k][x] = -coef_conv * (gradient_coef_ext[kplus1][x_ext] - gradient_coef_ext[kminus1][x_ext]) / (2 * h);

			}
		}
	}

	for (k = 0; k < height_ext; k++) {
		free(gradient_coef_ext[k]);
		free(extendedCoefPtr[k]);
	}
	free(gradient_coef_ext);
	free(extendedCoefPtr);

	return true;
}

bool generalizedGaussSeidelCoefficients(Segment_Image_Data inputImageData, dataType** edgeGradientPtr, Coefficient_Pointers CoefPtrs, Gradient_Pointers VPtrs, Segmentation_Parameters segParameters, dataType coef_dif)
{
	//checks if the memory was allocated
	if (inputImageData.segmentationFuntionPtr == NULL || edgeGradientPtr == NULL
		|| CoefPtrs.w_Ptr == NULL || CoefPtrs.n_Ptr == NULL || CoefPtrs.s_Ptr == NULL || CoefPtrs.t_Ptr == NULL || CoefPtrs.b_Ptr == NULL
		|| VPtrs.GePtr == NULL || VPtrs.GwPtr == NULL || VPtrs.GnPtr == NULL || VPtrs.GsPtr == NULL || VPtrs.GtPtr == NULL || VPtrs.GbPtr == NULL)
		return false;

	size_t i, j, k, x, x_ext; // length == xDim, width == yDim, height == zDim
	size_t kplus1, kminus1, iminus1, iplus1, jminus1, jplus1;
	size_t dim2D = inputImageData.length * inputImageData.width;
	size_t k_ext, j_ext, i_ext;
	size_t height = inputImageData.height, length = inputImageData.length, width = inputImageData.width;
	size_t height_ext = height + 2;
	size_t length_ext = length + 2;
	size_t width_ext = width + 2;
	dataType h = segParameters.h, quotient = (dataType)(4.0 * segParameters.h);
	dataType orig_ux, orig_uy, orig_uz; //change in x, y and z respectively
	dataType orig_u, orig_uN, orig_uS, orig_uE, orig_uW, orig_uNW, orig_uNE, orig_uSE, orig_uSW, orig_Tu, orig_TuN, orig_TuS,
		orig_TuE, orig_TuW, orig_TuNW, orig_TuNE, orig_TuSE, orig_TuSW, //current and surrounding voxel values
		orig_Bu, orig_BuN, orig_BuS, orig_BuE, orig_BuW, orig_BuNW, orig_BuNE, orig_BuSE, orig_BuSW;
	dataType orig_e, orig_w, orig_n, orig_s, orig_t, orig_b;
	dataType voxel_coef, average_face_coef;

	dataType** extendedCoefPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	for (k = 0; k < height_ext; k++) {
		extendedCoefPtr[k] = (dataType*)malloc(sizeof(dataType) * length_ext * width_ext);
	}
	if (extendedCoefPtr == NULL)
		return false;

	//copy data to extended area which will be used in each time step
	copyDataToExtendedArea(inputImageData.segmentationFuntionPtr, extendedCoefPtr, height, length, width);
	setBoundaryToZeroDirichletBC(extendedCoefPtr, length_ext, width_ext, height_ext);

	//calculation of coefficients
	for (k = 0, k_ext = 1; k < height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				// 2D to 1D representation for i, j
				x_ext = x_new(i_ext, j_ext, length_ext);
				x = x_new(i, j, length);
				iminus1 = i_ext - 1;
				iplus1 = i_ext + 1;
				jplus1 = j_ext + 1;
				jminus1 = j_ext - 1;
				kplus1 = k_ext + 1;
				kminus1 = k_ext - 1;

				//values of voxels in the extended data container for the original image
				orig_u = extendedCoefPtr[k_ext][x_ext];
				orig_uN = extendedCoefPtr[k_ext][x_new(i_ext, jminus1, length_ext)];
				orig_uS = extendedCoefPtr[k_ext][x_new(i_ext, jplus1, length_ext)];
				orig_uE = extendedCoefPtr[k_ext][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				orig_uW = extendedCoefPtr[k_ext][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				orig_uNW = extendedCoefPtr[k_ext][x_new(iminus1, jminus1, length_ext)];
				orig_uNE = extendedCoefPtr[k_ext][x_new(iplus1, jminus1, length_ext)];
				orig_uSE = extendedCoefPtr[k_ext][x_new(iplus1, jplus1, length_ext)];
				orig_uSW = extendedCoefPtr[k_ext][x_new(iminus1, jplus1, length_ext)];
				orig_Tu = extendedCoefPtr[kminus1][x_ext];
				orig_TuN = extendedCoefPtr[kminus1][x_new(i_ext, jminus1, length_ext)];
				orig_TuS = extendedCoefPtr[kminus1][x_new(i_ext, jplus1, length_ext)];
				orig_TuE = extendedCoefPtr[kminus1][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				orig_TuW = extendedCoefPtr[kminus1][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				orig_TuNW = extendedCoefPtr[kminus1][x_new(iminus1, jminus1, length_ext)];
				orig_TuNE = extendedCoefPtr[kminus1][x_new(iplus1, jminus1, length_ext)];
				orig_TuSE = extendedCoefPtr[kminus1][x_new(iplus1, jplus1, length_ext)];
				orig_TuSW = extendedCoefPtr[kminus1][x_new(iminus1, jplus1, length_ext)];
				orig_Bu = extendedCoefPtr[kplus1][x_ext];
				orig_BuN = extendedCoefPtr[kplus1][x_new(i_ext, jminus1, length_ext)];
				orig_BuS = extendedCoefPtr[kplus1][x_new(i_ext, jplus1, length_ext)];
				orig_BuE = extendedCoefPtr[kplus1][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				orig_BuW = extendedCoefPtr[kplus1][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				orig_BuNW = extendedCoefPtr[kplus1][x_new(iminus1, jminus1, length_ext)];
				orig_BuNE = extendedCoefPtr[kplus1][x_new(iplus1, jminus1, length_ext)];
				orig_BuSE = extendedCoefPtr[kplus1][x_new(iplus1, jplus1, length_ext)];
				orig_BuSW = extendedCoefPtr[kplus1][x_new(iminus1, jplus1, length_ext)];

				//calculation of coefficients in the original image data
				// Calculation of coefficients in east direction
				orig_ux = (orig_uE - orig_u) / h;
				orig_uy = ((orig_uN + orig_uNE) - (orig_uS + orig_uSE)) / quotient;
				orig_uz = ((orig_Tu + orig_TuE) - (orig_Bu + orig_BuE)) / quotient;
				orig_e = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in west direction
				orig_ux = (orig_uW - orig_u) / h;
				orig_uy = ((orig_uNW + orig_uN) - (orig_uSW + orig_uS)) / quotient;
				orig_uz = ((orig_TuW + orig_Tu) - (orig_BuW + orig_Bu)) / quotient;
				orig_w = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in north direction
				orig_ux = ((orig_uNE + orig_uE) - (orig_uNW + orig_uW)) / quotient;
				orig_uy = (orig_uN - orig_u) / h;
				orig_uz = ((orig_TuN + orig_Tu) - (orig_BuN + orig_Bu)) / quotient;
				orig_n = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in south direction
				orig_ux = ((orig_uE + orig_uSE) - (orig_uW + orig_uSW)) / quotient;
				orig_uy = (orig_uS - orig_u) / h;
				orig_uz = ((orig_TuS + orig_Tu) - (orig_BuS + orig_Bu)) / quotient;
				orig_s = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in top direction
				orig_ux = ((orig_TuE + orig_uE) - (orig_TuW + orig_uW)) / quotient;
				orig_uy = ((orig_TuN + orig_uN) - (orig_TuS + orig_uS)) / quotient;
				orig_uz = (orig_Tu - orig_u) / h;
				orig_t = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in bottom direction
				orig_ux = ((orig_BuW + orig_uW) - (orig_BuE + orig_uE)) / quotient;
				orig_uy = ((orig_BuN + orig_uN) - (orig_BuS + orig_uS)) / quotient;
				orig_uz = (orig_Bu - orig_u) / h;
				orig_b = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// evaluation of norm of gradient of image at each voxel
				average_face_coef = (dataType)(((orig_e + orig_w + orig_n + orig_s + orig_t + orig_b) / 6.0));

				voxel_coef = (dataType)sqrt(pow(average_face_coef, 2) + segParameters.eps2);

				/* evaluation of norm of gradient of image at each voxel, norm of gradient of presmoothed
				image at each voxel face and reciprocal of norm of gradient of image at each voxel face*/
				CoefPtrs.e_Ptr[k][x] = (dataType)(-min(VPtrs.GePtr[k][x],0) + coef_dif * voxel_coef * edgeGradientPtr[k][x] * (1.0 / orig_e));
				CoefPtrs.w_Ptr[k][x] = (dataType)(-min(VPtrs.GwPtr[k][x],0) + coef_dif * voxel_coef * edgeGradientPtr[k][x] * (1.0 / orig_w));
				CoefPtrs.n_Ptr[k][x] = (dataType)(-min(VPtrs.GnPtr[k][x],0) + coef_dif * voxel_coef * edgeGradientPtr[k][x] * (1.0 / orig_n));
				CoefPtrs.s_Ptr[k][x] = (dataType)(-min(VPtrs.GsPtr[k][x],0) + coef_dif * voxel_coef * edgeGradientPtr[k][x] * (1.0 / orig_s));
				CoefPtrs.t_Ptr[k][x] = (dataType)(-min(VPtrs.GtPtr[k][x],0) + coef_dif * voxel_coef * edgeGradientPtr[k][x] * (1.0 / orig_t));
				CoefPtrs.b_Ptr[k][x] = (dataType)(-min(VPtrs.GbPtr[k][x],0) + coef_dif * voxel_coef * edgeGradientPtr[k][x] * (1.0 / orig_b));
			}
		}
	}

	for (i = 0; i < height_ext; i++) {
		free(extendedCoefPtr[i]);
	}
	free(extendedCoefPtr);

	return true;
}

bool generalizedSubsurfSegmentationTimeStep(dataType** prevSol_extPtr, dataType** gauss_seidelPtr, Segment_Image_Data inputImageData,
	Segmentation_Parameters segParameters, Coefficient_Pointers CoefPtrs, Point3D* centers, size_t no_of_centers)
{
	//check if the memory was allocated successfully
	if (inputImageData.segmentationFuntionPtr == NULL || prevSol_extPtr == NULL || gauss_seidelPtr == NULL || CoefPtrs.e_Ptr == NULL || CoefPtrs.w_Ptr == NULL
		|| CoefPtrs.n_Ptr == NULL || CoefPtrs.s_Ptr == NULL || CoefPtrs.t_Ptr == NULL || CoefPtrs.b_Ptr == NULL)
		return false;

	size_t k, i, j;
	dataType hh = segParameters.h * segParameters.h;
	dataType tau = segParameters.tau;

	// Error value used to check iteration
	// sor - successive over relation value, used in Gauss-Seidel formula
	dataType mean_square_residue = 0.0, gauss_seidel = 0.0;

	// Prepare variables inputImageData.height, inputImageData.length, inputImageData.width
	size_t height = inputImageData.height, length = inputImageData.length, width = inputImageData.width;
	size_t height_ext = height + 2;
	size_t length_ext = length + 2;
	size_t width_ext = width + 2;
	size_t k_ext, j_ext, i_ext;
	size_t x; //x = x_new(i, j, length);
	size_t x_ext; //x_ext = x_new(i_ext, j_ext, length_ext);
	size_t z; // Steps counter

	const dataType coef_tauh = tau / hh;

	// The Implicit Scheme Evaluation
	z = 0;
	do
	{
		z = z + 1;
		for (k = 0, k_ext = 1; k < height; k++, k_ext++)
		{
			for (i = 0, i_ext = 1; i < length; i++, i_ext++)
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++)
				{
					// 2D to 1D representation for i, j
					x_ext = x_new(i_ext, j_ext, length_ext);
					x = x_new(i, j, length);

					// Gauss-Seidel Formula Evaluation
					gauss_seidel = (dataType)(((prevSol_extPtr[k_ext][x_ext] + coef_tauh * (CoefPtrs.e_Ptr[k][x] * gauss_seidelPtr[k_ext][x_ext + 1] 
					+ CoefPtrs.w_Ptr[k][x] * gauss_seidelPtr[k_ext][x_ext - 1] + CoefPtrs.s_Ptr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)]
					+ CoefPtrs.n_Ptr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)]
					+ CoefPtrs.b_Ptr[k][x] * gauss_seidelPtr[k_ext + 1][x_ext] + CoefPtrs.t_Ptr[k][x] * gauss_seidelPtr[k_ext - 1][x_ext]))
					/ (1 + coef_tauh * (CoefPtrs.e_Ptr[k][x] + CoefPtrs.w_Ptr[k][x] + CoefPtrs.s_Ptr[k][x] + CoefPtrs.n_Ptr[k][x] + CoefPtrs.b_Ptr[k][x] + CoefPtrs.t_Ptr[k][x]))));

					// SOR implementation using Gauss-Seidel
					if (segParameters.initialSegmentAsDirichletBoundaryCondition == true) {
						if (gauss_seidelPtr[k_ext][x_ext] != 1.0) {
							gauss_seidelPtr[k_ext][x_ext] = gauss_seidelPtr[k_ext][x_ext] + segParameters.omega_c * (gauss_seidel - gauss_seidelPtr[k_ext][x_ext]);
						}
					}
					else {
						gauss_seidelPtr[k_ext][x_ext] = gauss_seidelPtr[k_ext][x_ext] + segParameters.omega_c * (gauss_seidel - gauss_seidelPtr[k_ext][x_ext]);
					}
				}
			}
		}

		// Error Evaluation
		mean_square_residue = 0.0; // Initialize
		for (k = 0, k_ext = 1; k < height; k++, k_ext++)
		{
			for (i = 0, i_ext = 1; i < length; i++, i_ext++)
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++)
				{
					// 2D to 1D representation for i, j
					x_ext = x_new(i_ext, j_ext, length_ext);
					x = x_new(i, j, length);

					mean_square_residue += (dataType)(pow(gauss_seidelPtr[k_ext][x_ext] * (1 + coef_tauh * (CoefPtrs.e_Ptr[k][x] + CoefPtrs.w_Ptr[k][x] + CoefPtrs.s_Ptr[k][x] + CoefPtrs.n_Ptr[k][x] + CoefPtrs.b_Ptr[k][x] + CoefPtrs.t_Ptr[k][x]))
						- coef_tauh * (CoefPtrs.e_Ptr[k][x] * gauss_seidelPtr[k_ext][x_ext + 1] + CoefPtrs.w_Ptr[k][x] * gauss_seidelPtr[k_ext][x_ext - 1]
							+ CoefPtrs.s_Ptr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)]
							+ CoefPtrs.n_Ptr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)]
							+ CoefPtrs.b_Ptr[k][x] * gauss_seidelPtr[k_ext + 1][x_ext] + CoefPtrs.t_Ptr[k][x] * gauss_seidelPtr[k_ext - 1][x_ext])
						- prevSol_extPtr[k_ext][x_ext], 2) * hh);
				}
			}
		}

	} while (mean_square_residue > segParameters.gauss_seidelTolerance && z < segParameters.maxNoGSIteration);

	if (no_of_centers == 1)
	{
		rescaleToIntervalZeroOne(gauss_seidelPtr, length_ext, width_ext, height_ext);
	}
	else
	{
		for (i = 0; i < no_of_centers; i++)
		{
			rescaleLocallyToIntervalZeroOne(gauss_seidelPtr, length_ext, width_ext, height_ext, centers[i].x, centers[i].y, centers[i].z, 6., i);
		}
	}

	//Copy the current time step to original data array after timeStepsNum
	copyDataToReducedArea(inputImageData.segmentationFuntionPtr, gauss_seidelPtr, height, length, width);

	return true;
}
