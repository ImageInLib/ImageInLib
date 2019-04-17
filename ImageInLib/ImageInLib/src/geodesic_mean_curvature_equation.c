/*
* Author: Markjoe Olunna UBA
* Purpose: ImageInLife project - 4D Image Segmentation Methods
* Language:  C
*/
#include <stdio.h> // Standard lib for input and output functions
#include <stdlib.h>
#include <math.h> // Maths functions i.e. pow, sin, cos
#include <stdbool.h> // Boolean function bool
#include "heat_equation.h"
#include "non_linear_heat_equation.h"
#include "setting_boundary_values.h"
#include "common_functions.h"
#include "filter_params.h"

// Local Function Prototype

bool geodesicMeanCurvatureTimeStep(Image_Data inputImageData, Filter_Parameters filterParameters,
	const size_t maxNumberOfSolverIteration, dataType coef, dataType eps2, size_t numberOfTimeStep)
{
	//checks if the memory was allocated
	if (inputImageData.imageDataPtr == NULL)
		return false;

	size_t k, i, j;
	dataType hh = filterParameters.h * filterParameters.h;
	dataType tau = 600 * filterParameters.timeStepSize;
	// Error value used to check iteration
	// sor - successive over relation value, used in Gauss-Seidel formula
	dataType error, gauss_seidel;

	// Prepare variables toExplicitImage.height, toExplicitImage.length, toExplicitImage.width
	size_t height = inputImageData.height, length = inputImageData.length, width = inputImageData.width;
	size_t height_ext = height + 2;
	size_t length_ext = length + 2;
	size_t width_ext = width + 2;
	size_t k_ext, j_ext, i_ext;
	dataType ux, uy, uz, orig_ux, orig_uy, orig_uz; //change in x, y and z respectively
	size_t x; //x = x_new(i, j, length);
	size_t x_ext; //x_ext = x_new(i_ext, j_ext, length_ext);
	size_t z; // Steps counter

	const dataType coef_tauh = tau / hh;
	dataType u, uN, uS, uE, uW, uNW, uNE, uSE, uSW, Tu, TuN, TuS, TuE, TuW, TuNW, TuNE, TuSE, TuSW, //current and surrounding voxel values
		Bu, BuN, BuS, BuE, BuW, BuNW, BuNE, BuSE, BuSW;
	dataType orig_u, orig_uN, orig_uS, orig_uE, orig_uW, orig_uNW, orig_uNE, orig_uSE, orig_uSW, orig_Tu, orig_TuN, orig_TuS,
		orig_TuE, orig_TuW, orig_TuNW, orig_TuNE, orig_TuSE, orig_TuSW, //current and surrounding voxel values
		orig_Bu, orig_BuN, orig_BuS, orig_BuE, orig_BuW, orig_BuNW, orig_BuNE, orig_BuSE, orig_BuSW;
	dataType voxel_coef, average_face_coef;
	size_t kplus1, kminus1, iminus1, iplus1, jminus1, jplus1;

	Image_Data presmoothingData;
	presmoothingData.height = height_ext;
	presmoothingData.length = length_ext;
	presmoothingData.width = width_ext;

	// Create temporary Image Data holder for Previous time step data - with extended boundary because of boundary condition
	dataType ** prevSolPtr = (dataType **)malloc(sizeof(dataType *) * (height_ext));

	// Create temporary Image Data holder for Current time step data - with extended boundary because of boundary condition
	dataType ** gauss_seidelPtr = (dataType **)malloc(sizeof(dataType *) * (height_ext));

	/* Create tempporary Image Data holder for calculation of diffusion coefficients on presmoothed image
	- with extended boundary because of boundary condition*/
	dataType ** presmoothed_coefPtr = (dataType **)malloc(sizeof(dataType *) * (height_ext));

	//checks if the memory was allocated
	if (prevSolPtr == NULL || gauss_seidelPtr == NULL || presmoothed_coefPtr == NULL)
		return false;

	for (k = 0; k < height_ext; k++)
	{
		presmoothed_coefPtr[k] = malloc(sizeof(dataType)*(length_ext)*(width_ext));
		gauss_seidelPtr[k] = malloc(sizeof(dataType)*(length_ext)*(width_ext));
		prevSolPtr[k] = malloc(sizeof(dataType)*(length_ext)*(width_ext));
		//checks if the memory was allocated
		if (presmoothed_coefPtr[k] == NULL || gauss_seidelPtr[k] == NULL || prevSolPtr[k] == NULL)
			return false;
	}

	/* Create tempporary Image Data holder for diffusion coefficients
	- with extended boundary because of boundary condition*/
	dataType ** presmoot_e_coefPtr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** presmoot_w_coefPtr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** presmoot_n_coefPtr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** presmoot_s_coefPtr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** presmoot_t_coefPtr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** presmoot_b_coefPtr = (dataType **)malloc(sizeof(dataType *) * height);
	//checks if the memory was allocated
	if (presmoot_e_coefPtr == NULL || presmoot_w_coefPtr == NULL || presmoot_n_coefPtr == NULL ||
		presmoot_s_coefPtr == NULL || presmoot_t_coefPtr == NULL || presmoot_b_coefPtr == NULL)
		return false;

	for (k = 0; k < height; k++)
	{
		presmoot_e_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		presmoot_w_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		presmoot_n_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		presmoot_s_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		presmoot_t_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		presmoot_b_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		//checks if the memory was allocated
		if (presmoot_e_coefPtr[k] == NULL || presmoot_w_coefPtr[k] == NULL || presmoot_n_coefPtr[k] == NULL ||
			presmoot_s_coefPtr[k] == NULL || presmoot_t_coefPtr[k] == NULL || presmoot_b_coefPtr[k] == NULL)
			return false;
	}

	dataType ** orig_e_coefPtr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** orig_w_coefPtr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** orig_n_coefPtr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** orig_s_coefPtr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** orig_t_coefPtr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** orig_b_coefPtr = (dataType **)malloc(sizeof(dataType *) * height);
	//checks if the memory was allocated
	if (orig_e_coefPtr == NULL || orig_w_coefPtr == NULL || orig_n_coefPtr == NULL || orig_s_coefPtr == NULL ||
		orig_t_coefPtr == NULL || orig_b_coefPtr == NULL)
		return false;

	for (k = 0; k < height; k++)
	{
		orig_e_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		orig_w_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		orig_n_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		orig_s_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		orig_t_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		orig_b_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		//checks if the memory was allocated
		if (orig_e_coefPtr[k] == NULL || orig_w_coefPtr[k] == NULL || orig_n_coefPtr[k] == NULL || orig_s_coefPtr[k] == NULL
			|| orig_t_coefPtr[k] == NULL || orig_b_coefPtr[k] == NULL)
			return false;
	}

	dataType ** coefPtr_e = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** coefPtr_w = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** coefPtr_n = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** coefPtr_s = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** coefPtr_t = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** coefPtr_b = (dataType **)malloc(sizeof(dataType *) * height);
	//checks if the memory was allocated
	if (coefPtr_e == NULL || coefPtr_w == NULL || coefPtr_n == NULL || coefPtr_s == NULL || coefPtr_t == NULL ||
		coefPtr_b == NULL)
		return false;

	for (k = 0; k < height; k++)
	{
		coefPtr_e[k] = malloc(sizeof(dataType) * length * width);
		coefPtr_w[k] = malloc(sizeof(dataType) * length * width);
		coefPtr_n[k] = malloc(sizeof(dataType) * length * width);
		coefPtr_s[k] = malloc(sizeof(dataType) * length * width);
		coefPtr_t[k] = malloc(sizeof(dataType) * length * width);
		coefPtr_b[k] = malloc(sizeof(dataType) * length * width);
		//checks if the memory was allocated
		if (coefPtr_e[k] == NULL || coefPtr_w[k] == NULL || coefPtr_n[k] == NULL || coefPtr_s[k] == NULL ||
			coefPtr_t[k] == NULL || coefPtr_b[k] == NULL)
			return false;
	}

	presmoothingData.imageDataPtr = presmoothed_coefPtr;

	//copy data to extended area which will be used in each time step
	copyDataToExtendedArea(inputImageData.imageDataPtr, prevSolPtr, height, length, width);
	//copy data to extended area which will be used in each Gauss Seidel iteration
	copyDataToExtendedArea(inputImageData.imageDataPtr, gauss_seidelPtr, height, length, width);
	//copy data to extended area which will be used for calculation of diffusion coefficients
	copyDataToExtendedArea(inputImageData.imageDataPtr, presmoothed_coefPtr, height, length, width);

	//perform reflection of the extended area to ensure zero Neumann boundary condition (for LHE)
	reflection3D(presmoothed_coefPtr, height_ext, length_ext, width_ext);
	reflection3D(prevSolPtr, height_ext, length_ext, width_ext);
	reflection3D(gauss_seidelPtr, height_ext, length_ext, width_ext);

	//perfom presmoothing
	heatExplicitScheme(presmoothingData, filterParameters);

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
				u = presmoothed_coefPtr[k_ext][x_ext];
				uN = presmoothed_coefPtr[k_ext][x_new(i_ext, jminus1, length_ext)];
				uS = presmoothed_coefPtr[k_ext][x_new(i_ext, jplus1, length_ext)];
				uE = presmoothed_coefPtr[k_ext][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				uW = presmoothed_coefPtr[k_ext][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				uNW = presmoothed_coefPtr[k_ext][x_new(iminus1, jminus1, length_ext)];
				uNE = presmoothed_coefPtr[k_ext][x_new(iplus1, jminus1, length_ext)];
				uSE = presmoothed_coefPtr[k_ext][x_new(iplus1, jplus1, length_ext)];
				uSW = presmoothed_coefPtr[k_ext][x_new(iminus1, jplus1, length_ext)];
				Tu = presmoothed_coefPtr[kminus1][x_ext];
				TuN = presmoothed_coefPtr[kminus1][x_new(i_ext, jminus1, length_ext)];
				TuS = presmoothed_coefPtr[kminus1][x_new(i_ext, jplus1, length_ext)];
				TuE = presmoothed_coefPtr[kminus1][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				TuW = presmoothed_coefPtr[kminus1][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				TuNW = presmoothed_coefPtr[kminus1][x_new(iminus1, jminus1, length_ext)];
				TuNE = presmoothed_coefPtr[kminus1][x_new(iplus1, jminus1, length_ext)];
				TuSE = presmoothed_coefPtr[kminus1][x_new(iplus1, jplus1, length_ext)];
				TuSW = presmoothed_coefPtr[kminus1][x_new(iminus1, jplus1, length_ext)];
				Bu = presmoothed_coefPtr[kplus1][x_ext];
				BuN = presmoothed_coefPtr[kplus1][x_new(i_ext, jminus1, length_ext)];
				BuS = presmoothed_coefPtr[kplus1][x_new(i_ext, jplus1, length_ext)];
				BuE = presmoothed_coefPtr[kplus1][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				BuW = presmoothed_coefPtr[kplus1][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				BuNW = presmoothed_coefPtr[kplus1][x_new(iminus1, jminus1, length_ext)];
				BuNE = presmoothed_coefPtr[kplus1][x_new(iplus1, jminus1, length_ext)];
				BuSE = presmoothed_coefPtr[kplus1][x_new(iplus1, jplus1, length_ext)];
				BuSW = presmoothed_coefPtr[kplus1][x_new(iminus1, jplus1, length_ext)];

				//values of voxels in the extended data container for the original image
				orig_u = prevSolPtr[k_ext][x_ext];
				orig_uN = prevSolPtr[k_ext][x_new(i_ext, jminus1, length_ext)];
				orig_uS = prevSolPtr[k_ext][x_new(i_ext, jplus1, length_ext)];
				orig_uE = prevSolPtr[k_ext][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				orig_uW = prevSolPtr[k_ext][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				orig_uNW = prevSolPtr[k_ext][x_new(iminus1, jminus1, length_ext)];
				orig_uNE = prevSolPtr[k_ext][x_new(iplus1, jminus1, length_ext)];
				orig_uSE = prevSolPtr[k_ext][x_new(iplus1, jplus1, length_ext)];
				orig_uSW = prevSolPtr[k_ext][x_new(iminus1, jplus1, length_ext)];
				orig_Tu = prevSolPtr[kminus1][x_ext];
				orig_TuN = prevSolPtr[kminus1][x_new(i_ext, jminus1, length_ext)];
				orig_TuS = prevSolPtr[kminus1][x_new(i_ext, jplus1, length_ext)];
				orig_TuE = prevSolPtr[kminus1][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				orig_TuW = prevSolPtr[kminus1][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				orig_TuNW = prevSolPtr[kminus1][x_new(iminus1, jminus1, length_ext)];
				orig_TuNE = prevSolPtr[kminus1][x_new(iplus1, jminus1, length_ext)];
				orig_TuSE = prevSolPtr[kminus1][x_new(iplus1, jplus1, length_ext)];
				orig_TuSW = prevSolPtr[kminus1][x_new(iminus1, jplus1, length_ext)];
				orig_Bu = prevSolPtr[kplus1][x_ext];
				orig_BuN = prevSolPtr[kplus1][x_new(i_ext, jminus1, length_ext)];
				orig_BuS = prevSolPtr[kplus1][x_new(i_ext, jplus1, length_ext)];
				orig_BuE = prevSolPtr[kplus1][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				orig_BuW = prevSolPtr[kplus1][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				orig_BuNW = prevSolPtr[kplus1][x_new(iminus1, jminus1, length_ext)];
				orig_BuNE = prevSolPtr[kplus1][x_new(iplus1, jminus1, length_ext)];
				orig_BuSE = prevSolPtr[kplus1][x_new(iplus1, jplus1, length_ext)];
				orig_BuSW = prevSolPtr[kplus1][x_new(iminus1, jplus1, length_ext)];

				//calculation of coefficients in the presmooted image data

				// Calculation of coefficients in east direction
				ux = (uE - u) / filterParameters.h;
				uy = ((uN + uNE) - (uS + uSE))
					/ (4.0 * filterParameters.h);
				uz = ((Tu + TuE) - (Bu + BuE))
					/ (4.0 * filterParameters.h);
				presmoot_e_coefPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), coef);

				// Calculation of coefficients in west direction
				ux = (uW - u) / filterParameters.h;
				uy = ((uNW + uN) - (uSW + uS))
					/ (4.0 * filterParameters.h);
				uz = ((TuW + Tu) - (BuW + Bu))
					/ (4.0 * filterParameters.h);
				presmoot_w_coefPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), coef);

				// Calculation of coefficients in north direction
				ux = ((uNE + uE) - (uNW + uW))
					/ (4.0 * filterParameters.h);
				uy = (uN - u) / filterParameters.h;
				uz = ((TuN + Tu) - (BuN + Bu))
					/ (4.0 * filterParameters.h);
				presmoot_n_coefPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), coef);

				// Calculation of coefficients in south direction
				ux = ((uE + uSE) - (uW + uSW))
					/ (4.0 * filterParameters.h);
				uy = (uS - u) / filterParameters.h;
				uz = ((TuS + Tu) - (BuS + Bu))
					/ (4.0 * filterParameters.h);
				presmoot_s_coefPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), coef);

				// Calculation of coefficients in top direction
				ux = ((TuE + uE) - (TuW + uW))
					/ (4.0 * filterParameters.h);
				uy = ((TuN + uN) - (TuS + uS))
					/ (4.0 * filterParameters.h);
				uz = (Tu - u) / filterParameters.h;
				presmoot_t_coefPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), coef);

				// Calculation of coefficients in bottom direction
				ux = ((BuW + uW) - (BuE + uE))
					/ (4.0 * filterParameters.h);
				uy = ((BuN + uN) - (BuS + uS))
					/ (4.0 * filterParameters.h);
				uz = (Bu - u) / filterParameters.h;
				presmoot_b_coefPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), coef);

				//calculation of coefficients in the original image data
				// Calculation of coefficients in east direction
				orig_ux = (orig_uE - orig_u) / filterParameters.h;
				orig_uy = ((orig_uN + orig_uNE) - (orig_uS + orig_uSE))
					/ (4.0 * filterParameters.h);
				orig_uz = ((orig_Tu + orig_TuE) - (orig_Bu + orig_BuE))
					/ (4.0 * filterParameters.h);
				orig_e_coefPtr[k][x] = sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + eps2);

				// Calculation of coefficients in west direction
				orig_ux = (orig_uW - orig_u) / filterParameters.h;
				orig_uy = ((orig_uNW + orig_uN) - (orig_uSW + orig_uS))
					/ (4.0 * filterParameters.h);
				orig_uz = ((orig_TuW + orig_Tu) - (orig_BuW + orig_Bu))
					/ (4.0 * filterParameters.h);
				orig_w_coefPtr[k][x] = sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + eps2);

				// Calculation of coefficients in north direction
				orig_ux = ((orig_uNE + orig_uE) - (orig_uNW + orig_uW))
					/ (4.0 * filterParameters.h);
				orig_uy = (orig_uN - orig_u) / filterParameters.h;
				orig_uz = ((orig_TuN + orig_Tu) - (orig_BuN + orig_Bu))
					/ (4.0 * filterParameters.h);
				orig_n_coefPtr[k][x] = sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + eps2);

				// Calculation of coefficients in south direction
				orig_ux = ((orig_uE + orig_uSE) - (orig_uW + orig_uSW))
					/ (4.0 * filterParameters.h);
				orig_uy = (orig_uS - orig_u) / filterParameters.h;
				orig_uz = ((orig_TuS + orig_Tu) - (orig_BuS + orig_Bu))
					/ (4.0 * filterParameters.h);
				orig_s_coefPtr[k][x] = sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + eps2);

				// Calculation of coefficients in top direction
				orig_ux = ((orig_TuE + orig_uE) - (orig_TuW + orig_uW))
					/ (4.0 * filterParameters.h);
				orig_uy = ((orig_TuN + orig_uN) - (orig_TuS + orig_uS))
					/ (4.0 * filterParameters.h);
				orig_uz = (orig_Tu - orig_u) / filterParameters.h;
				orig_t_coefPtr[k][x] = sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + eps2);

				// Calculation of coefficients in bottom direction
				orig_ux = ((orig_BuW + orig_uW) - (orig_BuE + orig_uE))
					/ (4.0 * filterParameters.h);
				orig_uy = ((orig_BuN + orig_uN) - (orig_BuS + orig_uS))
					/ (4.0 * filterParameters.h);
				orig_uz = (orig_Bu - orig_u) / filterParameters.h;
				orig_b_coefPtr[k][x] = sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + eps2);

				// evaluation of norm of gradient of image at each voxel
				average_face_coef = ((orig_e_coefPtr[k][x] + orig_w_coefPtr[k][x] + orig_n_coefPtr[k][x] + orig_s_coefPtr[k][x]
					+ orig_t_coefPtr[k][x] + orig_b_coefPtr[k][x]) / 6.0);

				voxel_coef = sqrt(pow(average_face_coef, 2) + eps2);

				/* evaluation of norm of gradient of image at each voxel, norm of gradient of presmoothed
				image at each voxel face and reciprocal of norm of gradient of image at each voxel face*/
				coefPtr_e[k][x] = voxel_coef * presmoot_e_coefPtr[k][x] * (1.0 / orig_e_coefPtr[k][x]);//east coefficient
				coefPtr_w[k][x] = voxel_coef * presmoot_w_coefPtr[k][x] * (1.0 / orig_w_coefPtr[k][x]);//west coefficient
				coefPtr_n[k][x] = voxel_coef * presmoot_n_coefPtr[k][x] * (1.0 / orig_n_coefPtr[k][x]);//north coefficient
				coefPtr_s[k][x] = voxel_coef * presmoot_s_coefPtr[k][x] * (1.0 / orig_s_coefPtr[k][x]);//south coefficient
				coefPtr_t[k][x] = voxel_coef * presmoot_t_coefPtr[k][x] * (1.0 / orig_t_coefPtr[k][x]);//top coefficient
				coefPtr_b[k][x] = voxel_coef * presmoot_b_coefPtr[k][x] * (1.0 / orig_b_coefPtr[k][x]);//bottom coefficient
			}
		}
	}

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

					// Begin Gauss-Seidel Formula Evaluation
					gauss_seidel = (prevSolPtr[k_ext][x_ext] + coef_tauh * ((coefPtr_e[k][x] * gauss_seidelPtr[k_ext][x_ext + 1])
						+ (coefPtr_w[k][x] * gauss_seidelPtr[k_ext][x_ext - 1])
						+ (coefPtr_s[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)])
						+ (coefPtr_n[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)])
						+ (coefPtr_b[k][x] * gauss_seidelPtr[k_ext + 1][x_ext])
						+ (coefPtr_t[k][x] * gauss_seidelPtr[k_ext - 1][x_ext]))) /
						(1 + coef_tauh * (coefPtr_e[k][x] + coefPtr_w[k][x] + coefPtr_n[k][x]
							+ coefPtr_s[k][x] + coefPtr_t[k][x] + coefPtr_b[k][x]));

					// SOR implementation using Gauss-Seidel
					gauss_seidelPtr[k_ext][x_ext] = gauss_seidelPtr[k_ext][x_ext] +
						filterParameters.omega_c*(gauss_seidel - gauss_seidelPtr[k_ext][x_ext]);
				}
			}
		}

		// Error Evaluation
		error = 0.0; // Initialize
		for (k = 0, k_ext = 1; k < height; k++, k_ext++)
		{
			for (i = 0, i_ext = 1; i < length; i++, i_ext++)
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++)
				{
					// 2D to 1D representation for i, j
					x_ext = x_new(i_ext, j_ext, length_ext);
					x = x_new(i, j, length);

					error += pow(gauss_seidelPtr[k_ext][x_ext] * (1 + coef_tauh * (coefPtr_e[k][x]
						+ coefPtr_w[k][x] + coefPtr_n[k][x] + coefPtr_s[k][x]
						+ coefPtr_t[k][x] + coefPtr_b[k][x]))
						- coef_tauh * ((coefPtr_e[k][x] * gauss_seidelPtr[k_ext][x_ext + 1])
							+ (coefPtr_w[k][x] * gauss_seidelPtr[k_ext][x_ext - 1])
							+ (coefPtr_s[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)])
							+ (coefPtr_n[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)])
							+ (coefPtr_b[k][x] * gauss_seidelPtr[k_ext + 1][x_ext])
							+ (coefPtr_t[k][x] * gauss_seidelPtr[k_ext - 1][x_ext])) - prevSolPtr[k_ext][x_ext], 2);
				}
			}
		}
	} while (error > filterParameters.tolerance && z < maxNumberOfSolverIteration);
	printf("The number of iterations is %zd\n", z);
	printf("Error is %e\n", error);
	printf("Step is %zd\n", numberOfTimeStep);

	//Copy the current time step to original data holder after timeStepsNum
	copyDataToReducedArea(inputImageData.imageDataPtr, gauss_seidelPtr, height, length, width);

	// Freeing Memory after use
	for (i = 0; i < height_ext; i++)
	{
		free(presmoothed_coefPtr[i]);
		free(gauss_seidelPtr[i]);
		free(prevSolPtr[i]);
	}
	free(presmoothed_coefPtr);
	free(gauss_seidelPtr);
	free(prevSolPtr);

	// free _coefPtr pointers
	for (i = 0; i < height; i++)
	{
		free(presmoot_e_coefPtr[i]);
		free(presmoot_w_coefPtr[i]);
		free(presmoot_n_coefPtr[i]);
		free(presmoot_s_coefPtr[i]);
		free(presmoot_t_coefPtr[i]);
		free(presmoot_b_coefPtr[i]);
	}
	free(presmoot_e_coefPtr);
	free(presmoot_w_coefPtr);
	free(presmoot_n_coefPtr);
	free(presmoot_s_coefPtr);
	free(presmoot_t_coefPtr);
	free(presmoot_b_coefPtr);

	// free orig_ _coefPtr pointers
	for (i = 0; i < height; i++)
	{
		free(orig_e_coefPtr[i]);
		free(orig_w_coefPtr[i]);
		free(orig_n_coefPtr[i]);
		free(orig_s_coefPtr[i]);
		free(orig_t_coefPtr[i]);
		free(orig_b_coefPtr[i]);
	}
	free(orig_e_coefPtr);
	free(orig_w_coefPtr);
	free(orig_n_coefPtr);
	free(orig_s_coefPtr);
	free(orig_t_coefPtr);
	free(orig_b_coefPtr);

	// free coefPtr_ pointers
	for (i = 0; i < height; i++)
	{
		free(coefPtr_e[i]);
		free(coefPtr_w[i]);
		free(coefPtr_n[i]);
		free(coefPtr_s[i]);
		free(coefPtr_t[i]);
		free(coefPtr_b[i]);
	}
	free(coefPtr_e);
	free(coefPtr_w);
	free(coefPtr_n);
	free(coefPtr_s);
	free(coefPtr_t);
	free(coefPtr_b);

	return true;
}