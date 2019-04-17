#include <stdio.h> // Standard lib for input and output functions
#include <stdlib.h>
#include <math.h> // Maths functions i.e. pow, sin, cos
#include <stdbool.h> // Boolean function bool
#include "heat_equation.h"
#include "non_linear_heat_equation.h"
#include "filter_params.h"

// Local Function Prototype

bool nonLinearHeatImplicitScheme(Image_Data inputImageData, Filter_Parameters implicitParameters, size_t numberOfTimeStep)
{
	size_t k, i, j;
	dataType  hh = implicitParameters.h * implicitParameters.h;
	dataType  tau = 100 * implicitParameters.timeStepSize;

	// Error value used to check iteration
	// sor - successive over relation value, used in Gauss-Seidel formula
	dataType  error, gauss_seidel;

	// Perform Reflection of the tempPtr
	// Prepare variables toExplicitImage.height, toExplicitImage.length, toExplicitImage.width
	size_t height = inputImageData.height, length = inputImageData.length, width = inputImageData.width;
	size_t height_ext = height + 2;
	size_t length_ext = length + 2;
	size_t width_ext = width + 2;
	size_t x; //x = x_new(i, j, length);
	size_t x_ext; //x_ext = x_new(i_ext, j_ext, length_ext);
	size_t z; // Steps counter
	const dataType  coeff = tau / hh;
	dataType  u, uN, uS, uE, uW, uNW, uNE, uSE, uSW, Tu, TuN, TuS, TuE, TuW, TuNW, TuNE, TuSE, TuSW, //current and surrounding voxel values
		Bu, BuN, BuS, BuE, BuW, BuNW, BuNE, BuSE, BuSW;
	size_t kplus1, kminus1, iminus1, iplus1, jminus1, jplus1;
	size_t k_ext, j_ext, i_ext;
	dataType  ux, uy, uz;

	Image_Data presmoothingData;
	presmoothingData.height = height_ext;
	presmoothingData.length = length_ext;
	presmoothingData.width = width_ext;

	// Create temporary Image Data holder for Previous time step data - with extended boundary because of boundary condition
	dataType  ** prevSolPtr = (dataType  **)malloc(sizeof(dataType  *) * (height_ext));

	// Create temporary Image Data holder for Current time step data - with extended boundary because of boundary condition
	dataType  ** gauss_seidelPtr = (dataType  **)malloc(sizeof(dataType  *) * (height_ext));

	/* Create tempporary Image Data holder for calculation of diffusion coefficients on presmoothed image
	- with extended boundary because of boundary condition*/
	dataType  ** presmoothed_coefPtr = (dataType  **)malloc(sizeof(dataType  *) * (height_ext));

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
	dataType  ** e_coefPtr = (dataType  **)malloc(sizeof(dataType  *) * (height));
	dataType  ** w_coefPtr = (dataType  **)malloc(sizeof(dataType  *) * (height));
	dataType  ** n_coefPtr = (dataType  **)malloc(sizeof(dataType  *) * (height));
	dataType  ** s_coefPtr = (dataType  **)malloc(sizeof(dataType  *) * (height));
	dataType  ** t_coefPtr = (dataType  **)malloc(sizeof(dataType  *) * (height));
	dataType  ** b_coefPtr = (dataType  **)malloc(sizeof(dataType  *) * (height));

	//checks if the memory was allocated
	if (e_coefPtr == NULL || w_coefPtr == NULL || n_coefPtr == NULL || s_coefPtr == NULL || t_coefPtr == NULL || b_coefPtr == NULL)
		return false;

	for (k = 0; k < height; k++)
	{
		e_coefPtr[k] = malloc(sizeof(dataType)*(length)*(height));
		w_coefPtr[k] = malloc(sizeof(dataType)*(length)*(height));
		n_coefPtr[k] = malloc(sizeof(dataType)*(length)*(height));
		s_coefPtr[k] = malloc(sizeof(dataType)*(length)*(height));
		t_coefPtr[k] = malloc(sizeof(dataType)*(length)*(height));
		b_coefPtr[k] = malloc(sizeof(dataType)*(length)*(height));

		//checks if the memory was allocated
		if (e_coefPtr[k] == NULL || w_coefPtr[k] == NULL || n_coefPtr[k] == NULL || s_coefPtr[k] == NULL || t_coefPtr[k] == NULL || b_coefPtr[k] == NULL)
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

	//perform reflection of the extended area to solve border problem
	reflection3D(prevSolPtr, height_ext, length_ext, width_ext);
	reflection3D(gauss_seidelPtr, height_ext, length_ext, width_ext);

	//perfom presmoothing
	heatExplicitScheme(presmoothingData, implicitParameters);

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

				//calculation of coefficients in the presmooted image data

				// Calculation of coefficients in east direction
				ux = (uE - u) / implicitParameters.h;
				uy = ((uN + uNE) - (uS + uSE))
					/ (4.0 * implicitParameters.h);
				uz = ((Tu + TuE) - (Bu + BuE))
					/ (4.0 * implicitParameters.h);
				e_coefPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), implicitParameters.edge_detector_coefficient);

				// Calculation of coefficients in west direction
				ux = (uW - u) / implicitParameters.h;
				uy = ((uNW + uN) - (uSW + uS))
					/ (4.0 * implicitParameters.h);
				uz = ((TuW + Tu) - (BuW + Bu))
					/ (4.0 * implicitParameters.h);
				w_coefPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), implicitParameters.edge_detector_coefficient);

				// Calculation of coefficients in north direction
				ux = ((uNE + uE) - (uNW + uW))
					/ (4.0 * implicitParameters.h);
				uy = (uN - u) / implicitParameters.h;
				uz = ((TuN + Tu) - (BuN + Bu))
					/ (4.0 * implicitParameters.h);
				n_coefPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), implicitParameters.edge_detector_coefficient);

				// Calculation of coefficients in south direction
				ux = ((uE + uSE) - (uW + uSW))
					/ (4.0 * implicitParameters.h);
				uy = (uS - u) / implicitParameters.h;
				uz = ((TuS + Tu) - (BuS + Bu))
					/ (4.0 * implicitParameters.h);
				s_coefPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), implicitParameters.edge_detector_coefficient);

				// Calculation of coefficients in top direction
				ux = ((TuE + uE) - (TuW + uW))
					/ (4.0 * implicitParameters.h);
				uy = ((TuN + uN) - (TuS + uS))
					/ (4.0 * implicitParameters.h);
				uz = (Tu - u) / implicitParameters.h;
				t_coefPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), implicitParameters.edge_detector_coefficient);

				// Calculation of coefficients in bottom direction
				ux = ((BuW + uW) - (BuE + uE))
					/ (4.0 * implicitParameters.h);
				uy = ((BuN + uN) - (BuS + uS))
					/ (4.0 * implicitParameters.h);
				uz = (Bu - u) / implicitParameters.h;
				b_coefPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), implicitParameters.edge_detector_coefficient);
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
					gauss_seidel = (prevSolPtr[k_ext][x_ext] + coeff * ((e_coefPtr[k][x] * gauss_seidelPtr[k_ext][x_ext + 1])
						+ (w_coefPtr[k][x] * gauss_seidelPtr[k_ext][x_ext - 1])
						+ (s_coefPtr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)])
						+ (n_coefPtr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)])
						+ (b_coefPtr[k][x] * gauss_seidelPtr[k_ext + 1][x_ext])
						+ (t_coefPtr[k][x] * gauss_seidelPtr[k_ext - 1][x_ext]))) /
						(1 + coeff * (e_coefPtr[k][x] + w_coefPtr[k][x] + n_coefPtr[k][x]
							+ s_coefPtr[k][x] + t_coefPtr[k][x] + b_coefPtr[k][x]));

					// SOR implementation using Gauss-Seidel
					gauss_seidelPtr[k_ext][x_ext] = gauss_seidelPtr[k_ext][x_ext] +
						implicitParameters.omega_c*(gauss_seidel - gauss_seidelPtr[k_ext][x_ext]);
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

					error += pow(gauss_seidelPtr[k_ext][x_ext] * (1 + coeff * (e_coefPtr[k][x]
						+ w_coefPtr[k][x] + n_coefPtr[k][x] + s_coefPtr[k][x]
						+ t_coefPtr[k][x] + b_coefPtr[k][x]))
						- coeff * ((e_coefPtr[k][x] * gauss_seidelPtr[k_ext][x_ext + 1])
							+ (w_coefPtr[k][x] * gauss_seidelPtr[k_ext][x_ext - 1])
							+ (s_coefPtr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)])
							+ (n_coefPtr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)])
							+ (b_coefPtr[k][x] * gauss_seidelPtr[k_ext + 1][x_ext])
							+ (t_coefPtr[k][x] * gauss_seidelPtr[k_ext - 1][x_ext])) - prevSolPtr[k_ext][x_ext], 2);
				}
			}
		}
	} while (error > implicitParameters.tolerance && z < implicitParameters.maxNumberOfSolverIteration);
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

	for (i = 0; i < height; i++)
	{
		free(e_coefPtr[i]);
		free(w_coefPtr[i]);
		free(n_coefPtr[i]);
		free(s_coefPtr[i]);
		free(t_coefPtr[i]);
		free(b_coefPtr[i]);
	}
	free(e_coefPtr);
	free(w_coefPtr);
	free(n_coefPtr);
	free(s_coefPtr);
	free(t_coefPtr);
	free(b_coefPtr);

	return true;
}