#include <stdlib.h>
#include <stdio.h> // Standard lib for input and output functions
#include <math.h> // Maths functions i.e. pow, sin, cos
#include <stdbool.h> // Boolean function bool
#include "filter_params.h"
#include "non_linear_heat_equation.h"
#include "heat_equation.h"

// Local Function Prototype

bool nonLinearHeatExplicitScheme(Image_Data inputImageData, Filter_Parameters explicitParameters)
{
	size_t k, i, j;
	dataType  hh = explicitParameters.h * explicitParameters.h;
	dataType  tau = explicitParameters.timeStepSize;

	// Perform Reflection of the tempPtr
	// Prepare variables inputImageData.height, inputImageData.length, inputImageData.width
	size_t height = inputImageData.height, length = inputImageData.length, width = inputImageData.width;
	size_t height_ext = height + 2;
	size_t length_ext = length + 2;
	size_t width_ext = width + 2;
	const dataType  coeff = tau / hh;
	dataType  u, uN, uS, uE, uW, uNW, uNE, uSE, uSW, Tu, TuN, TuS, TuE, TuW, TuNW, TuNE, TuSE, TuSW, //current and surrounding voxel values
		Bu, BuN, BuS, BuE, BuW, BuNW, BuNE, BuSE, BuSW;
	size_t kplus1, kminus1, iminus1, iplus1, jminus1, jplus1;
	size_t x;
	size_t x_ext;
	size_t k_ext, j_ext, i_ext;
	dataType  ux, uy, uz;
	dataType  e_coef, w_coef, n_coef, s_coef, t_coef, b_coef, sum_coef;

	Image_Data presmoothingParamenters;
	presmoothingParamenters.height = height_ext;
	presmoothingParamenters.length = length_ext;
	presmoothingParamenters.width = width_ext;

	// Create temporary Image Data holder for Previous time step data - with extended boundary because of boundary condition
	dataType  ** prevSolPtr = (dataType  **)malloc(sizeof(dataType  *) * (height_ext));

	/* Create tempporary Image Data holder for calculation of diffusion coefficients on presmoothed image
	- with extended boundary because of boundary condition*/
	dataType  ** presmoothed_coefPtr = (dataType  **)malloc(sizeof(dataType  *) * (height_ext));

	//checks if the memory was allocated
	if (prevSolPtr == NULL || presmoothed_coefPtr == NULL)
		return false;

	for (k = 0; k < height_ext; k++)
	{
		presmoothed_coefPtr[k] = malloc(sizeof(dataType)*(length_ext)*(width_ext));
		prevSolPtr[k] = malloc(sizeof(dataType)*(length_ext)*(width_ext));
		//checks if the memory was allocated
		if (presmoothed_coefPtr[k] == NULL || prevSolPtr[k] == NULL)
			return false;
	}

	presmoothingParamenters.imageDataPtr = presmoothed_coefPtr;

	//copy data to extended area which will be used in each time step
	copyDataToExtendedArea(inputImageData.imageDataPtr, prevSolPtr, height, length, width);

	//copy data to extended area which will be used for presmoothing
	copyDataToExtendedArea(inputImageData.imageDataPtr, presmoothed_coefPtr, height, length, width);

	//perform reflection of the extended area to solve border problem
	reflection3D(presmoothed_coefPtr, height_ext, length_ext, width_ext);
	reflection3D(prevSolPtr, height_ext, length_ext, width_ext);

	//perfom presmoothing
	heatExplicitScheme(presmoothingParamenters, explicitParameters);

	// The Explicit Scheme Evaluation
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

				// Calculation of coefficients in east direction
				ux = (uE - u) / explicitParameters.h;
				uy = ((uN + uNE) - (uS + uSE))
					/ (4.0 * explicitParameters.h);
				uz = ((Tu + TuE) - (Bu + BuE))
					/ (4.0 * explicitParameters.h);
				e_coef = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), explicitParameters.edge_detector_coefficient);

				// Calculation of coefficients in west direction
				ux = (uW - u) / explicitParameters.h;
				uy = ((uNW + uN) - (uSW + uS))
					/ (4.0 * explicitParameters.h);
				uz = ((TuW + Tu) - (BuW + Bu))
					/ (4.0 * explicitParameters.h);
				w_coef = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), explicitParameters.edge_detector_coefficient);

				// Calculation of coefficients in north direction
				ux = ((uNE + uE) - (uNW + uW))
					/ (4.0 * explicitParameters.h);
				uy = (uN - u) / explicitParameters.h;
				uz = ((TuN + Tu) - (BuN + Bu))
					/ (4.0 * explicitParameters.h);
				n_coef = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), explicitParameters.edge_detector_coefficient);

				// Calculation of coefficients in south direction
				ux = ((uE + uSE) - (uW + uSW))
					/ (4.0 * explicitParameters.h);
				uy = (uS - u) / explicitParameters.h;
				uz = ((TuS + Tu) - (BuS + Bu))
					/ (4.0 * explicitParameters.h);
				s_coef = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), explicitParameters.edge_detector_coefficient);

				// Calculation of coefficients in top direction
				ux = ((TuE + uE) - (TuW + uW))
					/ (4.0 * explicitParameters.h);
				uy = ((TuN + uN) - (TuS + uS))
					/ (4.0 * explicitParameters.h);
				uz = (Tu - u) / explicitParameters.h;
				t_coef = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), explicitParameters.edge_detector_coefficient);

				// Calculation of coefficients in bottom direction
				ux = ((BuW + uW) - (BuE + uE))
					/ (4.0 * explicitParameters.h);
				uy = ((BuN + uN) - (BuS + uS))
					/ (4.0 * explicitParameters.h);
				uz = (Bu - u) / explicitParameters.h;
				b_coef = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), explicitParameters.edge_detector_coefficient);

				// Sum of all coefficients
				sum_coef = e_coef + w_coef + n_coef + s_coef + t_coef + b_coef;

				// Explicit formula
				inputImageData.imageDataPtr[k][x] = (1.0 - coeff * (sum_coef))*prevSolPtr[k_ext][x_ext]
					+ coeff * ((e_coef * prevSolPtr[k_ext][x_ext + 1])
						+ (w_coef * prevSolPtr[k_ext][x_ext - 1])
						+ (s_coef * prevSolPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)])
						+ (n_coef * prevSolPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)])
						+ (b_coef * prevSolPtr[k_ext + 1][x_ext])
						+ (t_coef * prevSolPtr[k_ext - 1][x_ext]));
			}
		}
	}

	// Freeing Memory after use
	for (i = 0; i < height_ext; i++)
	{
		free(prevSolPtr[i]);
		free(presmoothed_coefPtr[i]);
	}
	free(prevSolPtr);
	free(presmoothed_coefPtr);

	return true;
}