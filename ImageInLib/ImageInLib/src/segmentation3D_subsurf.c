/*
* Author: Markjoe Olunna UBA
* Purpose: ImageInLife project - 4D Image Segmentation Methods
* Language:  C
*/
#include <stdio.h> // Standard lib for input and output functions
#include <stdlib.h>
#include <time.h>
#include <math.h> // Maths functions i.e. pow, sin, cos
#include <stdbool.h> // Boolean function bool
#include <string.h>
#include "heat_equation.h"
#include "non_linear_heat_equation.h"
#include "segmentation3D_subsurf.h"
#include "image_norm.h"
#include "data_initialization.h"
#include "edgedetection.h"
#include "data_storage.h"
//#include "data_load.h"
#include "generate_3D_shapes.h"
#include "common_functions.h"
#include "Common_Math.h"
#include "setting_boundary_values.h"
#include "ctype.h"
#include "filter_params.h"
#include "vtk_params.h"
// Local Function Prototype

bool subsurfSegmentation(Image_Data inputImageData, Segmentation_Parameters segParameters, Filter_Parameters explicit_lhe_Parameters,
	Point3D * centers, size_t no_of_centers, unsigned char * outputPathPtr)//bool subsurfSegmentation()
{
	//const size_t length = 50, width = 57, height = 20;// length = 50, width = 57, height = 20;//length = 101, width = 101, height = 101
	size_t i; // length == xDim, width == yDim, height == zDim
	size_t dim2D = inputImageData.length * inputImageData.width;
	size_t height = inputImageData.height;
	size_t length = inputImageData.length;
	size_t width = inputImageData.width;
	size_t height_ext = inputImageData.height + 2;
	size_t length_ext = inputImageData.length + 2;
	size_t width_ext = inputImageData.width + 2;
	size_t dim2D_ext = length_ext * width_ext;
	dataType firstCpuTime, secondCpuTime, difference_btw_current_and_previous_sol;
	dataType ** gauss_seidelPtr = (dataType **)malloc(sizeof(dataType *) * height_ext);
	dataType ** prevSol_extPtr = (dataType **)malloc(sizeof(dataType *) * height_ext);
	dataType ** imageToBeSegPtr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** segmFuntionPtr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** GePtr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** GwPtr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** GnPtr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** GsPtr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** GtPtr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** GbPtr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** e_Ptr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** w_Ptr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** n_Ptr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** s_Ptr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** t_Ptr = (dataType **)malloc(sizeof(dataType *) * height);
	dataType ** b_Ptr = (dataType **)malloc(sizeof(dataType *) * height);

	//checks if the memory was allocated
	if (gauss_seidelPtr == NULL || prevSol_extPtr == NULL || imageToBeSegPtr == NULL || segmFuntionPtr == NULL || GePtr == NULL ||
		GwPtr == NULL || GnPtr == NULL || GsPtr == NULL || GtPtr == NULL || GbPtr == NULL || e_Ptr == NULL ||
		w_Ptr == NULL || n_Ptr == NULL || s_Ptr == NULL || t_Ptr == NULL || b_Ptr == NULL)
		return false;
	for (i = 0; i < height; i++)
	{
		imageToBeSegPtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
		segmFuntionPtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
		GePtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
		GwPtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
		GnPtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
		GsPtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
		GtPtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
		GbPtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
		e_Ptr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
		w_Ptr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
		n_Ptr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
		s_Ptr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
		t_Ptr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
		b_Ptr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
		//checks if the memory was allocated
		if (imageToBeSegPtr[i] == NULL || segmFuntionPtr[i] == NULL || GePtr[i] == NULL
			|| GwPtr[i] == NULL || GnPtr[i] == NULL || GsPtr[i] == NULL || GtPtr[i] == NULL || GbPtr[i] == NULL || e_Ptr[i] == NULL
			|| w_Ptr[i] == NULL || n_Ptr[i] == NULL || s_Ptr[i] == NULL || t_Ptr[i] == NULL || b_Ptr[i] == NULL)
			return false;
	}
	for (i = 0; i < height_ext; i++)
	{
		gauss_seidelPtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D_ext);
		prevSol_extPtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D_ext);
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

	Gradient_Pointers GPtrs;
	GPtrs.GePtr = GePtr;
	GPtrs.GwPtr = GwPtr;
	GPtrs.GnPtr = GnPtr;
	GPtrs.GsPtr = GsPtr;
	GPtrs.GtPtr = GtPtr;
	GPtrs.GbPtr = GbPtr;

	Coefficient_Pointers CoefPtrs;
	CoefPtrs.e_Ptr = e_Ptr;
	CoefPtrs.w_Ptr = w_Ptr;
	CoefPtrs.n_Ptr = n_Ptr;
	CoefPtrs.s_Ptr = s_Ptr;
	CoefPtrs.t_Ptr = t_Ptr;
	CoefPtrs.b_Ptr = b_Ptr;

	//generate initial segmentation function
	generateInitialSegmentationFunctionForMultipleCentres(segmFuntionPtr, length, width, height, centers, 0.5, 10, no_of_centers);

	//compute coefficients from presmoothed image
	gFunctionForImageToBeSegmented(inputImageData, prevSol_extPtr, GPtrs, segParameters, explicit_lhe_Parameters);

	//Array for name construction
	unsigned char name[350];
	unsigned char name_ending[100];

	//loop for segmentation time steps
	i = 1;
	do
	{
		segParameters.numberOfTimeStep = i;
		firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);

		//calcution of coefficients
		gaussSeidelCoefficients(prevSol_extPtr, imageData, GPtrs, CoefPtrs, segParameters);

		// Call to function that will evolve segmentation function in each discrete time step
		subsurfSegmentationTimeStep(prevSol_extPtr, gauss_seidelPtr, imageData, GPtrs, segParameters, CoefPtrs, centers, no_of_centers);

		secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);

		//Compute the L2 norm of the difference between the current and previous solutions
		difference_btw_current_and_previous_sol = l2normD(prevSol_extPtr, gauss_seidelPtr, length, width, height, segParameters.h);

		printf("mass is %e\n", difference_btw_current_and_previous_sol);
		printf("segTolerance is %e\n", segParameters.segTolerance);
		printf("CPU time: %e secs\n", secondCpuTime - firstCpuTime);

		//writing density.
		if ((i%segParameters.mod) == 0)
		{
			strcpy_s(name, sizeof name, outputPathPtr);
			sprintf_s(name_ending, sizeof(name_ending), "_seg_func_%03zd.vtk", i);
			strcat_s(name, sizeof(name), name_ending);
			Storage_Flags flags = { true, true };
			store3dDataVtkD(imageData.segmentationFuntionPtr, length, width, height, name, segParameters.h, flags);
		}
		i++;
	} while ((i <= segParameters.maxNoOfTimeSteps) && (difference_btw_current_and_previous_sol > segParameters.segTolerance));

	printf("finish: Segmentation tolerance is %lf\n", segParameters.segTolerance);

	for (i = 0; i < height; i++)
	{
		free(imageToBeSegPtr[i]);
		free(segmFuntionPtr[i]);
		free(GePtr[i]);
		free(GwPtr[i]);
		free(GnPtr[i]);
		free(GsPtr[i]);
		free(GtPtr[i]);
		free(GbPtr[i]);
		free(e_Ptr[i]);
		free(w_Ptr[i]);
		free(n_Ptr[i]);
		free(s_Ptr[i]);
		free(t_Ptr[i]);
		free(b_Ptr[i]);
	}
	free(imageToBeSegPtr);
	free(segmFuntionPtr);
	free(GePtr);
	free(GwPtr);
	free(GnPtr);
	free(GsPtr);
	free(GtPtr);
	free(GbPtr);
	free(e_Ptr);
	free(w_Ptr);
	free(n_Ptr);
	free(s_Ptr);
	free(t_Ptr);
	free(b_Ptr);

	for (i = 0; i < height_ext; i++)
	{
		free(prevSol_extPtr[i]);
		free(gauss_seidelPtr[i]);
	}
	free(prevSol_extPtr);
	free(gauss_seidelPtr);

	return true;
}

bool subsurfSegmentationTimeStep(dataType **prevSol_extPtr, dataType **gauss_seidelPtr, Segment_Image_Data inputImageData, Gradient_Pointers GPtrs,
	Segmentation_Parameters segParameters, Coefficient_Pointers CoefPtrs, Point3D * centers, size_t no_of_centers)
{
	//check if the memory was allocated successfully
	if (inputImageData.segmentationFuntionPtr == NULL || prevSol_extPtr == NULL || gauss_seidelPtr == NULL || GPtrs.GePtr == NULL || GPtrs.GwPtr == NULL
		|| GPtrs.GnPtr == NULL || GPtrs.GsPtr == NULL || GPtrs.GtPtr == NULL || GPtrs.GbPtr == NULL || CoefPtrs.e_Ptr == NULL || CoefPtrs.w_Ptr == NULL
		|| CoefPtrs.n_Ptr == NULL || CoefPtrs.s_Ptr == NULL || CoefPtrs.t_Ptr == NULL || CoefPtrs.b_Ptr == NULL)
		return false;

	size_t k, i, j;
	dataType hh = segParameters.h * segParameters.h, hhh = segParameters.h * segParameters.h * segParameters.h;
	dataType tau = segParameters.tau;

	// Error value used to check iteration
	// sor - successive over relation value, used in Gauss-Seidel formula
	dataType mean_square_residue, gauss_seidel;

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

	//copy data to extended area which will be used in each Gauss Seidel iteration
	copyDataToExtendedArea(inputImageData.segmentationFuntionPtr, gauss_seidelPtr, height, length, width);
	copyDataToExtendedArea(inputImageData.segmentationFuntionPtr, prevSol_extPtr, height, length, width);

	//set boundary values to ensure Zero Dirichlet boundary condition.
	setBoundaryToZeroDirichletBC(gauss_seidelPtr, length_ext, width_ext, height_ext);
	setBoundaryToZeroDirichletBC(prevSol_extPtr, length_ext, width_ext, height_ext);

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
					gauss_seidel = (prevSol_extPtr[k_ext][x_ext] + coef_tauh * ((CoefPtrs.e_Ptr[k][x] * gauss_seidelPtr[k_ext][x_ext + 1])
						+ (CoefPtrs.w_Ptr[k][x] * gauss_seidelPtr[k_ext][x_ext - 1])
						+ (CoefPtrs.s_Ptr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)])
						+ (CoefPtrs.n_Ptr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)])
						+ (CoefPtrs.b_Ptr[k][x] * gauss_seidelPtr[k_ext + 1][x_ext])
						+ (CoefPtrs.t_Ptr[k][x] * gauss_seidelPtr[k_ext - 1][x_ext]))) /
						(1 + coef_tauh * (CoefPtrs.e_Ptr[k][x] + CoefPtrs.w_Ptr[k][x] + CoefPtrs.n_Ptr[k][x]
							+ CoefPtrs.s_Ptr[k][x] + CoefPtrs.t_Ptr[k][x] + CoefPtrs.b_Ptr[k][x]));

					// SOR implementation using Gauss-Seidel
					gauss_seidelPtr[k_ext][x_ext] = gauss_seidelPtr[k_ext][x_ext] +
						segParameters.omega_c*(gauss_seidel - gauss_seidelPtr[k_ext][x_ext]);
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

					mean_square_residue += (pow(gauss_seidelPtr[k_ext][x_ext] * (1 + coef_tauh * (CoefPtrs.e_Ptr[k][x]
						+ CoefPtrs.w_Ptr[k][x] + CoefPtrs.n_Ptr[k][x] + CoefPtrs.s_Ptr[k][x]
						+ CoefPtrs.t_Ptr[k][x] + CoefPtrs.b_Ptr[k][x]))
						- coef_tauh * ((CoefPtrs.e_Ptr[k][x] * gauss_seidelPtr[k_ext][x_ext + 1])
							+ (CoefPtrs.w_Ptr[k][x] * gauss_seidelPtr[k_ext][x_ext - 1])
							+ (CoefPtrs.s_Ptr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)])
							+ (CoefPtrs.n_Ptr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)])
							+ (CoefPtrs.b_Ptr[k][x] * gauss_seidelPtr[k_ext + 1][x_ext])
							+ (CoefPtrs.t_Ptr[k][x] * gauss_seidelPtr[k_ext - 1][x_ext])) - prevSol_extPtr[k_ext][x_ext], 2) * hhh);
				}
			}
		}
	} while (mean_square_residue > segParameters.gauss_seidelTolerance && z < segParameters.maxNoGSIteration);
	printf("The number of iterations is %zd\n", z);
	printf("Residuum is %e\n", mean_square_residue);
	printf("Step is %zd\n", segParameters.numberOfTimeStep);

	//Copy the current time step to original data array after timeStepsNum
	copyDataToReducedArea(inputImageData.segmentationFuntionPtr, gauss_seidelPtr, height, length, width);

	//Rescale values of segmentation function and current time step to interval (0, 1)
	if (no_of_centers == 1)
	{
		rescaleToIntervalZeroOne(inputImageData.segmentationFuntionPtr, length, width, height);
		rescaleToIntervalZeroOne(gauss_seidelPtr, length_ext, width_ext, height_ext);
	}
	else
	{
		for (i = 0; i < no_of_centers; i++)
		{
			rescaleLocallyToIntervalZeroOne(inputImageData.segmentationFuntionPtr, length, width, height,
				centers[i].x, centers[i].y, centers[i].z, 6., i);
			rescaleLocallyToIntervalZeroOne(gauss_seidelPtr, length_ext, width_ext, height_ext,
				centers[i].x, centers[i].y, centers[i].z, 6., i);
		}
	}

	return true;
}
bool rescaleToIntervalZeroOne(dataType **imagePtr, size_t length, size_t width, size_t height)
{
	//check if the memory was allocated successfully
	if (imagePtr == NULL)
		return false;

	size_t k, i, j;
	dataType max = 0, min = 100000, quotient, offset;
	size_t x; //x = x_new(i, j, length);

			  //Determine minimum and maximum value
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				// 2D to 1D representation for i, j
				x = x_new(i, j, length);
				if (imagePtr[k][x] < min)
					min = imagePtr[k][x];
				if (imagePtr[k][x] > max)
					max = imagePtr[k][x];
			}
		}
	}
	quotient = 1. / (max - min);
	offset = min * quotient;
	//Rescale values to interval (0, 1)
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				// 2D to 1D representation for i, j
				x = x_new(i, j, length);

				imagePtr[k][x] = (quotient * imagePtr[k][x]) - offset;
			}
		}
	}
	return true;
}

bool rescaleLocallyToIntervalZeroOne(dataType **imagePtr, size_t length, size_t width, size_t height,
	dataType center_x, dataType center_y, dataType center_z, dataType R, size_t counter)
{
	//check if the memory was allocated successfully
	if (imagePtr == NULL)
		return false;

	size_t k, i, j;
	dataType max = 0, min = 100000, quotient, offset, dz, dy, dx, norm_of_distance, new_value;
	size_t x; //x = x_new(i, j, length);

			  //Determine minimum and maximum value
	for (k = 0; k < height; k++)
	{
		dz = k - center_z;
		for (i = 0; i < length; i++)
		{
			dx = i - center_x;
			for (j = 0; j < width; j++)
			{
				dy = j - center_y;
				// 2D to 1D representation for i, j
				x = x_new(i, j, length);

				//Find local minimum and maximum
				norm_of_distance = sqrt((dx * dx) + (dy * dy) + (dz * dz));
				if (norm_of_distance <= R)
				{
					if (imagePtr[k][x] < min)
						min = imagePtr[k][x];
					if (imagePtr[k][x] > max)
						max = imagePtr[k][x];
				}
			}
		}
	}
	quotient = 1. / (max - min);
	offset = min * quotient;
	//Rescale values to interval (0, 1)
	for (k = 0; k < height; k++)
	{
		dz = k - center_z;
		for (i = 0; i < length; i++)
		{
			dx = i - center_x;
			for (j = 0; j < width; j++)
			{
				dy = j - center_y;
				// 2D to 1D representation for i, j
				x = x_new(i, j, length);

				// recaling values
				norm_of_distance = sqrt((dx * dx) + (dy * dy) + (dz * dz));
				new_value = (quotient * imagePtr[k][x]) - offset;

				if (norm_of_distance <= R)
					imagePtr[k][x] = new_value;
			}
		}
	}
	return true;
}

bool generateInitialSegmentationFunctionForMultipleCentres(dataType **inputDataArrayPtr, size_t length, size_t width, size_t height,
	Point3D *centers, dataType v, dataType R, size_t no_of_centers)
{
	size_t i, j, k, s;//loop counter for z dimension
	dataType dx, dy, dz, norm_of_distance, new_value;
	//checks if the memory was allocated
	if (inputDataArrayPtr == NULL)
		return false;
	//Storage paths
	unsigned char pathArray1[] = "D:\\segmentation\\test for library release\\segFunction.vtk";

	// Construction of segmentation function
	for (s = 0; s < no_of_centers; s++)
	{
		for (k = 0; k < height; k++)// zDim
		{
			dz = k - centers[s].z;
			for (i = 0; i < length; i++)//xDim
			{
				dx = i - centers[s].x;
				for (j = 0; j < width; j++) //yDim
				{
					dy = j - centers[s].y;
					// 1D representation
					size_t x_n = x_new(i, j, length);
					// Set Value
					norm_of_distance = sqrt((dx * dx) + (dy * dy) + (dz * dz));
					new_value = (1. / (sqrt((dx * dx) + (dy * dy) + (dz * dz)) + v)) - (1. / (R + v));
					if (s == 0)
					{
						if (norm_of_distance > R)
							inputDataArrayPtr[k][x_n] = 0;
						else
							inputDataArrayPtr[k][x_n] = new_value;
					}
					else
					{
						if (norm_of_distance <= R)
							if (inputDataArrayPtr[k][x_n] < new_value)
								inputDataArrayPtr[k][x_n] = new_value;
					}
				}
			}
		}

		if (no_of_centers == 1)
		{
			rescaleToIntervalZeroOne(inputDataArrayPtr, length, width, height);
		}
		else
		{
			for (i = 0; i < no_of_centers; i++)
			{
				rescaleLocallyToIntervalZeroOne(inputDataArrayPtr, length, width, height, centers[s].x, centers[s].y, centers[s].z, 6., s);
			}
		}
	}
	Storage_Flags flags = { true, true };
	store3dDataVtkD(inputDataArrayPtr, length, width, height, pathArray1, (2.5 / (length)), flags);
	return true;
}

bool gFunctionForImageToBeSegmented(Image_Data inputImageData, dataType **extendedCoefPtr, Gradient_Pointers GPtrs,
	Segmentation_Parameters segParameters, Filter_Parameters explicit_lhe_Parameters)
{
	//checks if the memory was allocated
	if (inputImageData.imageDataPtr == NULL || extendedCoefPtr == NULL || GPtrs.GePtr == NULL || GPtrs.GwPtr == NULL
		|| GPtrs.GnPtr == NULL || GPtrs.GsPtr == NULL || GPtrs.GtPtr == NULL || GPtrs.GbPtr == NULL)
		return false;

	dataType **imageToBeSegPtr = inputImageData.imageDataPtr;
	size_t i, j, k, x, x_ext; // length == xDim, width == yDim, height == zDim
	size_t kplus1, kminus1, iminus1, iplus1, jminus1, jplus1;
	size_t dim2D = inputImageData.length * inputImageData.width;
	size_t k_ext, j_ext, i_ext;
	size_t height_ext = inputImageData.height + 2;
	size_t length_ext = inputImageData.length + 2;
	size_t width_ext = inputImageData.width + 2;
	dataType quotient = 4.0 * segParameters.h;
	dataType inverse = 1. / (VTK_MAX_HEADER_LINE_LENGTH - 1);
	dataType ux, uy, uz; //change in x, y and z respectively
	dataType u, uN, uS, uE, uW, uNW, uNE, uSE, uSW, Tu, TuN, TuS, TuE, TuW, TuNW, TuNE, TuSE, TuSW, //current and surrounding voxel values
		Bu, BuN, BuS, BuE, BuW, BuNW, BuNE, BuSE, BuSW;

	Image_Data presmoothingData;
	presmoothingData.height = height_ext;
	presmoothingData.length = length_ext;
	presmoothingData.width = width_ext;
	presmoothingData.imageDataPtr = extendedCoefPtr;

	//copy data to extended area which will be used for calculation of diffusion coefficients
	copyDataToExtendedArea(imageToBeSegPtr, extendedCoefPtr, inputImageData.height, inputImageData.length, inputImageData.width);

	//perform reflection of the extended area to ensure zero Neumann boundary condition (for LHE)
	reflection3D(extendedCoefPtr, height_ext, length_ext, width_ext);

	//perfom presmoothing
	heatExplicitScheme(presmoothingData, explicit_lhe_Parameters);
	copyDataToReducedArea(imageToBeSegPtr, extendedCoefPtr, inputImageData.height, inputImageData.length, inputImageData.width);

	//calculation of coefficients
	for (k = 0, k_ext = 1; k < inputImageData.height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < inputImageData.length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < inputImageData.width; j++, j_ext++)
			{
				// 2D to 1D representation for i, j
				x_ext = x_new(i_ext, j_ext, length_ext);
				x = x_new(i, j, inputImageData.length);
				iminus1 = i_ext - 1;
				iplus1 = i_ext + 1;
				jplus1 = j_ext + 1;
				jminus1 = j_ext - 1;
				kplus1 = k_ext + 1;
				kminus1 = k_ext - 1;

				//values of voxels in the extended data container for presmoothed image
				u = inverse * extendedCoefPtr[k_ext][x_ext];
				uN = inverse * extendedCoefPtr[k_ext][x_new(i_ext, jminus1, length_ext)];
				uS = inverse * extendedCoefPtr[k_ext][x_new(i_ext, jplus1, length_ext)];
				uE = inverse * extendedCoefPtr[k_ext][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				uW = inverse * extendedCoefPtr[k_ext][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				uNW = inverse * extendedCoefPtr[k_ext][x_new(iminus1, jminus1, length_ext)];
				uNE = inverse * extendedCoefPtr[k_ext][x_new(iplus1, jminus1, length_ext)];
				uSE = inverse * extendedCoefPtr[k_ext][x_new(iplus1, jplus1, length_ext)];
				uSW = inverse * extendedCoefPtr[k_ext][x_new(iminus1, jplus1, length_ext)];
				Tu = inverse * extendedCoefPtr[kminus1][x_ext];
				TuN = inverse * extendedCoefPtr[kminus1][x_new(i_ext, jminus1, length_ext)];
				TuS = inverse * extendedCoefPtr[kminus1][x_new(i_ext, jplus1, length_ext)];
				TuE = inverse * extendedCoefPtr[kminus1][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				TuW = inverse * extendedCoefPtr[kminus1][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				TuNW = inverse * extendedCoefPtr[kminus1][x_new(iminus1, jminus1, length_ext)];
				TuNE = inverse * extendedCoefPtr[kminus1][x_new(iplus1, jminus1, length_ext)];
				TuSE = inverse * extendedCoefPtr[kminus1][x_new(iplus1, jplus1, length_ext)];
				TuSW = inverse * extendedCoefPtr[kminus1][x_new(iminus1, jplus1, length_ext)];
				Bu = inverse * extendedCoefPtr[kplus1][x_ext];
				BuN = inverse * extendedCoefPtr[kplus1][x_new(i_ext, jminus1, length_ext)];
				BuS = inverse * extendedCoefPtr[kplus1][x_new(i_ext, jplus1, length_ext)];
				BuE = inverse * extendedCoefPtr[kplus1][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				BuW = inverse * extendedCoefPtr[kplus1][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				BuNW = inverse * extendedCoefPtr[kplus1][x_new(iminus1, jminus1, length_ext)];
				BuNE = inverse * extendedCoefPtr[kplus1][x_new(iplus1, jminus1, length_ext)];
				BuSE = inverse * extendedCoefPtr[kplus1][x_new(iplus1, jplus1, length_ext)];
				BuSW = inverse * extendedCoefPtr[kplus1][x_new(iminus1, jplus1, length_ext)];

				//calculation of coefficients in the presmooted image data

				// Calculation of coefficients in east direction
				ux = (uE - u) / segParameters.h;
				uy = ((uN + uNE) - (uS + uSE))
					/ quotient;
				uz = ((Tu + TuE) - (Bu + BuE))
					/ quotient;
				GPtrs.GePtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);

				// Calculation of coefficients in west direction
				ux = (uW - u) / segParameters.h;
				uy = ((uNW + uN) - (uSW + uS))
					/ quotient;
				uz = ((TuW + Tu) - (BuW + Bu))
					/ quotient;
				GPtrs.GwPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);

				// Calculation of coefficients in north direction
				ux = ((uNE + uE) - (uNW + uW))
					/ quotient;
				uy = (uN - u) / segParameters.h;
				uz = ((TuN + Tu) - (BuN + Bu))
					/ quotient;
				GPtrs.GnPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);

				// Calculation of coefficients in south direction
				ux = ((uE + uSE) - (uW + uSW))
					/ quotient;
				uy = (uS - u) / segParameters.h;
				uz = ((TuS + Tu) - (BuS + Bu))
					/ quotient;
				GPtrs.GsPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);

				// Calculation of coefficients in top direction
				ux = ((TuE + uE) - (TuW + uW))
					/ quotient;
				uy = ((TuN + uN) - (TuS + uS))
					/ quotient;
				uz = (Tu - u) / segParameters.h;
				GPtrs.GtPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);

				// Calculation of coefficients in bottom direction
				ux = ((BuW + uW) - (BuE + uE))
					/ quotient;
				uy = ((BuN + uN) - (BuS + uS))
					/ quotient;
				uz = (Bu - u) / segParameters.h;
				GPtrs.GbPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);
			}
		}
	}
	return true;
}

bool gaussSeidelCoefficients(dataType **extendedCoefPtr, Segment_Image_Data inputImageData, Gradient_Pointers GPtrs, Coefficient_Pointers CoefPtrs, Segmentation_Parameters segParameters)
{
	//checks if the memory was allocated
	if (extendedCoefPtr == NULL || inputImageData.segmentationFuntionPtr == NULL || GPtrs.GePtr == NULL || GPtrs.GwPtr == NULL
		|| GPtrs.GnPtr == NULL || GPtrs.GsPtr == NULL || GPtrs.GtPtr == NULL || GPtrs.GbPtr == NULL || CoefPtrs.e_Ptr == NULL
		|| CoefPtrs.w_Ptr == NULL || CoefPtrs.n_Ptr == NULL || CoefPtrs.s_Ptr == NULL || CoefPtrs.t_Ptr == NULL || CoefPtrs.b_Ptr == NULL)
		return false;

	size_t i, j, k, x, x_ext; // length == xDim, width == yDim, height == zDim
	size_t kplus1, kminus1, iminus1, iplus1, jminus1, jplus1;
	size_t dim2D = inputImageData.length * inputImageData.width;
	size_t k_ext, j_ext, i_ext;
	size_t height_ext = inputImageData.height + 2;
	size_t length_ext = inputImageData.length + 2;
	size_t width_ext = inputImageData.width + 2;
	dataType quotient = 4.0 * segParameters.h;
	dataType orig_ux, orig_uy, orig_uz; //change in x, y and z respectively
	dataType orig_u, orig_uN, orig_uS, orig_uE, orig_uW, orig_uNW, orig_uNE, orig_uSE, orig_uSW, orig_Tu, orig_TuN, orig_TuS,
		orig_TuE, orig_TuW, orig_TuNW, orig_TuNE, orig_TuSE, orig_TuSW, //current and surrounding voxel values
		orig_Bu, orig_BuN, orig_BuS, orig_BuE, orig_BuW, orig_BuNW, orig_BuNE, orig_BuSE, orig_BuSW;
	dataType orig_e, orig_w, orig_n, orig_s, orig_t, orig_b;
	dataType voxel_coef, average_face_coef;

	//copy data to extended area which will be used in each time step
	copyDataToExtendedArea(inputImageData.segmentationFuntionPtr, extendedCoefPtr, inputImageData.height, inputImageData.length, inputImageData.width);

	//set boundary values to ensure Zero Dirichlet boundary condition.
	setBoundaryToZeroDirichletBC(extendedCoefPtr, length_ext, width_ext, height_ext);

	//calculation of coefficients
	for (k = 0, k_ext = 1; k < inputImageData.height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < inputImageData.length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < inputImageData.width; j++, j_ext++)
			{
				// 2D to 1D representation for i, j
				x_ext = x_new(i_ext, j_ext, length_ext);
				x = x_new(i, j, inputImageData.length);
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
				orig_ux = (orig_uE - orig_u) / segParameters.h;
				orig_uy = ((orig_uN + orig_uNE) - (orig_uS + orig_uSE))
					/ quotient;
				orig_uz = ((orig_Tu + orig_TuE) - (orig_Bu + orig_BuE))
					/ quotient;
				orig_e = sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in west direction
				orig_ux = (orig_uW - orig_u) / segParameters.h;
				orig_uy = ((orig_uNW + orig_uN) - (orig_uSW + orig_uS))
					/ quotient;
				orig_uz = ((orig_TuW + orig_Tu) - (orig_BuW + orig_Bu))
					/ quotient;
				orig_w = sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in north direction
				orig_ux = ((orig_uNE + orig_uE) - (orig_uNW + orig_uW))
					/ quotient;
				orig_uy = (orig_uN - orig_u) / segParameters.h;
				orig_uz = ((orig_TuN + orig_Tu) - (orig_BuN + orig_Bu))
					/ quotient;
				orig_n = sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in south direction
				orig_ux = ((orig_uE + orig_uSE) - (orig_uW + orig_uSW))
					/ quotient;
				orig_uy = (orig_uS - orig_u) / segParameters.h;
				orig_uz = ((orig_TuS + orig_Tu) - (orig_BuS + orig_Bu))
					/ quotient;
				orig_s = sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in top direction
				orig_ux = ((orig_TuE + orig_uE) - (orig_TuW + orig_uW))
					/ quotient;
				orig_uy = ((orig_TuN + orig_uN) - (orig_TuS + orig_uS))
					/ quotient;
				orig_uz = (orig_Tu - orig_u) / segParameters.h;
				orig_t = sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in bottom direction
				orig_ux = ((orig_BuW + orig_uW) - (orig_BuE + orig_uE))
					/ quotient;
				orig_uy = ((orig_BuN + orig_uN) - (orig_BuS + orig_uS))
					/ quotient;
				orig_uz = (orig_Bu - orig_u) / segParameters.h;
				orig_b = sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// evaluation of norm of gradient of image at each voxel
				average_face_coef = ((orig_e + orig_w + orig_n + orig_s + orig_t + orig_b) / 6.0);

				voxel_coef = sqrt(pow(average_face_coef, 2) + segParameters.eps2);

				/* evaluation of norm of gradient of image at each voxel, norm of gradient of presmoothed
				image at each voxel face and reciprocal of norm of gradient of image at each voxel face*/
				CoefPtrs.e_Ptr[k][x] = voxel_coef * GPtrs.GePtr[k][x] * (1.0 / orig_e);//east coefficient
				CoefPtrs.w_Ptr[k][x] = voxel_coef * GPtrs.GwPtr[k][x] * (1.0 / orig_w);//west coefficient
				CoefPtrs.n_Ptr[k][x] = voxel_coef * GPtrs.GnPtr[k][x] * (1.0 / orig_n);//north coefficient
				CoefPtrs.s_Ptr[k][x] = voxel_coef * GPtrs.GsPtr[k][x] * (1.0 / orig_s);//south coefficient
				CoefPtrs.t_Ptr[k][x] = voxel_coef * GPtrs.GtPtr[k][x] * (1.0 / orig_t);//top coefficient
				CoefPtrs.b_Ptr[k][x] = voxel_coef * GPtrs.GbPtr[k][x] * (1.0 / orig_b);//bottom coefficient
			}
		}
	}

	return true;
}