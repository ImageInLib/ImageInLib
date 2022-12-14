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
#include <common_vtk.h>
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
	size_t i, j, k; // length == xDim, width == yDim, height == zDim
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
	generateInitialSegmentationFunctionForMultipleCentres(segmFuntionPtr, length, width, height, centers, 0.5, 15, no_of_centers);
	//copyDataToAnotherArray(initialSegment, segmFuntionPtr, height, length, width);

	//compute coefficients from presmoothed image
	gFunctionForImageToBeSegmented(inputImageData, prevSol_extPtr, GPtrs, segParameters, explicit_lhe_Parameters);

	Vtk_File_Info* vtkInfo = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	if (vtkInfo == NULL) return false;
	vtkInfo->spacing[0] = 1.0; vtkInfo->spacing[1] = 1.0; vtkInfo->spacing[2] = 1.0;
	vtkInfo->origin[0] = 0; vtkInfo->origin[1] = 0; vtkInfo->origin[2] = 0;
	vtkInfo->dimensions[0] = length; vtkInfo->dimensions[1] = width; vtkInfo->dimensions[2] = height;
	vtkInfo->vDataType = dta_Flt; vtkInfo->operation = copyTo; vtkDataForm dataForm = dta_binary;
	const char* pathsaveVTK;

	//Array for name construction
	unsigned char name[350];
	unsigned char name_ending[100];


	strcpy_s(name, sizeof name, outputPathPtr);
	sprintf_s(name_ending, sizeof(name_ending), "_g_est%03zd.vtk", i);
	strcat_s(name, sizeof(name), name_ending);
	pathsaveVTK = name;
	vtkInfo->dataPointer = GPtrs.GePtr;
	storeVtkFile(pathsaveVTK, vtkInfo, dataForm);

	//strcpy_s(name, sizeof name, outputPathPtr);
	//sprintf_s(name_ending, sizeof(name_ending), "_g_west.vtk", i);
	//strcat_s(name, sizeof(name), name_ending);
	//pathsaveVTK = name;
	//vtkInfo->dataPointer = GPtrs.GwPtr;
	//storeVtkFile(pathsaveVTK, vtkInfo, dataForm);

	//strcpy_s(name, sizeof name, outputPathPtr);
	//sprintf_s(name_ending, sizeof(name_ending), "_g_north.vtk", i);
	//strcat_s(name, sizeof(name), name_ending);
	//pathsaveVTK = name;
	//vtkInfo->dataPointer = GPtrs.GnPtr;
	//storeVtkFile(pathsaveVTK, vtkInfo, dataForm);

	//strcpy_s(name, sizeof name, outputPathPtr);
	//sprintf_s(name_ending, sizeof(name_ending), "_g_south.vtk", i);
	//strcat_s(name, sizeof(name), name_ending);
	//pathsaveVTK = name;
	//vtkInfo->dataPointer = GPtrs.GsPtr;
	//storeVtkFile(pathsaveVTK, vtkInfo, dataForm);

	//strcpy_s(name, sizeof name, outputPathPtr);
	//sprintf_s(name_ending, sizeof(name_ending), "_g_top.vtk", i);
	//strcat_s(name, sizeof(name), name_ending);
	//pathsaveVTK = name;
	//vtkInfo->dataPointer = GPtrs.GtPtr;
	//storeVtkFile(pathsaveVTK, vtkInfo, dataForm);

	//strcpy_s(name, sizeof name, outputPathPtr);
	//sprintf_s(name_ending, sizeof(name_ending), "_g_bottom.vtk", i);
	//strcat_s(name, sizeof(name), name_ending);
	//pathsaveVTK = name;
	//vtkInfo->dataPointer = GPtrs.GbPtr;
	//storeVtkFile(pathsaveVTK, vtkInfo, dataForm);


	//vtkInfo->dataPointer = inputDataArrayPtr;

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
			pathsaveVTK = name;
			vtkInfo->dataPointer = imageData.segmentationFuntionPtr;
			storeVtkFile(pathsaveVTK, vtkInfo, dataForm);
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

	free(vtkInfo);

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

					mean_square_residue += (dataType)(pow(gauss_seidelPtr[k_ext][x_ext] * (1 + coef_tauh * (CoefPtrs.e_Ptr[k][x]
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
	quotient = (dataType)1. / (max - min);
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
				norm_of_distance = (dataType)sqrt((dx * dx) + (dy * dy) + (dz * dz));
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
	quotient = (dataType)(1. / (max - min));
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
				norm_of_distance = (dataType)sqrt((dx * dx) + (dy * dy) + (dz * dz));
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

	Vtk_File_Info* vtkInfo = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	vtkInfo->spacing[0] = 1.0; vtkInfo->spacing[1] = 1.0; vtkInfo->spacing[2] = 1.0;
	vtkInfo->origin[0] = 0; vtkInfo->origin[1] = 0; vtkInfo->origin[2] = 0;
	vtkInfo->dimensions[0] = length; vtkInfo->dimensions[1] = width; vtkInfo->dimensions[2] = height;
	vtkInfo->vDataType = dta_Flt; vtkInfo->dataPointer = inputDataArrayPtr; vtkInfo->operation = copyTo;
	vtkDataForm dataForm = dta_binary;
	const char* pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/segFunction.vtk";
	
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
					norm_of_distance = (dataType)sqrt((dx * dx) + (dy * dy) + (dz * dz));
					new_value = (dataType)((1.0 / (sqrt((dx * dx) + (dy * dy) + (dz * dz)) + v)) - (1. / (R + v)));
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
	//Storage_Flags flags = { true, true };
	//store3dDataVtkD(inputDataArrayPtr, length, width, height, pathArray1, (2.5 / (length)), flags);
	storeVtkFile(pathsaveVTK, vtkInfo, dataForm);
	free(vtkInfo);
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
	dataType quotient = (dataType)(4.0 * segParameters.h);
	dataType inverse = 1;//(dataType)(1.0 / (VTK_MAX_HEADER_LINE_LENGTH - 1));
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
	//heatExplicitScheme(presmoothingData, explicit_lhe_Parameters);
	heatImplicitScheme(presmoothingData, explicit_lhe_Parameters);
	//geodesicMeanCurvatureTimeStep(presmoothingData, explicit_lhe_Parameters);
	//meanCurvatureTimeStep(presmoothingData, explicit_lhe_Parameters);

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
	dataType quotient = (dataType)(4.0 * segParameters.h);
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
				orig_e = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in west direction
				orig_ux = (orig_uW - orig_u) / segParameters.h;
				orig_uy = ((orig_uNW + orig_uN) - (orig_uSW + orig_uS))
					/ quotient;
				orig_uz = ((orig_TuW + orig_Tu) - (orig_BuW + orig_Bu))
					/ quotient;
				orig_w = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in north direction
				orig_ux = ((orig_uNE + orig_uE) - (orig_uNW + orig_uW))
					/ quotient;
				orig_uy = (orig_uN - orig_u) / segParameters.h;
				orig_uz = ((orig_TuN + orig_Tu) - (orig_BuN + orig_Bu))
					/ quotient;
				orig_n = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in south direction
				orig_ux = ((orig_uE + orig_uSE) - (orig_uW + orig_uSW))
					/ quotient;
				orig_uy = (orig_uS - orig_u) / segParameters.h;
				orig_uz = ((orig_TuS + orig_Tu) - (orig_BuS + orig_Bu))
					/ quotient;
				orig_s = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in top direction
				orig_ux = ((orig_TuE + orig_uE) - (orig_TuW + orig_uW))
					/ quotient;
				orig_uy = ((orig_TuN + orig_uN) - (orig_TuS + orig_uS))
					/ quotient;
				orig_uz = (orig_Tu - orig_u) / segParameters.h;
				orig_t = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in bottom direction
				orig_ux = ((orig_BuW + orig_uW) - (orig_BuE + orig_uE))
					/ quotient;
				orig_uy = ((orig_BuN + orig_uN) - (orig_BuS + orig_uS))
					/ quotient;
				orig_uz = (orig_Bu - orig_u) / segParameters.h;
				orig_b = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// evaluation of norm of gradient of image at each voxel
				average_face_coef = (dataType)(((orig_e + orig_w + orig_n + orig_s + orig_t + orig_b) / 6.0));

				voxel_coef = (dataType)sqrt(pow(average_face_coef, 2) + segParameters.eps2);

				/* evaluation of norm of gradient of image at each voxel, norm of gradient of presmoothed
				image at each voxel face and reciprocal of norm of gradient of image at each voxel face*/
				CoefPtrs.e_Ptr[k][x] = (dataType)(voxel_coef * GPtrs.GePtr[k][x] * (1.0 / orig_e));//east coefficient
				CoefPtrs.w_Ptr[k][x] = (dataType)(voxel_coef * GPtrs.GwPtr[k][x] * (1.0 / orig_w));//west coefficient
				CoefPtrs.n_Ptr[k][x] = (dataType)(voxel_coef * GPtrs.GnPtr[k][x] * (1.0 / orig_n));//north coefficient
				CoefPtrs.s_Ptr[k][x] = (dataType)(voxel_coef * GPtrs.GsPtr[k][x] * (1.0 / orig_s));//south coefficient
				CoefPtrs.t_Ptr[k][x] = (dataType)(voxel_coef * GPtrs.GtPtr[k][x] * (1.0 / orig_t));//top coefficient
				CoefPtrs.b_Ptr[k][x] = (dataType)(voxel_coef * GPtrs.GbPtr[k][x] * (1.0 / orig_b));//bottom coefficient
			}
		}
	}

	return true;
}


//bool gSubsurfSegmentation(Image_Data inputImageData, dataType** initialSegment, Segmentation_Parameters segParameters, Filter_Parameters explicit_lhe_Parameters, unsigned char* outputPathPtr, dataType alpha, dataType beta) {
//	
//	if (inputImageData.imageDataPtr == NULL || initialSegment == NULL)
//		return false;
//
//	size_t i, j, k, x, xd, x_ext;
//	size_t kplus1, kminus1, iminus1, iplus1, jminus1, jplus1;
//	size_t k_ext, j_ext, i_ext;
//	size_t dim2D = inputImageData.length * inputImageData.width;
//	size_t height = inputImageData.height;
//	size_t length = inputImageData.length;
//	size_t width = inputImageData.width;
//	size_t height_ext = inputImageData.height + 2;
//	size_t length_ext = inputImageData.length + 2;
//	size_t width_ext = inputImageData.width + 2;
//	size_t dim2D_ext = length_ext * width_ext;
//	dataType quotient = (dataType)(4.0 * segParameters.h);
//	dataType hhh = segParameters.h * segParameters.h * segParameters.h, coef_tauh = segParameters.tau / hhh;
//	dataType ux, uy, uz; 
//	dataType u, uN, uS, uE, uW, uNW, uNE, uSE, uSW, Tu, TuN, TuS, TuE, TuW, TuNW, TuNE, TuSE, TuSW, //current and surrounding voxel values
//		Bu, BuN, BuS, BuE, BuW, BuNW, BuNE, BuSE, BuSW;
//	dataType gradient_average;
//	dataType orig_ux, orig_uy, orig_uz;
//	dataType orig_u, orig_uN, orig_uS, orig_uE, orig_uW, orig_uNW, orig_uNE, orig_uSE, orig_uSW, orig_Tu, orig_TuN, orig_TuS,
//		orig_TuE, orig_TuW, orig_TuNW, orig_TuNE, orig_TuSE, orig_TuSW, //current and surrounding voxel values
//		orig_Bu, orig_BuN, orig_BuS, orig_BuE, orig_BuW, orig_BuNW, orig_BuNE, orig_BuSE, orig_BuSW;
//	dataType orig_e, orig_w, orig_n, orig_s, orig_t, orig_b;
//	dataType voxel_coef, average_face_coef;
//	dataType firstCpuTime, secondCpuTime, difference_btw_current_and_previous_sol, gauss_seidel;
//	
//	dataType** segmFuntionPtr = (dataType**)malloc(sizeof(dataType*) * height);
//
//	dataType** GePtr = (dataType**)malloc(sizeof(dataType*) * height);
//	dataType** GwPtr = (dataType**)malloc(sizeof(dataType*) * height);
//	dataType** GnPtr = (dataType**)malloc(sizeof(dataType*) * height);
//	dataType** GsPtr = (dataType**)malloc(sizeof(dataType*) * height);
//	dataType** GtPtr = (dataType**)malloc(sizeof(dataType*) * height);
//	dataType** GbPtr = (dataType**)malloc(sizeof(dataType*) * height);
//
//	dataType** e_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
//	dataType** w_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
//	dataType** n_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
//	dataType** s_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
//	dataType** t_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
//	dataType** b_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
//
//	dataType** VePtr = (dataType**)malloc(sizeof(dataType*) * height);
//	dataType** VwPtr = (dataType**)malloc(sizeof(dataType*) * height);
//	dataType** VnPtr = (dataType**)malloc(sizeof(dataType*) * height);
//	dataType** VsPtr = (dataType**)malloc(sizeof(dataType*) * height);
//	dataType** VtPtr = (dataType**)malloc(sizeof(dataType*) * height);
//	dataType** VbPtr = (dataType**)malloc(sizeof(dataType*) * height);
//
//	dataType** edgeGradient = (dataType**)malloc(sizeof(dataType*) * height);
//
//	for (k = 0; k < height; k++)
//	{
//		segmFuntionPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
//
//		GePtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
//		GwPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
//		GnPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
//		GsPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
//		GtPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
//		GbPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
//		e_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
//		w_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
//		n_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
//		s_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
//		t_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
//		b_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
//
//		VePtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
//		VwPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
//		VnPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
//		VsPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
//		VtPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
//		VbPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
//
//		edgeGradient[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
//	}
//	//checks if the memory was allocated
//	if (segmFuntionPtr == NULL ||
//		GePtr == NULL || GwPtr == NULL || GnPtr == NULL || GsPtr == NULL || GtPtr == NULL || GbPtr == NULL || 
//		e_Ptr == NULL || w_Ptr == NULL || n_Ptr == NULL || s_Ptr == NULL || t_Ptr == NULL || b_Ptr == NULL || 
//		VePtr == NULL || VwPtr == NULL || VnPtr == NULL || VsPtr == NULL || VtPtr == NULL || VnPtr == NULL ||
//		edgeGradient == NULL)
//		return false;
//
//	dataType** gauss_seidelPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
//	dataType** prevSol_extPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
//	dataType** extendedCoefPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
//	for (k = 0; k < height_ext; k++){
//		gauss_seidelPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
//		prevSol_extPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
//		extendedCoefPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
//	}
//	//checks if the memory was allocated
//	if (gauss_seidelPtr == NULL || prevSol_extPtr == NULL || extendedCoefPtr == NULL)
//		return false;
//
//	copyDataToExtendedArea(inputImageData.imageDataPtr, extendedCoefPtr, height, length, width);
//
//	//perform reflection of the extended area to ensure zero Neumann boundary condition(for LHE)
//	reflection3D(extendedCoefPtr, height_ext, length_ext, width_ext);
//
//	Image_Data presmoothingData;
//	presmoothingData.height = height_ext; presmoothingData.length = length_ext; presmoothingData.width = width_ext;
//	presmoothingData.imageDataPtr = extendedCoefPtr;
//	heatImplicitScheme(presmoothingData, explicit_lhe_Parameters);
//
//	//compute edges detector function and gradient
//	dataType norm_gradient_smoothed_e, norm_gradient_smoothed_w, norm_gradient_smoothed_n, norm_gradient_smoothed_s, norm_gradient_smoothed_b, norm_gradient_smoothed_t;
//	for (k = 0, k_ext = 1; k < height; k++, k_ext++){
//		for (i = 0, i_ext = 1; i < length; i++, i_ext++){
//			for (j = 0, j_ext = 1; j < width; j++, j_ext++){
//
//				// 2D to 1D representation for i, j
//				x_ext = x_new(i_ext, j_ext, length_ext);
//				x = x_new(i, j, length);
//				iminus1 = i_ext - 1;
//				iplus1 = i_ext + 1;
//				jplus1 = j_ext + 1;
//				jminus1 = j_ext - 1;
//				kplus1 = k_ext + 1;
//				kminus1 = k_ext - 1;
//
//				//values of voxels in the extended data container for presmoothed image
//				u = extendedCoefPtr[k_ext][x_ext];
//				uN = extendedCoefPtr[k_ext][x_new(i_ext, jminus1, length_ext)];
//				uS = extendedCoefPtr[k_ext][x_new(i_ext, jplus1, length_ext)];
//				uE = extendedCoefPtr[k_ext][x_new(iplus1, j_ext, length_ext)];
//				uW = extendedCoefPtr[k_ext][x_new(iminus1, j_ext, length_ext)];
//				uNW = extendedCoefPtr[k_ext][x_new(iminus1, jminus1, length_ext)];
//				uNE = extendedCoefPtr[k_ext][x_new(iplus1, jminus1, length_ext)];
//				uSE = extendedCoefPtr[k_ext][x_new(iplus1, jplus1, length_ext)];
//				uSW = extendedCoefPtr[k_ext][x_new(iminus1, jplus1, length_ext)];
//				Tu = extendedCoefPtr[kminus1][x_ext];
//				TuN = extendedCoefPtr[kminus1][x_new(i_ext, jminus1, length_ext)];
//				TuS = extendedCoefPtr[kminus1][x_new(i_ext, jplus1, length_ext)];
//				TuE = extendedCoefPtr[kminus1][x_new(iplus1, j_ext, length_ext)];
//				TuW = extendedCoefPtr[kminus1][x_new(iminus1, j_ext, length_ext)];
//				TuNW = extendedCoefPtr[kminus1][x_new(iminus1, jminus1, length_ext)];
//				TuNE = extendedCoefPtr[kminus1][x_new(iplus1, jminus1, length_ext)];
//				TuSE = extendedCoefPtr[kminus1][x_new(iplus1, jplus1, length_ext)];
//				TuSW = extendedCoefPtr[kminus1][x_new(iminus1, jplus1, length_ext)];
//				Bu = extendedCoefPtr[kplus1][x_ext];
//				BuN = extendedCoefPtr[kplus1][x_new(i_ext, jminus1, length_ext)];
//				BuS = extendedCoefPtr[kplus1][x_new(i_ext, jplus1, length_ext)];
//				BuE = extendedCoefPtr[kplus1][x_new(iplus1, j_ext, length_ext)];
//				BuW = extendedCoefPtr[kplus1][x_new(iminus1, j_ext, length_ext)];
//				BuNW = extendedCoefPtr[kplus1][x_new(iminus1, jminus1, length_ext)];
//				BuNE = extendedCoefPtr[kplus1][x_new(iplus1, jminus1, length_ext)];
//				BuSE = extendedCoefPtr[kplus1][x_new(iplus1, jplus1, length_ext)];
//				BuSW = extendedCoefPtr[kplus1][x_new(iminus1, jplus1, length_ext)];
//
//				//calculation of coefficients in the presmooted image data
//
//				// Calculation of coefficients in east direction
//				ux = (uE - u) / segParameters.h;
//				uy = ((uN + uNE) - (uS + uSE))
//					/ quotient;
//				uz = ((Tu + TuE) - (Bu + BuE))
//					/ quotient;
//				GePtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);
//
//				norm_gradient_smoothed_e = sqrt( (ux * ux) + (uy * uy) + (uz * uz) );
//
//				// Calculation of coefficients in west direction
//				ux = (uW - u) / segParameters.h;
//				uy = ((uNW + uN) - (uSW + uS))
//					/ quotient;
//				uz = ((TuW + Tu) - (BuW + Bu))
//					/ quotient;
//				GwPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);
//
//				norm_gradient_smoothed_w = sqrt( (ux * ux) + (uy * uy) + (uz * uz) );
//
//				// Calculation of coefficients in north direction
//				ux = ((uNE + uE) - (uNW + uW))
//					/ quotient;
//				uy = (uN - u) / segParameters.h;
//				uz = ((TuN + Tu) - (BuN + Bu))
//					/ quotient;
//				GnPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);
//
//				norm_gradient_smoothed_n = sqrt((ux * ux) + (uy * uy) + (uz * uz));
//
//				// Calculation of coefficients in south direction
//				ux = ((uE + uSE) - (uW + uSW))
//					/ quotient;
//				uy = (uS - u) / segParameters.h;
//				uz = ((TuS + Tu) - (BuS + Bu))
//					/ quotient;
//				GsPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);
//
//				norm_gradient_smoothed_s = sqrt((ux * ux) + (uy * uy) + (uz * uz));
//
//				// Calculation of coefficients in top direction
//				ux = ((TuE + uE) - (TuW + uW))
//					/ quotient;
//				uy = ((TuN + uN) - (TuS + uS))
//					/ quotient;
//				uz = (Tu - u) / segParameters.h;
//				GtPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);
//
//				norm_gradient_smoothed_t = sqrt((ux * ux) + (uy * uy) + (uz * uz) );
//
//				// Calculation of coefficients in bottom direction
//				ux = ((BuW + uW) - (BuE + uE))
//					/ quotient;
//				uy = ((BuN + uN) - (BuS + uS))
//					/ quotient;
//				uz = (Bu - u) / segParameters.h;
//				GbPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);
//
//				norm_gradient_smoothed_b = sqrt((ux * ux) + (uy * uy) + (uz * uz) );
//
//				//gradient_average = (dataType)((GePtr[k][x] + GwPtr[k][x] + GnPtr[k][x] + GsPtr[k][x] + GtPtr[k][x] + GbPtr[k][x]) / 6.0);
//				gradient_average = (norm_gradient_smoothed_e + norm_gradient_smoothed_w + norm_gradient_smoothed_n +
//					norm_gradient_smoothed_s + norm_gradient_smoothed_t + norm_gradient_smoothed_b) / 6.0;
//
//				edgeGradient[k][x] = gradientFunction(gradient_average, segParameters.coef);
//
//				VePtr[k][x] = alpha * (GePtr[k][x] - edgeGradient[k][x]) / (2 * segParameters.h); 
//				VwPtr[k][x] = alpha * (GwPtr[k][x] - edgeGradient[k][x]) / (2 * segParameters.h); 
//				VnPtr[k][x] = alpha * (GnPtr[k][x] - edgeGradient[k][x]) / (2 * segParameters.h); 
//				VsPtr[k][x] = alpha * (GsPtr[k][x] - edgeGradient[k][x]) / (2 * segParameters.h);
//				VtPtr[k][x] = alpha * (GtPtr[k][x] - edgeGradient[k][x]) / (2 * segParameters.h);
//				VbPtr[k][x] = alpha * (GbPtr[k][x] - edgeGradient[k][x]) / (2 * segParameters.h);
//			}
//		}
//	}
//
//	//Preparation of saving in vtk
//	Vtk_File_Info* vtkInfo = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
//	if (vtkInfo == NULL) return false;
//	vtkInfo->spacing[0] = 1.0; vtkInfo->spacing[1] = 1.0; vtkInfo->spacing[2] = 1.0;
//	vtkInfo->origin[0] = 0; vtkInfo->origin[1] = 0; vtkInfo->origin[2] = 0;
//	vtkInfo->dimensions[0] = length; vtkInfo->dimensions[1] = width; vtkInfo->dimensions[2] = height;
//	vtkInfo->vDataType = dta_Flt; vtkInfo->operation = copyTo; vtkDataForm dataForm = dta_binary;
//	const char* pathsaveVTK;
//
//	//Array for name construction
//	unsigned char name[350];
//	unsigned char name_ending[100];
//
//	
//	strcpy_s(name, sizeof name, outputPathPtr);
//	sprintf_s(name_ending, sizeof(name_ending), "_Gest_%03zd.vtk", i);
//	strcat_s(name, sizeof(name), name_ending);
//	pathsaveVTK = name;
//	vtkInfo->dataPointer = GePtr;
//	storeVtkFile(pathsaveVTK, vtkInfo, dataForm);
//
//	//strcpy_s(name, sizeof name, outputPathPtr);
//	//sprintf_s(name_ending, sizeof(name_ending), "_Gwest_%03zd.vtk", i);
//	//strcat_s(name, sizeof(name), name_ending);
//	//pathsaveVTK = name;
//	//vtkInfo->dataPointer = GwPtr;
//	//storeVtkFile(pathsaveVTK, vtkInfo, dataForm);
//
//	//strcpy_s(name, sizeof name, outputPathPtr);
//	//sprintf_s(name_ending, sizeof(name_ending), "_Gnorth_%03zd.vtk", i);
//	//strcat_s(name, sizeof(name), name_ending);
//	//pathsaveVTK = name;
//	//vtkInfo->dataPointer = GnPtr;
//	//storeVtkFile(pathsaveVTK, vtkInfo, dataForm);
//
//	//strcpy_s(name, sizeof name, outputPathPtr);
//	//sprintf_s(name_ending, sizeof(name_ending), "_Gsouth_%03zd.vtk", i);
//	//strcat_s(name, sizeof(name), name_ending);
//	//pathsaveVTK = name;
//	//vtkInfo->dataPointer = GsPtr;
//	//storeVtkFile(pathsaveVTK, vtkInfo, dataForm);
//
//	//strcpy_s(name, sizeof name, outputPathPtr);
//	//sprintf_s(name_ending, sizeof(name_ending), "_Gtop_%03zd.vtk", i);
//	//strcat_s(name, sizeof(name), name_ending);
//	//pathsaveVTK = name;
//	//vtkInfo->dataPointer = GtPtr;
//	//storeVtkFile(pathsaveVTK, vtkInfo, dataForm);
//
//	//strcpy_s(name, sizeof name, outputPathPtr);
//	//sprintf_s(name_ending, sizeof(name_ending), "_Gbottom_%03zd.vtk", i);
//	//strcat_s(name, sizeof(name), name_ending);
//	//pathsaveVTK = name;
//	//vtkInfo->dataPointer = GbPtr;
//	//storeVtkFile(pathsaveVTK, vtkInfo, dataForm);
//
//	//generate initial segmentation function
//	//generateInitialSegmentationFunctionForMultipleCentres(segmFuntionPtr, length, width, height, centers, 0.5, 5, no_of_centers);
//	copyDataToAnotherArray(initialSegment, segmFuntionPtr, height, length, width);
//
//	copyDataToExtendedArea(segmFuntionPtr, gauss_seidelPtr, height, length, width);
//	setBoundaryToZeroDirichletBC(gauss_seidelPtr, length_ext, width_ext, height_ext);
//
//	//loop for segmentation
//	size_t cpt = 1, z;
//	dataType mean_square_residue;
//	do
//	{
//		segParameters.numberOfTimeStep = cpt;
//		firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
//
//		//copy data to extended area which will be used in each time step
//		copyDataToExtendedArea(segmFuntionPtr, extendedCoefPtr, height, length, width);
//		setBoundaryToZeroDirichletBC(extendedCoefPtr, length_ext, width_ext, height_ext);
//		copyDataToAnotherArray(extendedCoefPtr, prevSol_extPtr, height_ext, length_ext, width_ext);
//
//		//compute coefficients
//		for (k = 0, k_ext = 1; k < height; k++, k_ext++){
//			for (i = 0, i_ext = 1; i < length; i++, i_ext++){
//				for (j = 0, j_ext = 1; j < width; j++, j_ext++){
//
//					// 2D to 1D representation for i, j
//					x_ext = x_new(i_ext, j_ext, length_ext);
//					x = x_new(i, j, length);
//					iminus1 = i_ext - 1;
//					iplus1 = i_ext + 1;
//					jplus1 = j_ext + 1;
//					jminus1 = j_ext - 1;
//					kplus1 = k_ext + 1;
//					kminus1 = k_ext - 1;
//
//					//values of voxels in the extended data container for the original image
//					orig_u = extendedCoefPtr[k_ext][x_ext];
//					orig_uN = extendedCoefPtr[k_ext][x_new(i_ext, jminus1, length_ext)];
//					orig_uS = extendedCoefPtr[k_ext][x_new(i_ext, jplus1, length_ext)];
//					orig_uE = extendedCoefPtr[k_ext][x_new(iplus1, j_ext, length_ext)];
//					orig_uW = extendedCoefPtr[k_ext][x_new(iminus1, j_ext, length_ext)];
//					orig_uNW = extendedCoefPtr[k_ext][x_new(iminus1, jminus1, length_ext)];
//					orig_uNE = extendedCoefPtr[k_ext][x_new(iplus1, jminus1, length_ext)];
//					orig_uSE = extendedCoefPtr[k_ext][x_new(iplus1, jplus1, length_ext)];
//					orig_uSW = extendedCoefPtr[k_ext][x_new(iminus1, jplus1, length_ext)];
//					orig_Tu = extendedCoefPtr[kminus1][x_ext];
//					orig_TuN = extendedCoefPtr[kminus1][x_new(i_ext, jminus1, length_ext)];
//					orig_TuS = extendedCoefPtr[kminus1][x_new(i_ext, jplus1, length_ext)];
//					orig_TuE = extendedCoefPtr[kminus1][x_new(iplus1, j_ext, length_ext)];
//					orig_TuW = extendedCoefPtr[kminus1][x_new(iminus1, j_ext, length_ext)];
//					orig_TuNW = extendedCoefPtr[kminus1][x_new(iminus1, jminus1, length_ext)];
//					orig_TuNE = extendedCoefPtr[kminus1][x_new(iplus1, jminus1, length_ext)];
//					orig_TuSE = extendedCoefPtr[kminus1][x_new(iplus1, jplus1, length_ext)];
//					orig_TuSW = extendedCoefPtr[kminus1][x_new(iminus1, jplus1, length_ext)];
//					orig_Bu = extendedCoefPtr[kplus1][x_ext];
//					orig_BuN = extendedCoefPtr[kplus1][x_new(i_ext, jminus1, length_ext)];
//					orig_BuS = extendedCoefPtr[kplus1][x_new(i_ext, jplus1, length_ext)];
//					orig_BuE = extendedCoefPtr[kplus1][x_new(iplus1, j_ext, length_ext)];
//					orig_BuW = extendedCoefPtr[kplus1][x_new(iminus1, j_ext, length_ext)];
//					orig_BuNW = extendedCoefPtr[kplus1][x_new(iminus1, jminus1, length_ext)];
//					orig_BuNE = extendedCoefPtr[kplus1][x_new(iplus1, jminus1, length_ext)];
//					orig_BuSE = extendedCoefPtr[kplus1][x_new(iplus1, jplus1, length_ext)];
//					orig_BuSW = extendedCoefPtr[kplus1][x_new(iminus1, jplus1, length_ext)];
//
//					//calculation of coefficients in the original image data
//					// Calculation of coefficients in east direction
//					orig_ux = (orig_uE - orig_u) / segParameters.h;
//					orig_uy = ((orig_uN + orig_uNE) - (orig_uS + orig_uSE))
//						/ quotient;
//					orig_uz = ((orig_Tu + orig_TuE) - (orig_Bu + orig_BuE))
//						/ quotient;
//					orig_e = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);
//
//					// Calculation of coefficients in west direction
//					orig_ux = (orig_uW - orig_u) / segParameters.h;
//					orig_uy = ((orig_uNW + orig_uN) - (orig_uSW + orig_uS))
//						/ quotient;
//					orig_uz = ((orig_TuW + orig_Tu) - (orig_BuW + orig_Bu))
//						/ quotient;
//					orig_w = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);
//
//					// Calculation of coefficients in north direction
//					orig_ux = ((orig_uNE + orig_uE) - (orig_uNW + orig_uW))
//						/ quotient;
//					orig_uy = (orig_uN - orig_u) / segParameters.h;
//					orig_uz = ((orig_TuN + orig_Tu) - (orig_BuN + orig_Bu))
//						/ quotient;
//					orig_n = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);
//
//					// Calculation of coefficients in south direction
//					orig_ux = ((orig_uE + orig_uSE) - (orig_uW + orig_uSW))
//						/ quotient;
//					orig_uy = (orig_uS - orig_u) / segParameters.h;
//					orig_uz = ((orig_TuS + orig_Tu) - (orig_BuS + orig_Bu))
//						/ quotient;
//					orig_s = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);
//
//					// Calculation of coefficients in top direction
//					orig_ux = ((orig_TuE + orig_uE) - (orig_TuW + orig_uW))
//						/ quotient;
//					orig_uy = ((orig_TuN + orig_uN) - (orig_TuS + orig_uS))
//						/ quotient;
//					orig_uz = (orig_Tu - orig_u) / segParameters.h;
//					orig_t = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);
//
//					// Calculation of coefficients in bottom direction
//					orig_ux = ((orig_BuW + orig_uW) - (orig_BuE + orig_uE))
//						/ quotient;
//					orig_uy = ((orig_BuN + orig_uN) - (orig_BuS + orig_uS))
//						/ quotient;
//					orig_uz = (orig_Bu - orig_u) / segParameters.h;
//					orig_b = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);
//
//					// evaluation of norm of gradient of image at each voxel
//					average_face_coef = (dataType)(((orig_e + orig_w + orig_n + orig_s + orig_t + orig_b) / 6.0));
//
//					voxel_coef = (dataType)sqrt(pow(average_face_coef, 2) + segParameters.eps2);
//
//					e_Ptr[k][x] = (dataType)(beta * voxel_coef * GePtr[k][x] * (1.0 / orig_e));//east coefficient
//					w_Ptr[k][x] = (dataType)(beta * voxel_coef * GwPtr[k][x] * (1.0 / orig_w));//west coefficient
//					n_Ptr[k][x] = (dataType)(beta * voxel_coef * GnPtr[k][x] * (1.0 / orig_n));//north coefficient
//					s_Ptr[k][x] = (dataType)(beta * voxel_coef * GsPtr[k][x] * (1.0 / orig_s));//south coefficient
//					t_Ptr[k][x] = (dataType)(beta * voxel_coef * GtPtr[k][x] * (1.0 / orig_t));//top coefficient
//					b_Ptr[k][x] = (dataType)(beta * voxel_coef * GbPtr[k][x] * (1.0 / orig_b));//bottom coefficient
//				}
//			}
//		}
//
//		setBoundaryToZeroDirichletBC(gauss_seidelPtr, length_ext, width_ext, height_ext);
//
//		// The Implicit Scheme Evaluation
//		z = 0;
//		do
//		{
//			z = z + 1;
//			for (k = 0, k_ext = 1; k < height; k++, k_ext++){
//				for (i = 0, i_ext = 1; i < length; i++, i_ext++){
//					for (j = 0, j_ext = 1; j < width; j++, j_ext++){
//
//						// 2D to 1D representation for i, j
//						x_ext = x_new(i_ext, j_ext, length_ext);
//						x = x_new(i, j, length);
//
//						// Begin Gauss-Seidel Formula Evaluation
//						gauss_seidel = ( prevSol_extPtr[k_ext][x_ext] - coef_tauh * ( (prevSol_extPtr[k_ext][x_ext + 1] - prevSol_extPtr[k_ext][x_ext]) * VePtr[k][x]
//							+ (prevSol_extPtr[k_ext][x_ext - 1] - prevSol_extPtr[k_ext][x_ext]) * VwPtr[k][x]
//							+ (prevSol_extPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)] - prevSol_extPtr[k_ext][x_ext]) * VsPtr[k][x]
//							+ (prevSol_extPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)] - prevSol_extPtr[k_ext][x_ext]) * VnPtr[k][x]
//							+ (prevSol_extPtr[k_ext + 1][x_ext] - prevSol_extPtr[k_ext][x_ext]) * VbPtr[k][x]
//							+ (prevSol_extPtr[k_ext - 1][x_ext] - prevSol_extPtr[k_ext][x_ext]) * VtPtr[k][x] )
//							+ coef_tauh * ( e_Ptr[k][x] * gauss_seidelPtr[k_ext][x_ext + 1]
//							+ w_Ptr[k][x] * gauss_seidelPtr[k_ext][x_ext - 1]
//							+ s_Ptr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)]
//							+ n_Ptr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)]
//							+ b_Ptr[k][x] * gauss_seidelPtr[k_ext + 1][x_ext]
//							+ t_Ptr[k][x] * gauss_seidelPtr[k_ext - 1][x_ext]) ) /
//							(1 + coef_tauh * (e_Ptr[k][x] + w_Ptr[k][x] + n_Ptr[k][x]
//								+ s_Ptr[k][x] + t_Ptr[k][x] + b_Ptr[k][x]) );
//
//						// SOR implementation using Gauss-Seidel
//						gauss_seidelPtr[k_ext][x_ext] = gauss_seidelPtr[k_ext][x_ext] +
//							segParameters.omega_c * (gauss_seidel - gauss_seidelPtr[k_ext][x_ext]);
//					}
//				}
//			}
//
//			// Error Evaluation
//			mean_square_residue = 0.0; // Initialize
//			for (k = 0, k_ext = 1; k < height; k++, k_ext++){
//				for (i = 0, i_ext = 1; i < length; i++, i_ext++){
//					for (j = 0, j_ext = 1; j < width; j++, j_ext++){
//
//						// 2D to 1D representation for i, j
//						x_ext = x_new(i_ext, j_ext, length_ext);
//						x = x_new(i, j, length);
//
//						mean_square_residue += (dataType)( pow( gauss_seidelPtr[k_ext][x_ext] * (1 + coef_tauh * (e_Ptr[k][x]
//							+ w_Ptr[k][x] + n_Ptr[k][x] + s_Ptr[k][x] + t_Ptr[k][x] + b_Ptr[k][x]))
//							+ coef_tauh * ( (prevSol_extPtr[k_ext][x_ext + 1] - prevSol_extPtr[k_ext][x_ext]) * VePtr[k][x]
//								+ (prevSol_extPtr[k_ext][x_ext - 1] - prevSol_extPtr[k_ext][x_ext]) * VwPtr[k][x]
//								+ (prevSol_extPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)] - prevSol_extPtr[k_ext][x_ext]) * VsPtr[k][x]
//								+ (prevSol_extPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)] - prevSol_extPtr[k_ext][x_ext]) * VnPtr[k][x]
//								+ (prevSol_extPtr[k_ext + 1][x_ext] - prevSol_extPtr[k_ext][x_ext]) * VbPtr[k][x]
//								+ (prevSol_extPtr[k_ext - 1][x_ext] - prevSol_extPtr[k_ext][x_ext]) * VtPtr[k][x] )
//							- coef_tauh * ( (e_Ptr[k][x] * gauss_seidelPtr[k_ext][x_ext + 1])
//								+ (w_Ptr[k][x] * gauss_seidelPtr[k_ext][x_ext - 1])
//								+ (s_Ptr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)])
//								+ (n_Ptr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)])
//								+ (b_Ptr[k][x] * gauss_seidelPtr[k_ext + 1][x_ext])
//								+ (t_Ptr[k][x] * gauss_seidelPtr[k_ext - 1][x_ext]) ) - prevSol_extPtr[k_ext][x_ext], 2) * hhh);
//					}
//				}
//			}
//		} while (mean_square_residue > segParameters.gauss_seidelTolerance && z < segParameters.maxNoGSIteration);
//		printf("The number of iterations is %zd\n", z);
//		printf("Residuum is %e\n", mean_square_residue);
//		printf("Step is %zd\n", segParameters.numberOfTimeStep);
//
//		//Copy the current time step to original data array after timeStepsNum
//		copyDataToReducedArea(segmFuntionPtr, gauss_seidelPtr, height, length, width);
//
//		//Rescale values of segmentation function and current time step to interval (0, 1)
//		rescaleToIntervalZeroOne(segmFuntionPtr, length, width, height);
//		rescaleToIntervalZeroOne(gauss_seidelPtr, length_ext, width_ext, height_ext);
//
//		secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
//
//		//Compute the L2 norm of the difference between the current and previous solutions
//		difference_btw_current_and_previous_sol = l2normD(prevSol_extPtr, gauss_seidelPtr, length_ext, width_ext, height_ext, segParameters.h);
//
//		printf("mass is %e\n", difference_btw_current_and_previous_sol);
//		printf("segTolerance is %e\n", segParameters.segTolerance);
//		printf("CPU time: %e secs\n", secondCpuTime - firstCpuTime);
//
//		//writing density.
//		if ((cpt % segParameters.mod) == 0)
//		{
//			strcpy_s(name, sizeof name, outputPathPtr);
//			sprintf_s(name_ending, sizeof(name_ending), "_seg_func_%03zd.vtk", cpt);
//			strcat_s(name, sizeof(name), name_ending);
//			pathsaveVTK = name;
//			vtkInfo->dataPointer = segmFuntionPtr;
//			storeVtkFile(pathsaveVTK, vtkInfo, dataForm);
//		}
//
//		cpt++;
//	} while ((cpt <= segParameters.maxNoOfTimeSteps) && (difference_btw_current_and_previous_sol > segParameters.segTolerance));
//
//
//	for (k = 0; k < height; k++)
//	{
//		free(segmFuntionPtr[k]);
//
//		free(GePtr[k]);
//		free(GwPtr[k]);
//		free(GnPtr[k]);
//		free(GsPtr[k]);
//		free(GtPtr[k]);
//		free(GbPtr[k]);
//
//		free(e_Ptr[k]);
//		free(w_Ptr[k]);
//		free(n_Ptr[k]);
//		free(s_Ptr[k]);
//		free(t_Ptr[k]);
//		free(b_Ptr[k]);
//
//		free(VePtr[k]);
//		free(VwPtr[k]);
//		free(VnPtr[k]);
//		free(VsPtr[k]);
//		free(VtPtr[k]);
//		free(VbPtr[k]);
//
//		free(edgeGradient[k]);
//	}
//	free(segmFuntionPtr);
//
//	free(GePtr);
//	free(GwPtr);
//	free(GnPtr);
//	free(GsPtr);
//	free(GtPtr);
//	free(GbPtr);
//
//	free(e_Ptr);
//	free(w_Ptr);
//	free(n_Ptr);
//	free(s_Ptr);
//	free(t_Ptr);
//	free(b_Ptr);
//
//	free(VePtr);
//	free(VwPtr);
//	free(VnPtr);
//	free(VsPtr);
//	free(VtPtr);
//	free(VbPtr);
//
//	free(edgeGradient);
//
//	for (k = 0; k < height_ext; k++)
//	{
//		free(prevSol_extPtr[k]);
//		free(gauss_seidelPtr[k]);
//		free(extendedCoefPtr[k]);
//	}
//	free(prevSol_extPtr);
//	free(gauss_seidelPtr);
//	free(extendedCoefPtr);
//
//	free(vtkInfo);
//}


bool generalizedSubsurfSegmentation(Image_Data inputImageData, dataType** initialSegment, Segmentation_Parameters segParameters, Filter_Parameters explicit_lhe_Parameters, unsigned char* outputPathPtr, dataType a, dataType b) {

	//test differentes branches
	size_t i, j, k;
	size_t dim2D = inputImageData.length * inputImageData.width;
	size_t height = inputImageData.height;
	size_t length = inputImageData.length;
	size_t width = inputImageData.width;
	size_t height_ext = inputImageData.height + 2;
	size_t length_ext = inputImageData.length + 2;
	size_t width_ext = inputImageData.width + 2;
	size_t dim2D_ext = length_ext * width_ext;
	dataType firstCpuTime, secondCpuTime, difference_btw_current_and_previous_sol;
	dataType** gauss_seidelPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	dataType** prevSol_extPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	dataType** imageToBeSegPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** segmFuntionPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** edgeGradientPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** GePtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** GwPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** GnPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** GsPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** GtPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** GbPtr = (dataType**)malloc(sizeof(dataType*) * height);
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
		GePtr == NULL || GwPtr == NULL || GnPtr == NULL || GsPtr == NULL || GtPtr == NULL || GbPtr == NULL || 
		e_Ptr == NULL || w_Ptr == NULL || n_Ptr == NULL || s_Ptr == NULL || t_Ptr == NULL || b_Ptr == NULL ||
		VePtr == NULL || VwPtr == NULL || VnPtr == NULL || VsPtr == NULL || VtPtr == NULL || VbPtr == NULL)
		return false;
	for (i = 0; i < height; i++)
	{
		imageToBeSegPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segmFuntionPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		edgeGradientPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		GePtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		GwPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		GnPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		GsPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		GtPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		GbPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
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
		if (imageToBeSegPtr[i] == NULL || segmFuntionPtr[i] == NULL || GePtr[i] == NULL || edgeGradientPtr[i] == NULL
			|| GwPtr[i] == NULL || GnPtr[i] == NULL || GsPtr[i] == NULL || GtPtr[i] == NULL || GbPtr[i] == NULL || e_Ptr[i] == NULL
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

	Gradient_Pointers VPtrs;
	VPtrs.GePtr = VePtr;
	VPtrs.GwPtr = VwPtr;
	VPtrs.GnPtr = VnPtr;
	VPtrs.GsPtr = VsPtr;
	VPtrs.GtPtr = VtPtr;
	VPtrs.GbPtr = VbPtr;

	//generate initial segmentation function
	//generateInitialSegmentationFunctionForMultipleCentres(segmFuntionPtr, length, width, height, centers, 0.5, 5, no_of_centers);
	copyDataToAnotherArray(initialSegment, segmFuntionPtr, height, length, width);

	//compute coefficients from presmoothed image
	generalizedGFunctionForImageToBeSegmented(inputImageData, prevSol_extPtr, edgeGradientPtr, GPtrs, VPtrs, segParameters, explicit_lhe_Parameters);

	Vtk_File_Info* vtkInfo = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	if (vtkInfo == NULL) return false;
	vtkInfo->spacing[0] = 1.0; vtkInfo->spacing[1] = 1.0; vtkInfo->spacing[2] = 1.0;
	vtkInfo->origin[0] = 0; vtkInfo->origin[1] = 0; vtkInfo->origin[2] = 0;
	vtkInfo->dimensions[0] = length; vtkInfo->dimensions[1] = width; vtkInfo->dimensions[2] = height;
	vtkInfo->vDataType = dta_Flt; vtkInfo->operation = copyTo; vtkDataForm dataForm = dta_binary;
	const char* pathsaveVTK;

	//Array for name construction
	unsigned char name[350];
	unsigned char name_ending[100];


	strcpy_s(name, sizeof name, outputPathPtr);
	sprintf_s(name_ending, sizeof(name_ending), "_g_est%03zd.vtk", i);
	strcat_s(name, sizeof(name), name_ending);
	pathsaveVTK = name;
	vtkInfo->dataPointer = GPtrs.GePtr;
	storeVtkFile(pathsaveVTK, vtkInfo, dataForm);

	//strcpy_s(name, sizeof name, outputPathPtr);
	//sprintf_s(name_ending, sizeof(name_ending), "_g_west.vtk", i);
	//strcat_s(name, sizeof(name), name_ending);
	//pathsaveVTK = name;
	//vtkInfo->dataPointer = GPtrs.GwPtr;
	//storeVtkFile(pathsaveVTK, vtkInfo, dataForm);

	//strcpy_s(name, sizeof name, outputPathPtr);
	//sprintf_s(name_ending, sizeof(name_ending), "_g_north.vtk", i);
	//strcat_s(name, sizeof(name), name_ending);
	//pathsaveVTK = name;
	//vtkInfo->dataPointer = GPtrs.GnPtr;
	//storeVtkFile(pathsaveVTK, vtkInfo, dataForm);

	//strcpy_s(name, sizeof name, outputPathPtr);
	//sprintf_s(name_ending, sizeof(name_ending), "_g_south.vtk", i);
	//strcat_s(name, sizeof(name), name_ending);
	//pathsaveVTK = name;
	//vtkInfo->dataPointer = GPtrs.GsPtr;
	//storeVtkFile(pathsaveVTK, vtkInfo, dataForm);

	//strcpy_s(name, sizeof name, outputPathPtr);
	//sprintf_s(name_ending, sizeof(name_ending), "_g_top.vtk", i);
	//strcat_s(name, sizeof(name), name_ending);
	//pathsaveVTK = name;
	//vtkInfo->dataPointer = GPtrs.GtPtr;
	//storeVtkFile(pathsaveVTK, vtkInfo, dataForm);

	//strcpy_s(name, sizeof name, outputPathPtr);
	//sprintf_s(name_ending, sizeof(name_ending), "_g_bottom.vtk", i);
	//strcat_s(name, sizeof(name), name_ending);
	//pathsaveVTK = name;
	//vtkInfo->dataPointer = GPtrs.GbPtr;
	//storeVtkFile(pathsaveVTK, vtkInfo, dataForm);

	//loop for segmentation time steps
	i = 1;
	do
	{
		segParameters.numberOfTimeStep = i;
		firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);

		//calcution of coefficients
		generalizedGaussSeidelCoefficients(prevSol_extPtr, imageData, edgeGradientPtr, CoefPtrs, segParameters);

		// Call to function that will evolve segmentation function in each discrete time step
		generalizedSubsurfSegmentationTimeStep(prevSol_extPtr, gauss_seidelPtr, imageData, GPtrs, VPtrs, segParameters, CoefPtrs, a, b);

		secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);

		//Compute the L2 norm of the difference between the current and previous solutions
		difference_btw_current_and_previous_sol = l2normD(prevSol_extPtr, gauss_seidelPtr, length, width, height, segParameters.h);

		printf("mass is %e\n", difference_btw_current_and_previous_sol);
		printf("segTolerance is %e\n", segParameters.segTolerance);
		printf("CPU time: %e secs\n", secondCpuTime - firstCpuTime);

		//writing density.
		if ((i % segParameters.mod) == 0)
		{
			strcpy_s(name, sizeof name, outputPathPtr);
			sprintf_s(name_ending, sizeof(name_ending), "_seg_func_%03zd.vtk", i);
			strcat_s(name, sizeof(name), name_ending);
			pathsaveVTK = name;
			vtkInfo->dataPointer = imageData.segmentationFuntionPtr;
			storeVtkFile(pathsaveVTK, vtkInfo, dataForm);
		}
		i++;
	} while ((i <= segParameters.maxNoOfTimeSteps) && (difference_btw_current_and_previous_sol > segParameters.segTolerance));

	printf("finish: Segmentation tolerance is %lf\n", segParameters.segTolerance);

	for (i = 0; i < height; i++)
	{
		free(imageToBeSegPtr[i]);
		free(segmFuntionPtr[i]);
		free(edgeGradientPtr[i]);
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

	free(vtkInfo);

	return true;
}

bool generalizedGFunctionForImageToBeSegmented(Image_Data inputImageData, dataType** extendedCoefPtr, dataType** edgeGradientPtr, Gradient_Pointers GPtrs, Gradient_Pointers VPtrs,
	Segmentation_Parameters segParameters, Filter_Parameters explicit_lhe_Parameters)
{
	//checks if the memory was allocated
	if (inputImageData.imageDataPtr == NULL || extendedCoefPtr == NULL || edgeGradientPtr == NULL || GPtrs.GePtr == NULL || GPtrs.GwPtr == NULL
		|| GPtrs.GnPtr == NULL || GPtrs.GsPtr == NULL || GPtrs.GtPtr == NULL || GPtrs.GbPtr == NULL || VPtrs.GePtr == NULL || VPtrs.GwPtr == NULL
		|| VPtrs.GnPtr == NULL || VPtrs.GsPtr == NULL || VPtrs.GtPtr == NULL || VPtrs.GbPtr == NULL)
		return false;

	dataType** imageToBeSegPtr = inputImageData.imageDataPtr;
	size_t i, j, k, x, x_ext; // length == xDim, width == yDim, height == zDim
	size_t kplus1, kminus1, iminus1, iplus1, jminus1, jplus1;
	size_t dim2D = inputImageData.length * inputImageData.width;
	size_t k_ext, j_ext, i_ext;
	size_t height_ext = inputImageData.height + 2;
	size_t length_ext = inputImageData.length + 2;
	size_t width_ext = inputImageData.width + 2;
	dataType quotient = (dataType)(4.0 * segParameters.h);
	dataType ux, uy, uz; //change in x, y and z respectively
	dataType u, uN, uS, uE, uW, uNW, uNE, uSE, uSW, Tu, TuN, TuS, TuE, TuW, TuNW, TuNE, TuSE, TuSW, //current and surrounding voxel values
		Bu, BuN, BuS, BuE, BuW, BuNW, BuNE, BuSE, BuSW;
	dataType norm_image_smoothed_e, norm_image_smoothed_w, norm_image_smoothed_n, norm_image_smoothed_s, norm_image_smoothed_t, norm_image_smoothed_b;
	dataType norm_image_smoothed_average;

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
	//heatExplicitScheme(presmoothingData, explicit_lhe_Parameters);
	heatImplicitScheme(presmoothingData, explicit_lhe_Parameters);
	//geodesicMeanCurvatureTimeStep(presmoothingData, explicit_lhe_Parameters);
	//meanCurvatureTimeStep(presmoothingData, explicit_lhe_Parameters);

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

				// Calculation of coefficients in east direction
				ux = (uE - u) / segParameters.h;
				uy = ((uN + uNE) - (uS + uSE))
					/ quotient;
				uz = ((Tu + TuE) - (Bu + BuE))
					/ quotient;
				GPtrs.GePtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);

				norm_image_smoothed_e = sqrt((ux * ux) + (uy * uy) + (uz * uz));

				// Calculation of coefficients in west direction
				ux = (uW - u) / segParameters.h;
				uy = ((uNW + uN) - (uSW + uS))
					/ quotient;
				uz = ((TuW + Tu) - (BuW + Bu))
					/ quotient;
				GPtrs.GwPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);

				norm_image_smoothed_w = sqrt((ux * ux) + (uy * uy) + (uz * uz));

				// Calculation of coefficients in north direction
				ux = ((uNE + uE) - (uNW + uW))
					/ quotient;
				uy = (uN - u) / segParameters.h;
				uz = ((TuN + Tu) - (BuN + Bu))
					/ quotient;
				GPtrs.GnPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);

				norm_image_smoothed_n = sqrt((ux * ux) + (uy * uy) + (uz * uz));

				// Calculation of coefficients in south direction
				ux = ((uE + uSE) - (uW + uSW))
					/ quotient;
				uy = (uS - u) / segParameters.h;
				uz = ((TuS + Tu) - (BuS + Bu))
					/ quotient;
				GPtrs.GsPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);

				norm_image_smoothed_s = sqrt((ux * ux) + (uy * uy) + (uz * uz));

				// Calculation of coefficients in top direction
				ux = ((TuE + uE) - (TuW + uW))
					/ quotient;
				uy = ((TuN + uN) - (TuS + uS))
					/ quotient;
				uz = (Tu - u) / segParameters.h;
				GPtrs.GtPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);

				norm_image_smoothed_t = sqrt((ux * ux) + (uy * uy) + (uz * uz));

				// Calculation of coefficients in bottom direction
				ux = ((BuW + uW) - (BuE + uE))
					/ quotient;
				uy = ((BuN + uN) - (BuS + uS))
					/ quotient;
				uz = (Bu - u) / segParameters.h;
				GPtrs.GbPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);

				norm_image_smoothed_b = sqrt((ux * ux) + (uy * uy) + (uz * uz));

				norm_image_smoothed_average = (norm_image_smoothed_e + norm_image_smoothed_w + norm_image_smoothed_n + norm_image_smoothed_s +
					norm_image_smoothed_t + norm_image_smoothed_b) / 6.0;

				edgeGradientPtr[k][x] = gradientFunction(norm_image_smoothed_average, segParameters.coef);

				VPtrs.GePtr[k][x] = (GPtrs.GePtr[k][x] - edgeGradientPtr[k][x]) / (2 * segParameters.h);
				VPtrs.GwPtr[k][x] = (GPtrs.GwPtr[k][x] - edgeGradientPtr[k][x]) / (2 * segParameters.h);
				VPtrs.GnPtr[k][x] = (GPtrs.GnPtr[k][x] - edgeGradientPtr[k][x]) / (2 * segParameters.h);
				VPtrs.GsPtr[k][x] = (GPtrs.GsPtr[k][x] - edgeGradientPtr[k][x]) / (2 * segParameters.h);
				VPtrs.GtPtr[k][x] = (GPtrs.GtPtr[k][x] - edgeGradientPtr[k][x]) / (2 * segParameters.h);
				VPtrs.GbPtr[k][x] = (GPtrs.GbPtr[k][x] - edgeGradientPtr[k][x]) / (2 * segParameters.h);
			
			}
		}
	}
	return true;
}

bool generalizedGaussSeidelCoefficients(dataType** extendedCoefPtr, Segment_Image_Data inputImageData, dataType** edgeGradientPtr, Coefficient_Pointers CoefPtrs, Segmentation_Parameters segParameters)
{
	//checks if the memory was allocated
	if (extendedCoefPtr == NULL || inputImageData.segmentationFuntionPtr == NULL || edgeGradientPtr == NULL 
		|| CoefPtrs.w_Ptr == NULL || CoefPtrs.n_Ptr == NULL || CoefPtrs.s_Ptr == NULL || CoefPtrs.t_Ptr == NULL || CoefPtrs.b_Ptr == NULL)
		return false;

	size_t i, j, k, x, x_ext; // length == xDim, width == yDim, height == zDim
	size_t kplus1, kminus1, iminus1, iplus1, jminus1, jplus1;
	size_t dim2D = inputImageData.length * inputImageData.width;
	size_t k_ext, j_ext, i_ext;
	size_t height_ext = inputImageData.height + 2;
	size_t length_ext = inputImageData.length + 2;
	size_t width_ext = inputImageData.width + 2;
	dataType quotient = (dataType)(4.0 * segParameters.h);
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
				orig_e = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in west direction
				orig_ux = (orig_uW - orig_u) / segParameters.h;
				orig_uy = ((orig_uNW + orig_uN) - (orig_uSW + orig_uS))
					/ quotient;
				orig_uz = ((orig_TuW + orig_Tu) - (orig_BuW + orig_Bu))
					/ quotient;
				orig_w = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in north direction
				orig_ux = ((orig_uNE + orig_uE) - (orig_uNW + orig_uW))
					/ quotient;
				orig_uy = (orig_uN - orig_u) / segParameters.h;
				orig_uz = ((orig_TuN + orig_Tu) - (orig_BuN + orig_Bu))
					/ quotient;
				orig_n = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in south direction
				orig_ux = ((orig_uE + orig_uSE) - (orig_uW + orig_uSW))
					/ quotient;
				orig_uy = (orig_uS - orig_u) / segParameters.h;
				orig_uz = ((orig_TuS + orig_Tu) - (orig_BuS + orig_Bu))
					/ quotient;
				orig_s = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in top direction
				orig_ux = ((orig_TuE + orig_uE) - (orig_TuW + orig_uW))
					/ quotient;
				orig_uy = ((orig_TuN + orig_uN) - (orig_TuS + orig_uS))
					/ quotient;
				orig_uz = (orig_Tu - orig_u) / segParameters.h;
				orig_t = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in bottom direction
				orig_ux = ((orig_BuW + orig_uW) - (orig_BuE + orig_uE))
					/ quotient;
				orig_uy = ((orig_BuN + orig_uN) - (orig_BuS + orig_uS))
					/ quotient;
				orig_uz = (orig_Bu - orig_u) / segParameters.h;
				orig_b = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// evaluation of norm of gradient of image at each voxel
				average_face_coef = (dataType)(((orig_e + orig_w + orig_n + orig_s + orig_t + orig_b) / 6.0));

				voxel_coef = (dataType)sqrt(pow(average_face_coef, 2) + segParameters.eps2);

				/* evaluation of norm of gradient of image at each voxel, norm of gradient of presmoothed
				image at each voxel face and reciprocal of norm of gradient of image at each voxel face*/
				CoefPtrs.e_Ptr[k][x] = (dataType)(voxel_coef * edgeGradientPtr[k][x] * (1.0 / orig_e));//east coefficient
				CoefPtrs.w_Ptr[k][x] = (dataType)(voxel_coef * edgeGradientPtr[k][x] * (1.0 / orig_w));//west coefficient
				CoefPtrs.n_Ptr[k][x] = (dataType)(voxel_coef * edgeGradientPtr[k][x] * (1.0 / orig_n));//north coefficient
				CoefPtrs.s_Ptr[k][x] = (dataType)(voxel_coef * edgeGradientPtr[k][x] * (1.0 / orig_s));//south coefficient
				CoefPtrs.t_Ptr[k][x] = (dataType)(voxel_coef * edgeGradientPtr[k][x] * (1.0 / orig_t));//top coefficient
				CoefPtrs.b_Ptr[k][x] = (dataType)(voxel_coef * edgeGradientPtr[k][x] * (1.0 / orig_b));//bottom coefficient
			}
		}
	}

	return true;
}

bool generalizedSubsurfSegmentationTimeStep(dataType** prevSol_extPtr, dataType** gauss_seidelPtr, Segment_Image_Data inputImageData, Gradient_Pointers GPtrs, Gradient_Pointers VPtrs,
	Segmentation_Parameters segParameters, Coefficient_Pointers CoefPtrs, dataType a, dataType b)
{
	//check if the memory was allocated successfully
	if (inputImageData.segmentationFuntionPtr == NULL || prevSol_extPtr == NULL || gauss_seidelPtr == NULL || GPtrs.GePtr == NULL || GPtrs.GwPtr == NULL
		|| GPtrs.GnPtr == NULL || GPtrs.GsPtr == NULL || GPtrs.GtPtr == NULL || GPtrs.GbPtr == NULL || CoefPtrs.e_Ptr == NULL || CoefPtrs.w_Ptr == NULL
		|| CoefPtrs.n_Ptr == NULL || CoefPtrs.s_Ptr == NULL || CoefPtrs.t_Ptr == NULL || CoefPtrs.b_Ptr == NULL
		|| VPtrs.GePtr == NULL || VPtrs.GwPtr == NULL || VPtrs.GnPtr == NULL || VPtrs.GsPtr == NULL || VPtrs.GtPtr == NULL || VPtrs.GbPtr == NULL)
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

	const dataType coef_tauh = tau / hhh;

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

					// Gauss-Seidel Formula Evaluation
					gauss_seidel = (prevSol_extPtr[k_ext][x_ext] -
						coef_tauh * a * ( (prevSol_extPtr[k_ext][x_ext + 1] - prevSol_extPtr[k_ext][x_ext]) * VPtrs.GePtr[k][x] +
						(prevSol_extPtr[k_ext][x_ext - 1] - prevSol_extPtr[k_ext][x_ext]) * VPtrs.GwPtr[k][x] +
						(prevSol_extPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)] - prevSol_extPtr[k_ext][x_ext]) * VPtrs.GsPtr[k][x] +
						(prevSol_extPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)] - prevSol_extPtr[k_ext][x_ext]) * VPtrs.GnPtr[k][x] +
						(prevSol_extPtr[k_ext + 1][x_ext] - prevSol_extPtr[k_ext][x_ext]) * VPtrs.GbPtr[k][x] +
						(prevSol_extPtr[k_ext - 1][x_ext] - prevSol_extPtr[k_ext][x_ext]) * VPtrs.GtPtr[k][x] ) +
						coef_tauh * b * ( (CoefPtrs.e_Ptr[k][x] * gauss_seidelPtr[k_ext][x_ext + 1])
						+ (CoefPtrs.w_Ptr[k][x] * gauss_seidelPtr[k_ext][x_ext - 1])
						+ (CoefPtrs.s_Ptr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)])
						+ (CoefPtrs.n_Ptr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)])
						+ (CoefPtrs.b_Ptr[k][x] * gauss_seidelPtr[k_ext + 1][x_ext])
						+ (CoefPtrs.t_Ptr[k][x] * gauss_seidelPtr[k_ext - 1][x_ext])) ) /
						(1 + coef_tauh * b * (CoefPtrs.e_Ptr[k][x] + CoefPtrs.w_Ptr[k][x] + CoefPtrs.n_Ptr[k][x]
							+ CoefPtrs.s_Ptr[k][x] + CoefPtrs.t_Ptr[k][x] + CoefPtrs.b_Ptr[k][x]));

					// SOR implementation using Gauss-Seidel
					gauss_seidelPtr[k_ext][x_ext] = gauss_seidelPtr[k_ext][x_ext] +
						segParameters.omega_c * (gauss_seidel - gauss_seidelPtr[k_ext][x_ext]);
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

					mean_square_residue += (dataType)(pow(gauss_seidelPtr[k_ext][x_ext] * (1 + coef_tauh * b * (CoefPtrs.e_Ptr[k][x]
						+ CoefPtrs.w_Ptr[k][x] + CoefPtrs.n_Ptr[k][x] + CoefPtrs.s_Ptr[k][x]
						+ CoefPtrs.t_Ptr[k][x] + CoefPtrs.b_Ptr[k][x]))
						+ coef_tauh * a * ( (prevSol_extPtr[k_ext][x_ext + 1] - prevSol_extPtr[k_ext][x_ext]) * VPtrs.GePtr[k][x] 
							+ (prevSol_extPtr[k_ext][x_ext - 1] - prevSol_extPtr[k_ext][x_ext]) * VPtrs.GwPtr[k][x]
							+ (prevSol_extPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)] - prevSol_extPtr[k_ext][x_ext]) * VPtrs.GsPtr[k][x]
							+ (prevSol_extPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)] - prevSol_extPtr[k_ext][x_ext]) * VPtrs.GnPtr[k][x]
							+ (prevSol_extPtr[k_ext + 1][x_ext] - prevSol_extPtr[k_ext][x_ext]) * VPtrs.GbPtr[k][x]
							+ (prevSol_extPtr[k_ext - 1][x_ext] - prevSol_extPtr[k_ext][x_ext]) * VPtrs.GtPtr[k][x] )
						- coef_tauh * b * ((CoefPtrs.e_Ptr[k][x] * gauss_seidelPtr[k_ext][x_ext + 1])
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

	////Rescale values of segmentation function and current time step to interval (0, 1)
	rescaleToIntervalZeroOne(inputImageData.segmentationFuntionPtr, length, width, height);
	rescaleToIntervalZeroOne(gauss_seidelPtr, length_ext, width_ext, height_ext);
	
	//if (no_of_centers == 1)
	//{
	//	rescaleToIntervalZeroOne(inputImageData.segmentationFuntionPtr, length, width, height);
	//	rescaleToIntervalZeroOne(gauss_seidelPtr, length_ext, width_ext, height_ext);
	//}
	//else
	//{
	//	for (i = 0; i < no_of_centers; i++)
	//	{
	//		rescaleLocallyToIntervalZeroOne(inputImageData.segmentationFuntionPtr, length, width, height,
	//			centers[i].x, centers[i].y, centers[i].z, 6., i);
	//		rescaleLocallyToIntervalZeroOne(gauss_seidelPtr, length_ext, width_ext, height_ext,
	//			centers[i].x, centers[i].y, centers[i].z, 6., i);
	//	}
	//}

	return true;
}
