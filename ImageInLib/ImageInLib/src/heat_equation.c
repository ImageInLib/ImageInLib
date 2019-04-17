#include <stdlib.h>
#include<stdio.h>
#include <stdbool.h> // Boolean function bool
#include <math.h> // Maths functions i.e. pow, sin, cos
#include "heat_equation.h"

// Local Function Prototype

void heatExplicitScheme(Image_Data toExplicitImage, const Filter_Parameters explicitParameters)
{
	size_t k, i, j;
	dataType hh = explicitParameters.h*explicitParameters.h;
	dataType tau = explicitParameters.timeStepSize;

	// Perform Reflection of the tempPtr
	// Prepare variables toExplicitImage.height, toExplicitImage.length, toExplicitImage.width
	// Less the borders because in the loops we add back the border p
	size_t height = toExplicitImage.height, length = toExplicitImage.length, width = toExplicitImage.width;
	size_t height_ext = height + 2;
	size_t length_ext = length + 2;
	size_t width_ext = width + 2;

	// Create temp Image Data holder for Previous time step data - with extended boundary because of boundary condition

	dataType ** tempPtr = (dataType **)malloc(sizeof(dataType *) * (height_ext));
	for (k = 0; k < ((height + 2)); k++)
	{
		tempPtr[k] = malloc(sizeof(dataType)*(length_ext)*(width_ext));
	}

	size_t k_ext, j_ext, i_ext;
	size_t sliceBound = (length_ext - 1)* width_ext;
	size_t i_d = 0;

	copyDataToExtendedArea(toExplicitImage.imageDataPtr, tempPtr, height, length, width);
	reflection3D(tempPtr, height_ext, length_ext, width_ext);

	size_t x;
	size_t x_ext;

	const dataType coeff = tau / hh;

	// The Explicit Scheme Evaluation
	for (size_t t = 0; t < explicitParameters.timeStepsNum; t++)
	{
		for (k = 0, k_ext = 1; k < height; k++, k_ext++)
		{
			for (i = 0, i_ext = 1; i < length; i++, i_ext++)
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++)
				{
					// 2D to 1D representation for i, j
					x_ext = x_new(i_ext, j_ext, length_ext);
					x = x_new(i, j, length);

					// Explicit formula
					toExplicitImage.imageDataPtr[k][x] = (1.0 - 6.0 * coeff)*tempPtr[k_ext][x_ext]
						+ coeff * (tempPtr[k_ext][x_ext + 1]
							+ tempPtr[k_ext][x_ext - 1]
							+ tempPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)]
							+ tempPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)]
							+ tempPtr[k_ext + 1][x_ext]
							+ tempPtr[k_ext - 1][x_ext]);
				}
			}
		}

		copyDataToExtendedArea(toExplicitImage.imageDataPtr, tempPtr, height, length, width);
		reflection3D(tempPtr, height_ext, length_ext, width_ext);
	}
	// Freeing Memory after use
	for (i = 0; i < height_ext; i++)
	{
		free(tempPtr[i]);
	}
	free(tempPtr);
}

void heatImplicitScheme(Image_Data toImplicitImage, const Filter_Parameters implicitParameters)
{
	size_t k, i, j, z, steps = implicitParameters.timeStepsNum, p = implicitParameters.p;
	dataType hhh = implicitParameters.h*implicitParameters.h*implicitParameters.h;
	dataType coeff = implicitParameters.timeStepSize / hhh;
	// Error value used to check iteration
	// sor - successive over relation value, used in Gauss-Seidel formula
	dataType error, sor;
	// Prepare variables toExplicitImage.height, toExplicitImage.length, toExplicitImage.width
	// Less the borders because in the loops we add back the border p
	size_t height = toImplicitImage.height, length = toImplicitImage.length, width = toImplicitImage.width;
	size_t height_ext = height + 2;
	size_t length_ext = length + 2;
	size_t width_ext = width + 2;// Create temp Image Data holder for Previous time step data
	dataType ** tempPtr = (dataType **)malloc(sizeof(dataType *) * (height_ext));
	dataType ** currentPtr = (dataType **)malloc(sizeof(dataType *) * (height_ext)); // holds current
	for (i = 0; i < ((height_ext)); i++)
	{
		tempPtr[i] = malloc(sizeof(dataType)*(length_ext)*(width_ext));
		currentPtr[i] = malloc(sizeof(dataType)*(length_ext)*(width_ext));
	}
	copyDataToExtendedArea(toImplicitImage.imageDataPtr, tempPtr, height, length, width);
	copyDataToExtendedArea(toImplicitImage.imageDataPtr, currentPtr, height, length, width);
	// Perform Reflection of the tempPtr
	reflection3D(tempPtr, height_ext, length_ext, width_ext);
	reflection3D(currentPtr, height_ext, length_ext, width_ext);
	//reflection3DB(tempPtr, height, length, width,p);
	// The Gauss-Seidel Implicit Scheme
	size_t k_ext, j_ext, i_ext, x_ext, x;
	z = 0; // Steps counter
	do
	{
		z = z + 1;
		// Gauss-Seidel Method
		for (k = 0, k_ext = 1; k < height; k++, k_ext++)
		{
			for (i = 0, i_ext = 1; i < length; i++, i_ext++)
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++)
				{
					// 2D to 1D representation for i, j
					x_ext = x_new(i_ext, j_ext, length_ext);
					// Begin Gauss-Seidel Formula Evaluation
					sor = (tempPtr[k_ext][x_ext] + coeff * (currentPtr[k_ext][x_new(i_ext + 1, j_ext, length_ext)]
						+ currentPtr[k_ext][x_new(i_ext - 1, j_ext, length_ext)]
						+ currentPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)]
						+ currentPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)]
						+ currentPtr[k_ext + 1][x_ext] + currentPtr[k_ext - 1][x_ext])) / (1 + 6.0 * coeff);
					// Gauss-Seidel
					currentPtr[k_ext][x_ext] = currentPtr[k_ext][x_ext] + implicitParameters.omega_c*(sor - currentPtr[k_ext][x_ext]);
				}
			}
		}
		// Error Evaluation
		error = 0.0; // Initialize
		//reflection3DB(tempPtr, height, length, width, p);
		for (k = 0, k_ext = 1; k < height; k++, k_ext++)
		{
			for (i = 0, i_ext = 1; i < length; i++, i_ext++)
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++)
				{
					// 2D to 1D representation for i, j
					x_ext = x_new(i_ext, j_ext, length_ext);
					// Begin Error Calculation
					error += pow(currentPtr[k_ext][x_ext] * (1 + 6.0*coeff)
						- coeff * (currentPtr[k_ext][x_new(i_ext + 1, j_ext, length_ext)]
							+ currentPtr[k_ext][x_new(i_ext - 1, j_ext, length_ext)] + currentPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)]
							+ currentPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)] + currentPtr[k_ext + 1][x_ext]
							+ currentPtr[k_ext - 1][x_ext]) - tempPtr[k_ext][x_ext], 2);
				}
			}
		}
	} while (error > implicitParameters.tolerance && z < steps);

	printf("The number of iterations is %zd\n", z);
	printf("Error is %e\n", error);

	// Copy back to original after filter
	for (k = 0, k_ext = 1; k < height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				// 2D to 1D representation for i, j
				x_ext = x_new(i_ext, j_ext, length_ext);
				x = x_new(i, j, length);
				toImplicitImage.imageDataPtr[k][x] = currentPtr[k_ext][x_ext];
			}
		}
	}
	// Freeing Memory after use
	for (i = 0; i < ((height + 2)); i++)
	{
		free(tempPtr[i]);
		free(currentPtr[i]);
	}
	free(tempPtr);
	free(currentPtr);
}