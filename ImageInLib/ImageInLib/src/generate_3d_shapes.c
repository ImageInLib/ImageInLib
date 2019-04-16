#include "generate_3d_shapes.h"

#include <stdio.h> // Standard lib for input and output functions
#include <stdbool.h> // Boolean Function
#include <math.h>
#include "common_math.h"
#include "common_functions.h"
#include "data_storage.h"
#include "data_initialization.h"

/*
* Function to fill the points of Cuboid (block)
*/
void fillBlock3D(double **inputDataArrayPtr, size_t inputHeight, size_t inputLength, size_t inputWidth, Point3D blockCorner, double *fillBlockDimension, double fillValue)
{
	if (inputHeight < 1 || inputWidth < 1 || inputLength < 1)
		return;

	// Variables to be used to loop
	size_t k, i, j;
	// Begin points
	size_t k_0 = (size_t)blockCorner.z, i_0 = (size_t)blockCorner.x, j_0 = (size_t)blockCorner.y;
	// Block Distance
	double gK = fillBlockDimension[2], gI = fillBlockDimension[0], gJ = fillBlockDimension[1];
	// Begin Fill

	const size_t k_max = (size_t)fmin((double)k_0 + gK - 1, (double)inputHeight - 1);
	const size_t j_max = (size_t)fmin((double)j_0 + gJ - 1, (double)inputWidth - 1);
	const size_t i_max = (size_t)fmin((double)i_0 + gI - 1, (double)inputLength - 1);

	for (k = k_0; k <= k_max; k++)
	{
		for (i = i_0; i <= i_max; i++)
		{
			for (j = j_0; j <= j_max; j++)
			{
				// 1D representation
				size_t x = x_new(i, j, inputLength);
				// Fill Value
				inputDataArrayPtr[k][x] = fillValue;
			}
		}
	}
}
/*
* Function to fill 3D Ball Shape Points
*/
void fillBall3D(double **inputDataArrayPtr, size_t inputHeight, size_t inputLength, size_t inputWidth, double sphereRadius, Point3D sphereCenter, double fillValue)
{
	// Variables to be used to loop
	size_t k, i, j;
	// Center points
	size_t k_0 = (size_t)sphereCenter.z, i_0 = (size_t)sphereCenter.x, j_0 = (size_t)sphereCenter.y;
	// Difference
	size_t k_d, i_d, j_d;
	// Begin Fill
	for (k = 0; k < inputHeight; k++)
	{
		k_d = k - k_0;
		for (i = 0; i < inputLength; i++)
		{
			i_d = i - i_0;
			for (j = 0; j < inputWidth; j++)
			{
				j_d = j - j_0;
				// Calculates distance
				double distance = sqrt((double)(k_d*k_d + i_d * i_d + j_d * j_d));
				// Check if within Given radius
				if (distance <= sphereRadius)
				{
					// 1D representation
					size_t x = x_new(i, j, inputLength);
					// Fill Value
					inputDataArrayPtr[k][x] = fillValue;
				}
			}
		}
	}
}

bool generateSphereWithSixHoles(dataType ** dataArray3D, Point3D center, size_t length, size_t width, size_t height,
	dataType sphereRadius, dataType smallRadius, dataType fillValue, unsigned char * outputPathPtr)
{
	//checks if the memory was allocated
	if (dataArray3D == NULL)
		return false;

	size_t i, j, k, s;
	dataType phi, theta, r = sphereRadius;
	dataType point_in_circleX, point_in_circleY, point_in_circleZ;

	center.x = (dataType)(length / 2);
	center.y = (dataType)(width / 2);
	center.z = (dataType)(height / 2);

	//Initialization of arrays
	if (fillValue == 255)
	{
		initialize3dArrayD(dataArray3D, length, width, height, 0.);
	}
	else
	{
		initialize3dArrayD(dataArray3D, length, width, height, 255.);
	}

	double step = M_PI / 1000.0;
	for (phi = 0.0; phi < 2 * M_PI; phi += step) {
		for (theta = 0.0; theta < M_PI; theta += step) {
			i = (size_t)(center.x + r * cos(phi)*sin(theta));// + 0.5
			j = (size_t)(center.y + r * sin(phi)*sin(theta));
			k = (size_t)(center.z + r * cos(theta));

			// Calculates distance
			point_in_circleX = sqrt(((i - center.x)*(i - center.x) + (j - center.y)*(j - center.y)));
			point_in_circleY = sqrt(((i - center.x)*(i - center.x) + (k - center.z)*(k - center.z)));
			point_in_circleZ = sqrt(((k - center.z)*(k - center.z) + (j - center.y)*(j - center.y)));
			// 2D to 1D representation
			s = x_new(i, j, length);
			// Insertion of tube to through the sphere
			if (fillValue == 255)
			{
				if (((point_in_circleX <= smallRadius) || (point_in_circleY <= smallRadius) || (point_in_circleZ <= smallRadius)))
					dataArray3D[k][s] = 0.;
				else
					dataArray3D[k][s] = fillValue;
			}
			else
			{
				if (((point_in_circleX <= smallRadius) || (point_in_circleY <= smallRadius) || (point_in_circleZ <= smallRadius)))
					dataArray3D[k][s] = 255.;
				else
					dataArray3D[k][s] = fillValue;
			}
		}
	}

	// Store generated object to file
	Storage_Flags flags = { true, true };
	store3dDataVtkD(dataArray3D, length, width, height, outputPathPtr, (2.5 / (length)), flags);

	return true;
}

bool generateSphere(dataType ** dataArray3D, Point3D center, size_t length, size_t width, size_t height,
	dataType sphereRadius, dataType fillValue, unsigned char * outputPathPtr)
{
	//checks if the memory was allocated
	if (dataArray3D == NULL)
		return false;

	size_t i, j, k, s;
	dataType phi, theta, r = sphereRadius;

	center.x = (dataType)(length / 2);
	center.y = (dataType)(width / 2);
	center.z = (dataType)(height / 2);

	//Initialization of arrays
	if (fillValue == 255)
	{
		initialize3dArrayD(dataArray3D, length, width, height, 0.);
	}
	else
	{
		initialize3dArrayD(dataArray3D, length, width, height, 255.);
	}

	double step = M_PI / 1000.0;
	for (phi = 0.0; phi < 2 * M_PI; phi += step) {
		for (theta = 0.0; theta < M_PI; theta += step) {
			i = (size_t)(center.x + r * cos(phi)*sin(theta));// + 0.5
			j = (size_t)(center.y + r * sin(phi)*sin(theta));
			k = (size_t)(center.z + r * cos(theta));

			// 2D to 1D representation
			s = x_new(i, j, length);

			dataArray3D[k][s] = fillValue;
		}
	}

	// Store generated object to file
	Storage_Flags flags = { true, true };
	store3dDataVtkD(dataArray3D, length, width, height, outputPathPtr, (2.5 / (length)), flags);

	return true;
}