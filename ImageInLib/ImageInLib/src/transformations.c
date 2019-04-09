#include <stdlib.h> // Malloc and other functions
#include <math.h> // Mathematical functions
#include <stdbool.h>
#include "transformations.h"
//==============================================================================
// Local Functions Prototype
/*
* Rotated indices
*/
dataType x_rotate(dataType z, dataType x, dataType y, dataType theta, dataType psi);
dataType y_rotate(dataType z, dataType x, dataType y, dataType theta, dataType psi, dataType phi);
dataType z_rotate(dataType z, dataType x, dataType y, dataType theta, dataType psi, dataType phi);
//==============================================================================
// Inverse
dataType x_rotateInv(dataType z, dataType x, dataType y, dataType theta, dataType psi, dataType phi);
dataType y_rotateInv(dataType z, dataType x, dataType y, dataType theta, dataType psi, dataType phi);
dataType z_rotateInv(dataType z, dataType x, dataType y, dataType theta, dataType psi, dataType phi);
//==============================================================================
/*
* Function to calculate and return interpolated values
*/
dataType interpolated(dataType k_t, dataType i_t, dataType j_t, int top, int bottom, int left, int right, int begin, int end, dataType ** imageDataPtr, size_t imageWidth);
/*
* Transform Function for imageDataPtr
*/
void transform3DImage(dataType ** sourceDataPtr, dataType ** imageDataPtr, Point3D translation, Point3D scaling, Point3D rotation, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType bgValue, dataType centroid[3])
{
	size_t k, i, j;
	// Creates a new Pointer to fill the transformed values
	dataType ** transformPointsPtr = (dataType **)malloc(imageHeight * sizeof(dataType *));

	dataType k_a, i_a, j_a; // Affine indices
	// Rotation Angles -
	//dataType theta = (rotation.y*M_PI) / 180, psi = (rotation.z*M_PI) / 180, phi = (rotation.x*M_PI) / 180;
	dataType theta = (rotation.y), psi = (rotation.z), phi = (rotation.x);
	// Center Points
	//int hcenter = floor(imageHeight / 2), lcenter = floor(imageLength / 2), wcenter = floor(imageWidth / 2);
	dataType cz = centroid[2], cx = centroid[0], cy = centroid[1];

	// Dimension variable
	// dimWidth = X*Y
	size_t dimWidth = imageLength * imageWidth;

	// Temporary parameters
	dataType tmpX, tmpY, tmpZ, tmp;
	// Transformed
	dataType k_t, i_t, j_t; // Transformed indices
	for (k = 0; k < imageHeight; k++)
	{
		transformPointsPtr[k] = (dataType *)malloc(dimWidth * sizeof(dataType));
		k_a = k - cz; // Move to origin Z
		for (i = 0; i < imageLength; i++)
		{
			i_a = i - cx; // Move to origin x
			for (j = 0; j < imageWidth; j++)
			{
				// 2D to 1D representation for i, j
				size_t x = x_new(i, j, imageLength);

				j_a = j - cy; // Move to origin Y

							  // Apply scaling
				tmpZ = k_a / scaling.z;
				tmpX = i_a / scaling.x;
				tmpY = j_a / scaling.y;

				// Apply Rotation

				// 3. Rotation - No rotation
				i_t = x_rotate(tmpZ, tmpX, tmpY, theta, psi);
				j_t = y_rotate(tmpZ, tmpX, tmpY, theta, psi, phi);
				k_t = z_rotate(tmpZ, tmpX, tmpY, theta, psi, phi);

				// Move back to centroid
				tmpX = i_t + cx;
				tmpY = j_t + cy;
				tmpZ = k_t + cz;

				// Set the values
				i_t = tmpX;
				j_t = tmpY;
				k_t = tmpZ;

				// Add translation
				i_t = i_t - translation.x;
				j_t = j_t - translation.y;
				k_t = k_t - translation.z;

				// Use Interpolation to get the values
				// Locations for Tri-linear Interpolation
				// Z
				int bottom = (int)floor(k_t);
				int top = bottom + 1;
				// X
				int left = (int)floor(i_t);
				int right = left + 1;
				// Y
				int begin = (int)floor(j_t);
				int end = begin + 1;
				// Check if within limits
				if (bottom >= 0 && top < imageHeight && left >= 0 && right < imageLength && begin >= 0 && end < imageWidth)
				{
					tmp = interpolated(k_t, i_t, j_t, top, bottom, left, right, begin, end, sourceDataPtr, imageLength);
					transformPointsPtr[k][x] = tmp;
				}
				else
				{
					transformPointsPtr[k][x] = bgValue; // Background value
				}
			}
		}
	}
	// Copy Transformed back to Image Data Pointer
	for (k = 0; k < imageHeight; k++)
	{
		for (i = 0; i < dimWidth; i++)
		{
			imageDataPtr[k][i] = transformPointsPtr[k][i];
		}
	}
	// Free
	for (k = 0; k < imageHeight; k++)
	{
		free(transformPointsPtr[k]);
	}
	free(transformPointsPtr);
	transformPointsPtr = NULL;
}
//==============================================================================
void transformInverse3DImage(dataType ** sourceDataPtr, dataType ** imageDataPtr, Point3D translation, Point3D scaling, Point3D rotation, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType bgValue, dataType centroid[3])
{
	size_t k, i, j;
	// Creates a new Pointer to fill the transformed values
	dataType ** transformPointsPtr = (dataType **)malloc(imageHeight * sizeof(dataType *));

	dataType k_a, i_a, j_a; // Affine indices
						  // Rotation Angles -
	dataType theta = (rotation.y*M_PI) / 180, psi = (rotation.z*M_PI) / 180, phi = (rotation.x*M_PI) / 180;
	// Center Points
	//int hcenter = floor(imageHeight / 2), lcenter = floor(imageLength / 2), wcenter = floor(imageWidth / 2);
	dataType cz = centroid[2], cx = centroid[0], cy = centroid[1];

	// Dimension variable
	// dimWidth = X*Y
	size_t dimWidth = imageLength * imageWidth;
	// Temporary parameters
	dataType tmpX, tmpY, tmpZ, tmp;
	// Transformed
	dataType k_t, i_t, j_t; // Transformed indices
	for (k = 0; k < imageHeight; k++)
	{
		transformPointsPtr[k] = (dataType *)malloc(dimWidth * sizeof(dataType));
		k_a = k - cz; // Move to origin Z
		for (i = 0; i < imageLength; i++)
		{
			i_a = i - cx; // Move to origin x
			for (j = 0; j < imageWidth; j++)
			{
				// 2D to 1D representation for i, j
				size_t x = x_new(i, j, imageLength);

				j_a = j - cy; // Move to origin Y

							  // Apply scaling
				tmpZ = k_a * scaling.z;
				tmpX = i_a * scaling.x;
				tmpY = j_a * scaling.y;

				// Apply Rotation

				// 3. Rotation - No rotation
				i_t = x_rotateInv(tmpZ, tmpX, tmpY, theta, psi, phi);
				j_t = y_rotateInv(tmpZ, tmpX, tmpY, theta, psi, phi);
				k_t = z_rotateInv(tmpZ, tmpX, tmpY, theta, psi, phi);

				// Move back to centroid
				tmpX = i_t + cx;
				tmpY = j_t + cy;
				tmpZ = k_t + cz;

				// Set the values
				i_t = tmpX;
				j_t = tmpY;
				k_t = tmpZ;

				// Add translation
				i_t = i_t + translation.x;
				j_t = j_t + translation.y;
				k_t = k_t + translation.z;

				// Use Interpolation to get the values
				// Locations for Tri-linear Interpolation
				// Z
				int bottom = (int)floor(k_t);
				int top = bottom + 1;
				// X
				int left = (int)floor(i_t);
				int right = left + 1;
				// Y
				int begin = (int)floor(j_t);
				int end = begin + 1;
				// Check if within limits
				if (bottom >= 0 && top < imageHeight && left >= 0 && right < imageLength && begin >= 0 && end < imageWidth)
				{
					tmp = interpolated(k_t, i_t, j_t, top, bottom, left, right, begin, end, sourceDataPtr, imageLength);
					transformPointsPtr[k][x] = tmp;
				}
				else
				{
					transformPointsPtr[k][x] = bgValue; // Background value
				}
			}
		}
	}
	// Copy Transformed back to Image Data Pointer
	for (k = 0; k < imageHeight; k++)
	{
		for (i = 0; i < dimWidth; i++)
		{
			imageDataPtr[k][i] = transformPointsPtr[k][i];
		}
	}
	// Free
	for (k = 0; k < imageHeight; k++)
	{
		free(transformPointsPtr[k]);
	}
	free(transformPointsPtr);
	transformPointsPtr = NULL;
}
//==============================================================================
/*
* Rotated indices
*/
dataType x_rotate(dataType z, dataType x, dataType y, dataType theta, dataType psi)
{
	return (x*cos(psi)*cos(theta) - y * cos(theta)*sin(psi) + z * sin(theta));
}
// Inverse
dataType x_rotateInv(dataType z, dataType x, dataType y, dataType theta, dataType psi, dataType phi)
{
	return (cos(theta)*cos(psi))*(x)+((cos(phi))*sin(psi) + sin(phi)*sin(theta)*cos(psi))*(y)+(sin(phi)*sin(psi) + cos(phi)*(-sin(theta))*cos(psi))*(z);
}
//==============================================================================
dataType y_rotate(dataType z, dataType x, dataType y, dataType theta, dataType psi, dataType phi)
{
	return (((cos(phi))*sin(psi) + sin(phi)*sin(theta)*cos(psi))*(x)+(cos(phi)*cos(psi) + sin(phi)*sin(theta)*(-sin(psi)))*(y)+((-sin(phi))*cos(theta))*(z));
}
//==============================================================================
dataType y_rotateInv(dataType z, dataType x, dataType y, dataType theta, dataType psi, dataType phi)
{
	return ((cos(theta)*(-sin(psi)))*(x)+(cos(phi)*cos(psi) + sin(phi)*sin(theta)*(-sin(psi)))*(y)+((sin(phi))*cos(psi) + cos(phi)*sin(theta)*sin(psi))*(z));
}
//==============================================================================
dataType z_rotate(dataType z, dataType x, dataType y, dataType theta, dataType psi, dataType phi)
{
	return ((sin(phi)*sin(psi) + cos(phi)*(-sin(theta))*cos(psi))*(x)+((sin(phi))*cos(psi) + cos(phi)*sin(theta)*sin(psi))*(y)+(cos(phi)*cos(theta))*(z));
}
//==============================================================================
dataType z_rotateInv(dataType z, dataType x, dataType y, dataType theta, dataType psi, dataType phi)
{
	return ((sin(theta))*(x)+((-sin(phi))*cos(theta))*(y)+(cos(phi)*cos(theta))*(z));
}
//==============================================================================
dataType interpolated(dataType k_t, dataType i_t, dataType j_t, int top, int bottom, int left, int right, int begin, int end, dataType ** imageDataPtr, size_t imageLength)
{
	// 8 Corner Values
	dataType c000 = imageDataPtr[top][x_new(left, begin, imageLength)];
	dataType c001 = imageDataPtr[bottom][x_new(left, begin, imageLength)];
	dataType c010 = imageDataPtr[top][x_new(left, end, imageLength)];
	dataType c011 = imageDataPtr[bottom][x_new(left, end, imageLength)];
	dataType c100 = imageDataPtr[top][x_new(right, begin, imageLength)];
	dataType c101 = imageDataPtr[bottom][x_new(right, begin, imageLength)];
	dataType c110 = imageDataPtr[top][x_new(right, end, imageLength)];
	dataType c111 = imageDataPtr[bottom][x_new(right, end, imageLength)];

	return trilinearInterpolation(i_t, left, right, j_t, begin, end, c000, c001, c010, c011, c100, c101, c110, c111, k_t, bottom, top);
}
//==============================================================================