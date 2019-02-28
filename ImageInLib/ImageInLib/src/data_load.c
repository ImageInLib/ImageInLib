/*
* Author: Markjoe Olunna UBA
* Purpose: ImageInLife project - 4D Image Segmentation Methods
* Language:  C
*/
#include <stdio.h>
#include "data_load.h"
#include "common_functions.h"

bool load3dDataArrayVTK(unsigned char ** imageDataPtr, const size_t imageLength, const size_t imageWidth,
	const size_t imageHeight, unsigned char * pathPtr, VTKHeaderLines * lines)
{
	//checks if the memory was allocated
	if (imageDataPtr == NULL)
		return false;

	if (pathPtr == NULL)
		return false;

	// imageLength == xDim, imageWidth == yDim, imageHeight == zDim
	size_t i, j, k;
	size_t x; //x = x_new(i, j, length);
	int value;
	FILE *file;

	//Reading data from file
	if (fopen_s(&file, pathPtr, "rb") != 0) {
		return false;
	}
	else {
		//Read header
		fgets(lines->line1, 100, file);
		fgets(lines->line2, 100, file);
		fgets(lines->line3, 100, file);
		fgets(lines->line4, 100, file);
		fgets(lines->line5, 100, file);
		fgets(lines->line6, 100, file);
		fgets(lines->line7, 100, file);
		fgets(lines->line8, 100, file);
		fgets(lines->line9, 100, file);
		fgets(lines->line10, 100, file);

		//Read other data
		for (k = 0; k < imageHeight; k++)
		{
			for (i = 0; i < imageLength; i++)
			{
				for (j = 0; j < imageWidth; j++)
				{
					// 2D to 1D representation for i, j
					x = x_new(i, j, imageLength);

					//fscanf(file, "%d", &value);
					value = getc(file);
					imageDataPtr[k][x] = value;//  / 255.
				}
			}
		}
	}
	fclose(file);
	return true;
}

bool load3dDataArrayD(dataType ** imageDataPtr, const size_t imageLength, const size_t imageWidth,
	const size_t imageHeight, unsigned char * pathPtr)
{
	//checks if the memory was allocated
	if (imageDataPtr == NULL)
		return false;
	if (pathPtr == NULL)
		return false;

	size_t k;
	const size_t dim2D = imageWidth * imageLength;
	FILE *file;

	//writing binary data to file
	if (fopen_s(&file, pathPtr, "rb") != 0) {
		return false;
	}
	else {
		for (k = 0; k < imageHeight; k++)
		{
			fread(imageDataPtr[k], sizeof(double), dim2D, file);
		}
	}
	fclose(file);
	return true;
}

bool load3dDataArrayA(dataType ** imageDataPtr, const size_t imageLength, const size_t imageWidth,
	const size_t imageHeight, unsigned char * pathPtr)
{
	//checks if the memory was allocated
	if (imageDataPtr == NULL)
		return false;

	if (pathPtr == NULL)
		return false;

	size_t i, j, k;
	size_t x; //x = x_new(i, j, length);
	int value;
	FILE *file;

	//writing binary data to file
	if (fopen_s(&file, pathPtr, "r") != 0) {
		return false;
	}
	else {
		for (k = 0; k < imageHeight; k++)
		{
			for (i = 0; i < imageLength; i++)
			{
				for (j = 0; j < imageWidth; j++)
				{
					// 2D to 1D representation for i, j
					x = x_new(i, j, imageLength);

					// Old
					//fscanf(file, "%d", &value);
					// New
					fscanf_s(file, "%d", &value, sizeof(&value));
					imageDataPtr[k][x] = value / 255.;
				}
			}
		}
	}
	fclose(file);
	return true;
}