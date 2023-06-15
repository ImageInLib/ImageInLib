#pragma warning(disable : 4996)
#pragma warning(disable : 6386)
#pragma warning(disable : 6031)
#pragma warning(disable : 6387)

/*
* Author: Markjoe Olunna UBA
* Purpose: ImageInLife project - 4D Image Segmentation Methods
* Language:  C
*/
#include "common_functions.h"
#include "data_load.h"
#include "endianity_bl.h"
#include <stdio.h>
#include <string.h>
#include "endianity_bl.h"

bool load3dDataArrayVTK(unsigned char ** imageDataPtr, const size_t imageLength, const size_t imageWidth,
	const size_t imageHeight, unsigned char * pathPtr, VTK_Header_Lines * lines)
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

bool load3dDataArrayRAW(dataType ** imageDataPtr, const size_t imageLength, const size_t imageWidth,
	const size_t imageHeight, unsigned char * pathPtr, LoadDataType dType)
{
	//checks if the memory was allocated
	if (imageDataPtr == NULL)
		return false;

	if (pathPtr == NULL)
		return false;

	size_t i, j, k, xd;
	dataType value;
	FILE *file;
	char rmode[4];

	if (dType == BINARY_DATA)
	{
		strcpy_s(rmode, sizeof(rmode), "rb");
	}
	else if (dType == ASCII_DATA)
	{
		strcpy_s(rmode, sizeof(rmode), "r");
	}
	else
	{
		return false; // Unindentified read mode, exiting the function
	}
	// Reading data file
	if (fopen_s(&file, pathPtr, rmode) != 0) {
		fprintf(stderr, "Error: Unable to open file %s\n\n", pathPtr);
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
					xd = x_new(i, j, imageLength);

					if (dType == BINARY_DATA)
					{
						fread(&value, sizeof(dataType), 1, file);
						imageDataPtr[k][xd] = (dataType)value;
					}
					else if (dType == ASCII_DATA)
					{
						fscanf_s(file, "%f", &value);
						imageDataPtr[k][xd] = (dataType)value;
					}
				}
			}
		}
	}
  
	//change from little endian to big endian
	for (k = 0; k < imageHeight; k++)
	{
		for (i = 0; i < imageLength * imageWidth; i++)
		{
			revertBytes(&imageDataPtr[k][i], sizeof(dataType));
		}
	}
  
	fclose(file);
	return true;
}

//==================================
//Load 2D .pgm (ascii) image
bool load2dPGM(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const char* pathPtr)
{
	int intensity;
	size_t i, j;

	char line1[5];
	char line2[80];

	FILE* file;
	if (fopen_s(&file, pathPtr, "r") != 0) {
		printf("File not found");
		return false;
	}

	fgets(line1, 10, file);

	do {
		fgets(line2, 80, file);
	} while (line2[0] == '#');

	sscanf(line2, "%d %d", &xDim, &yDim);

	fgets(line2, 10, file);

	for (i = 0; i < xDim; i++) {
		for (j = 0;j < yDim; j++) {
			fscanf(file, "%d", &intensity);
			imageDataPtr[i][j] = (float)intensity;
		}
	}
	
	fclose(file);

	return true;
}