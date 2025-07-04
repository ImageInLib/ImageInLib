#pragma warning(disable : 4996)

#include "common_functions.h"
#include "data_load.h"
#include "endianity_bl.h"
#include <stdio.h>
#include <string.h>
#include "endianity_bl.h"
#include <stdlib.h>
#include "imageInterpolation.h"


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

	const size_t pointsInSlice = imageLength * imageWidth;

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
		if (dType == BINARY_DATA) {
			for (k = 0; k < imageHeight; k++) {
				fread(imageDataPtr[k], sizeof(dataType), pointsInSlice, file);
			}
		}
		else {
			if (dType == ASCII_DATA) {
				for (k = 0; k < imageHeight; k++)
				{
					for (i = 0; i < imageLength; i++)
					{
						for (j = 0; j < imageWidth; j++)
						{
							// 2D to 1D representation for i, j
							xd = x_new(i, j, imageLength);
							fscanf_s(file, "%f", &value);
							imageDataPtr[k][xd] = (dataType)value;
						}
					}
				}
			}
		}
	}
  
	fclose(file);
	return true;
}

//==================================
//Load 2D .pgm (ascii) image
bool load2dPGM(dataType* imageDataPtr, const size_t xDim, const size_t yDim, const char* pathPtr)
{
    int intensity;
    size_t i, j;

    char line1[4];
	const int dataSize = (int)(xDim * yDim);
	unsigned char* line2 = malloc(dataSize);

    FILE* file;
    if (fopen_s(&file, pathPtr, "r") != 0) {
        printf("File not found");
        return false;
    }

    int pgmVersion;
    fgets(line1, 4, file);
	if (sscanf(line1, "P%d", &pgmVersion) != 1) {
		fclose(file);
		free(line2);
		return false; // Error reading PGM version
	}

	//filtering out potential comment
	do {
		fgets(line2, dataSize, file);
	} while (line2[0] == '#');

	size_t tmpX, tmpY;
	if(sscanf(line2, "%zu %zu", &tmpX, &tmpY) != 2) {
		fclose(file);
		free(line2);
		return false; // Error reading dimensions
	}

	if (xDim != tmpX || yDim != tmpY) {
		//dimensions of used array is not compatible with loaded image
		fclose(file);
		return false;
	}

	fgets(line2, dataSize, file);

    if (pgmVersion == 2) //ascii
	{

		size_t xd = 0;
        for (i = 0; i < yDim; i++) {
            for (j = 0; j < xDim; j++) {
				if (fscanf(file, "%d", &intensity) != 1) {
					fclose(file);
					return false; // Error reading intensity
				}
				// 2D to 1D representation for i, j
				xd = x_new(j, i, xDim);
                imageDataPtr[xd] = (dataType)intensity;
            }
        }
    }     
	if (pgmVersion == 5) //raw
	{
		fgets(line2, dataSize, file);
		for (i = 0; i < dataSize; i++) {
				imageDataPtr[i] = (dataType)line2[i];
		}
	}
    else
    {
		free(line2);
        fclose(file);
        return false;
    }

	free(line2);
    fclose(file);
    return true;
}

//==================================
//Load 2D .raw image

bool load2dArrayRAW(dataType* imageDataPtr, const size_t length, const size_t width, const char* pathPtr, LoadDataType dType){
	
	//checks if the memory was allocated
	if (imageDataPtr == NULL || pathPtr == NULL) {
		printf("Not initialized data \n");
		return false;
	}

	FILE* file;
	char rmode[4];
	const size_t pointsInSlice = length * width;
	dataType value = 0.0;

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
		if (dType == BINARY_DATA) {
			fread(imageDataPtr, sizeof(dataType), pointsInSlice, file);
		}
		else {
			if (dType == ASCII_DATA) {
				for (size_t i = 0; i < length; i++) {
					for (size_t j = 0; j < width; j++) {
						fscanf_s(file, "%f", &value);
						imageDataPtr[x_new(i, j, length)] = (dataType)value;
					}
				}			
			}
		}
	}

	fclose(file);
	return true;
}

//==================================

bool loadListof3dPoints(Image_Data image, Curve3D* pCurve, const char* filePath, CoordinateSystem cSystem)
{
	
	FILE* file;
	if (fopen_s(&file, filePath, "r") != 0) {
		printf("Enable to open");
		return false;
	}

	Point3D previous_point;
	double dist = 0;
	
	dataType x = 0, y = 0, z = 0;

	for (size_t i = 0; i < pCurve->numPoints; i++) {
		
		fscanf_s(file, "%f", &x);
		fscanf_s(file, ",");
		fscanf_s(file, "%f", &y);
		fscanf_s(file, ",");
		fscanf_s(file, "%f", &z);
		fscanf_s(file, "\n");

		Point3D current_point = { x , y , z };

		if (cSystem == REAL) {
			//Get image coordinate
			current_point = getImageCoordFromRealCoord3D(current_point, image.origin, image.spacing, image.orientation);
			pCurve->pPoints[i].x = current_point.x;
			pCurve->pPoints[i].y = current_point.y;
			pCurve->pPoints[i].z = current_point.z;
		}
		else {
			pCurve->pPoints[i].x = x;
			pCurve->pPoints[i].y = y;
			pCurve->pPoints[i].z = z;
		}

		previous_point = current_point;

	}
	fclose(file);

	return true;
}