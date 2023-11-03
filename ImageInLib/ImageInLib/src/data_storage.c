#pragma warning(disable : 4996)
#pragma warning(disable : 6386)
#pragma warning(disable : 6031)
#pragma warning(disable : 6387)


#include <stdio.h>
#include "data_storage.h"
#include "endianity_bl.h"
#include <stdlib.h>

//function for storage of data to 3D array.
//xDim is the x dimension, yDim is the y dimension and zDim is the z dimension
//pathPtr is pointer to file path
bool store3dDataArrayUC(unsigned char ** array3DPtr, const size_t xDim, const size_t yDim,
	const size_t zDim, unsigned char * pathPtr, bool flagMode)
{
	const size_t dimXY = xDim * yDim;
	size_t k;//loop counter for z dimension
	FILE *cfPtr;

	//checks if the memory was allocated
	if (array3DPtr == NULL)
		return false;

	if (flagMode == true)
	{
		//writing binary data to file
		if ((fopen_s(&cfPtr, pathPtr, "ab")) != 0)
			return false;
	}
	else
	{
		//writing binary data to file
		if ((fopen_s(&cfPtr, pathPtr, "wb")) != 0)
			return false;
	}
	for (k = 0; k < zDim; k++)
	{
		fwrite(array3DPtr[k], sizeof(unsigned char), dimXY, cfPtr);
	}

	fclose(cfPtr);

	return true;
}

//function for initialization of 3D array with some constant value.
//xDim is the x dimension, yDim is the y dimension and zDim is the z dimension
//value is the initial constant value
bool store3dDataVtkUC(unsigned char ** array3DPtr, const size_t xDim, const size_t yDim,
	const size_t zDim, unsigned char * pathPtr, double h)
{
	FILE * outputfile; //file stream
	size_t dimXYZ = xDim * yDim * zDim;
	double sx = h, sy = h, sz = h;
	//checks if the memory was allocated
	if (array3DPtr == NULL)
		return false;

	//checks if the file was sucessfully opened
	if ((fopen_s(&outputfile, pathPtr, "w")) != 0) {
		puts("File could not be opened.");
	}
	else
	{
		fprintf(outputfile, "# vtk DataFile Version 3.0\n");
		fprintf(outputfile, "file in binary format\n");
		fprintf(outputfile, "BINARY\n");
		fprintf(outputfile, "DATASET STRUCTURED_POINTS\n");
		fprintf(outputfile, "DIMENSIONS %zd %zd %zd\n", xDim, yDim, zDim);

		fprintf(outputfile, "ORIGIN %f %f %f\n", (-1.25 + h / 2.), (-1.25 + h / 2.), (-1.25 + h / 2.));
		fprintf(outputfile, "SPACING %f %f %f\n", sx, sy, sz);
		fprintf(outputfile, "POINT_DATA %zd\n", dimXYZ);
		fprintf(outputfile, "SCALARS scalars unsigned_char\n");
		fprintf(outputfile, "LOOKUP_TABLE default\n");
	}
	fclose(outputfile);
	// writing data to vtk file
	store3dDataArrayUC(array3DPtr, xDim, yDim, zDim, pathPtr, true);
	return true;
}

//function for storage of data in 3D binary format.
bool store3dDataArrayD(dataType ** array3DPtr, const size_t xDim, const size_t yDim,
	const size_t zDim, unsigned char * pathPtr, Storage_Flags flags)
{
	size_t i, k;
	const size_t dimXY = xDim * yDim;
	FILE *cfPtr;

	//checks if the memory was allocated
	if (array3DPtr == NULL)
		return false;

	if (flags.appendToFile == true)
	{
		//writing binary data to file
		if ((fopen_s(&cfPtr, pathPtr, "ab")) != 0)
			return false;
	}
	else
	{
		//writing binary data to file
		if ((fopen_s(&cfPtr, pathPtr, "wb")) != 0)
			return false;
	}

	if (flags.revertDataBytes)
	{
		for (k = 0; k < zDim; k++)
		{
			for (i = 0; i < dimXY; i++)
			{
				dataType tmp = array3DPtr[k][i];
				revertBytes(&tmp, sizeof(dataType));
				fwrite(&tmp, sizeof(dataType), dimXY, cfPtr);
			}
		}
	}
	else
	{
		const size_t pointsInSlice = xDim * yDim;

		for (k = 0; k < zDim; k++)
		{
			fwrite(array3DPtr[k], sizeof(dataType), pointsInSlice, cfPtr);
		}
	}

	fclose(cfPtr);

	return true;
}

//function for storage of data in 3D ASCII format.
bool store3dDataArrayASCII(dataType ** array3DPtr, const size_t xDim, const size_t yDim,
	const size_t zDim, unsigned char * pathPtr, Storage_Flags flags)
{
	const size_t dim2D = xDim * yDim;
	FILE *cfPtr;

	//checks if the memory was allocated
	if (array3DPtr == NULL)
		return false;

	if (flags.appendToFile == true)
	{
		//writing binary data to file
		if ((fopen_s(&cfPtr, pathPtr, "a")) != 0)
			return false;
	}
	else
	{
		//writing binary data to file
		if ((fopen_s(&cfPtr, pathPtr, "w")) != 0)
			return false;
	}

	if (flags.revertDataBytes)
	{
		for (size_t k = 0; k < zDim; k++)
		{
			for (size_t j = 0; j < dim2D; j++)
			{
				double tmp = array3DPtr[k][j];
				revertBytes(&tmp, sizeof(double));
				fprintf(cfPtr, "%.16lf \n", tmp);
			}
		}
	}
	else
	{
		for (size_t k = 0; k < zDim; k++)
		{
			for (size_t i = 0; i < dim2D; i++)
			{
				fprintf(cfPtr, "%.16lf \n", array3DPtr[k][i]);
			}
		}
	}

	fclose(cfPtr);

	return true;
}

//xDim is the x dimension, yDim is the y dimension and zDim is the z dimension
//value is the initial constant value
bool store3dDataVtkD(dataType ** array3DPtr, const size_t xDim, const size_t yDim,
	const size_t zDim, unsigned char * pathPtr, double h, Storage_Flags flags)
{
	FILE * outputfile; //file stream
	size_t dimXYZ = xDim * yDim * zDim;
	double sx = h, sy = h, sz = h;
	//checks if the memory was allocated
	if (array3DPtr == NULL)
		return false;

	//checks if the file was sucessfully opened
	if ((fopen_s(&outputfile, pathPtr, "w")) != 0) {
		return false;
	}
	else
	{
		fprintf(outputfile, "# vtk DataFile Version 3.0\n");
		fprintf(outputfile, "file in binary format\n");
		fprintf(outputfile, "BINARY\n");
		fprintf(outputfile, "DATASET STRUCTURED_POINTS\n");
		fprintf(outputfile, "DIMENSIONS %zd %zd %zd\n", xDim, yDim, zDim);

		//fprintf(outputfile, "ORIGIN %f %f %f\n", (-1.25 + h / 2.), (-1.25 + h / 2.), (-1.25 + h / 2.));
		fprintf(outputfile, "ORIGIN %f %f %f\n", (double)0, (double)0, (double)0);
		fprintf(outputfile, "SPACING %f %f %f\n", sx, sy, sz);
		fprintf(outputfile, "POINT_DATA %zd\n", dimXYZ);
		fprintf(outputfile, "SCALARS scalars double\n");
		fprintf(outputfile, "LOOKUP_TABLE default\n");
	}
	fclose(outputfile);
	// writing data to vtk file
	store3dDataArrayD(array3DPtr, xDim, yDim, zDim, pathPtr, flags);
	return true;
}

bool store3dRealDataVtkD(dataType ** array3DPtr, const size_t imageLength, const size_t imageWidth,
	const size_t imageHeight, unsigned char * pathPtr, VTK_Header_Lines * lines, Storage_Flags flags)
{
	//checks if the memory was allocated
	if (array3DPtr == NULL)
		return false;

	// imageLength == xDim, imageWidth == yDim, imageHeight == zDim
	FILE * outputfile; //file stream

					   //checks if the file was sucessfully opened
	if ((fopen_s(&outputfile, pathPtr, "w")) != 0) {
		return false;
	}
	else
	{
		//Write header
		fprintf(outputfile, "%s", lines->line1);
		fprintf(outputfile, "%s", lines->line2);
		fprintf(outputfile, "%s", lines->line3);
		fprintf(outputfile, "%s", lines->line4);
		fprintf(outputfile, "%s", lines->line5);
		fprintf(outputfile, "%s", lines->line6);
		fprintf(outputfile, "%s", lines->line7);
		fprintf(outputfile, "%s", lines->line8);
		fprintf(outputfile, "SCALARS scalars double\n");
		fprintf(outputfile, "%s", lines->line10);
		//fprintf(outputfile, "SCALARS scalars double\n");
		//fprintf(outputfile, "LOOKUP_TABLE default\n");
	}
	fclose(outputfile);
	// writing data to vtk file
	store3dDataArrayD(array3DPtr, imageLength, imageWidth, imageHeight, pathPtr, flags);
	return true;
}

bool store3dRealDataVtkUC(unsigned char ** array3DPtr, const size_t imageLength, const size_t imageWidth,
	const size_t imageHeight, unsigned char * pathPtr, VTK_Header_Lines * lines)
{
	//checks if the memory was allocated
	if (array3DPtr == NULL)
		return false;

	// imageLength == xDim, imageWidth == yDim, imageHeight == zDim
	FILE * outputfile; //file stream

					   //checks if the file was sucessfully opened
	if ((fopen_s(&outputfile, pathPtr, "w")) != 0) {
		return false;
	}
	else
	{
		//Write header
		fprintf(outputfile, "%s", lines->line1);
		fprintf(outputfile, "%s", lines->line2);
		fprintf(outputfile, "%s", lines->line3);
		fprintf(outputfile, "%s", lines->line4);
		fprintf(outputfile, "%s", lines->line5);
		fprintf(outputfile, "%s", lines->line6);
		fprintf(outputfile, "%s", lines->line7);
		fprintf(outputfile, "%s", lines->line8);
		fprintf(outputfile, "%s", lines->line9);
		fprintf(outputfile, "%s", lines->line10);
	}
	fclose(outputfile);
	// writing data to vtk file
	store3dDataArrayUC(array3DPtr, imageLength, imageWidth, imageHeight, pathPtr, true);
	return true;
}

//==================================
//function for storage of data in 2D PGM. Used format is defined by a flag writeRawData (true = raw, false = ascii).
bool store2dPGM(dataType* imageDataPtr, const size_t xDim, const size_t yDim, const char* pathPtr, const bool writeRawData)
{
	FILE* pgmimg;
	pgmimg = fopen(pathPtr, "w");

	if (writeRawData) //pgm format
	{
		fprintf(pgmimg, "P5\n");
	}
	else
	{
		fprintf(pgmimg, "P2\n");
	}

	// Writing Width and Height
	fprintf(pgmimg, "%zu %zu\n", xDim, yDim);

	// Writing the maximum gray value
	fprintf(pgmimg, "255\n");

	const int dataSize = (int)(xDim * yDim);

	if (writeRawData)
	{
		unsigned char * rawData = (unsigned char *)malloc(dataSize);
		
		if (rawData) {
			for (size_t i = 0; i < dataSize; i++) {
				rawData[i] = (unsigned char)imageDataPtr[i];
			}
		}

		fwrite(rawData, sizeof(unsigned char), dataSize, pgmimg);
		free(rawData);
	}
	else
	{
		for (size_t i = 0; i < dataSize; i++) {
			// Writing the gray values in the 2D array to the file
			fprintf(pgmimg, "%d \n", (unsigned char)imageDataPtr[i]);
		}
	}

	fclose(pgmimg);
	return true;
}

//==================================
bool storeCSV(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const char* pathPtr)
{
	FILE* pgmimg;
	pgmimg = fopen(pathPtr, "w");

    for (size_t i = 0; i < xDim; i++) {
        for (size_t j = 0; j < yDim; j++) {
            // Writing the gray values in the 2D array to the file
            fprintf(pgmimg, "%f", imageDataPtr[i][j]);

			if (j == yDim - 1)
			{
				fprintf(pgmimg, "\n");
			}
			else
			{
				fprintf(pgmimg, ",");
			}
        }
    }

	fclose(pgmimg);
	return true;
}

bool loadCSV(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const char* pathPtr)
{
	FILE* pgmimg;
	pgmimg = fopen(pathPtr, "r");

	for (size_t i = 0; i < xDim; i++) {
		for (size_t j = 0; j < yDim; j++) {
			// Writing the gray values in the 2D array to the file
			fscanf(pgmimg, "%f", &(imageDataPtr[i][j]));

			if (j == yDim - 1)
			{
				fscanf(pgmimg, "\n");
			}
			else
			{
				fscanf(pgmimg, ",");
			}
		}
	}

	fclose(pgmimg);
	return true;
}

//==================================
bool store2dRawData(dataType* array2DPtr, const size_t xDim, const size_t yDim, const char * pathPtr, Storage_Flags flags) {

	size_t i;

	//checks if the memory was allocated
	if (array2DPtr == NULL)
		return false;

	FILE* cfPtr;

	if (flags.appendToFile == true) {
		//writing binary data to file
		if ((fopen_s(&cfPtr, pathPtr, "ab")) != 0)
			return false;
	}
	else
	{
		//writing binary data to file
		if ((fopen_s(&cfPtr, pathPtr, "wb")) != 0)
			return false;
	}

	if (flags.revertDataBytes == true) {
		for (i = 0; i < xDim * yDim; i++) {
			dataType tmp = array2DPtr[i];
			revertBytes(&tmp, sizeof(dataType));
			fwrite(&tmp, sizeof(dataType), 1, cfPtr);
		}
	}
	else {
		fwrite(array2DPtr, sizeof(dataType), xDim * yDim, cfPtr);
	}

	return true;
}