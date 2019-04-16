#include <stdio.h>
#include "data_storage.h"
#include "endianity_bl.h"

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
bool store3dDataArrayD(double ** array3DPtr, const size_t xDim, const size_t yDim,
	const size_t zDim, unsigned char * pathPtr, Storage_Flags flags)
{
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
		for (size_t k = 0; k < zDim; k++)
		{
			for (size_t j = 0; j < dimXY; j++)
			{
				double tmp = array3DPtr[k][j];
				revertBytes(&tmp, sizeof(double));
				fwrite(&tmp, sizeof(double), 1, cfPtr);
			}
		}
	}
	else
	{
		for (size_t k = 0; k < zDim; k++)
		{
			fwrite(array3DPtr[k], sizeof(double), dimXY, cfPtr);
		}
	}

	fclose(cfPtr);

	return true;
}

//function for storage of data in 3D ASCII format.
bool store3dDataArrayASCII(double ** array3DPtr, const size_t xDim, const size_t yDim,
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
bool store3dDataVtkD(double ** array3DPtr, const size_t xDim, const size_t yDim,
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

		fprintf(outputfile, "ORIGIN %f %f %f\n", (-1.25 + h / 2.), (-1.25 + h / 2.), (-1.25 + h / 2.));
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

bool store3dRealDataVtkD(double ** array3DPtr, const size_t imageLength, const size_t imageWidth,
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