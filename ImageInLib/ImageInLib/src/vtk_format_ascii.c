#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "vtk_format_ascii.h"

//function for initialization of 3D array with some constant value.
//xDim is the x dimension, yDim is the y dimension and zDim is the z dimension
//value is the initial constant value
bool vtkformat_3D(unsigned char ** array3DPtr, const size_t xDim, const size_t yDim,
	const size_t zDim, unsigned char * pathPtr)
{
	size_t k;//loop counter for z dimension
	FILE * outputfile; //file stream
	size_t dim2D = xDim * yDim;
	size_t dim3D = xDim * yDim * zDim;
	double sx = 1, sy = 1, sz = 1;
    //checks if the memory was allocated
	if (array3DPtr == NULL)
		return false;

	//checks if the file was sucessfully opened
	if (fopen_s(&outputfile, pathPtr, "w") != 0) {
		puts("File could not be opened.");
	}
	else
	{
		fprintf(outputfile, "# vtk DataFile Version 3.0\n");
		fprintf(outputfile, "file in ascii format\n");
		fprintf(outputfile, "ASCII\n");
		fprintf(outputfile, "DATASET STRUCTURED_POINTS\n");
		fprintf(outputfile, "DIMENSIONS %zd %zd %zd\n", xDim, yDim, zDim);

		fprintf(outputfile, "ORIGIN 0 0  0\n");
		fprintf(outputfile, "SPACING %f %f %f\n", sx, sy, sz);
		fprintf(outputfile, "POINT_DATA %zd\n", dim3D);
		fprintf(outputfile, "SCALARS scalars float\n");
		fprintf(outputfile, "LOOKUP_TABLE default\n");	
	}
	fclose(outputfile);
	// writing data to vtk file
	for (k = 0; k < zDim; k++)
	{
		vtkformat_2D(array3DPtr[k], yDim, dim2D, pathPtr);
	}
	return true;
}

//function for initialization of 2D array with some constant value.
bool vtkformat_2D(unsigned char * array2DPtr, const size_t yDim, size_t dim2D,
	              unsigned char * pathPtr)
{
	size_t i;
	FILE * outputfile; //file stream

	//checks if the memory was allocated
	if (array2DPtr == NULL)
		return false;

	//checks if the file was sucessfully opened
	if (fopen_s(&outputfile, pathPtr, "a") != 0) {
		puts("File could not be opened.");
	}
	else
	{
		// writing data to vtk file
		for (i = 0; i < dim2D; i++) {
			fprintf(outputfile, "%f ", (float)array2DPtr[i]);
			//putc(array2DPtr[i], outputfile);
			/*if (scale == 0)
				fprintf(outputfile, "%f ", (float)array2DPtr[i]);
			else
				fprintf(outputfile, "%f ", (float)(scale*array2DPtr[i]));*/
			
			if (((i + 1) % (yDim)) == 0)
				fprintf(outputfile, "\n");
		}
		fclose(outputfile);
	}

	return true;
}