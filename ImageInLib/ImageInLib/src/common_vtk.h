#pragma once
//==============================================================================
// ImageInLib Imports
#include "common_functions.h"
#include <stdbool.h>
//==============================================================================
// VTK Imports
// vtkImageData
#include "vtkImageData.h"
//==============================================================================
typedef enum {
	dta_UChar,
	dta_Int,
	dta_UInt,
	dta_Flt,
	dta_Dbl
} VtkDataType;
//==============================================================================
typedef enum {
	copyTo,
	copyFrom
} vtkOperation;
typedef struct {
	double spacing[3]; // distance between neighboring pixels
	double origin[3]; // world coordinate position of the lower left hand corner of the data
	int dimensions[3]; // length, width. height
	VtkDataType vDataType; // datatype
	dataType ** dataPointer; // Pointer with data values
	vtkOperation operation;
} Vtk_File_Info;
//==============================================================================
int readVtkFile(const char * inputFilePath, Vtk_File_Info * vtkMetaInfo);
int storeVtkFile(const char * outputFilePath, Vtk_File_Info * vtkMetaInfo);
//==============================================================================
/* Manual Creation of vtkImageData*/
int createVtkImageData(vtkImageData * imageData, Vtk_File_Info * vtkMetaInfo);
//==============================================================================
