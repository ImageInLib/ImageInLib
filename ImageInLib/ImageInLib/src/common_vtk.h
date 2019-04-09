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
typedef enum{
	dta_UChar,
	dta_Int,
	dta_UInt,
	dta_Flt,
	dta_Dbl
} vtkDataType;
//==============================================================================
typedef enum {
	copyTo,
	copyFrom
} vtkOperation;
typedef struct {
	double spacing[3]; // distance between neighboring pixels
	double origin[3]; // world coordinate position of the lower left hand corner of the data
	int dimensions[3]; // length, width. height
	vtkDataType vDataType; // datatype
	dataType ** dataPointer; // Pointer with data values
	vtkOperation operation;
} vtkFileInfo;
//==============================================================================
int readVtkFile(const char * inputFilePath, vtkFileInfo * vtkMetaInfo);
int storeVtkFile(const char * outputFilePath, vtkFileInfo * vtkMetaInfo);
//==============================================================================
/* Manual Creation of vtkImageData*/
int createVtkImageData(vtkImageData * imageData, vtkFileInfo * vtkMetaInfo);
//==============================================================================
