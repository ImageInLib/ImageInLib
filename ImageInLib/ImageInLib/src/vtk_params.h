#ifdef __cplusplus
extern "C" {
#endif

#pragma once
	//==============================================================================
#include <stdbool.h>
//==============================================================================
#define VTK_MAX_HEADER_LINE_LENGTH 256

	typedef char vtk_header_line[VTK_MAX_HEADER_LINE_LENGTH];

	typedef struct
	{
		vtk_header_line line1;
		vtk_header_line line2;
		vtk_header_line line3;
		vtk_header_line line4;
		vtk_header_line line5;
		vtk_header_line line6;
		vtk_header_line line7;
		vtk_header_line line8;
		vtk_header_line line9;
		vtk_header_line line10;
		vtk_header_line line11;
	} VTK_Header_Lines;
	//==============================================================================
	typedef enum {
		BINARY_DATA = 1,
		ASCII_DATA
	} LoadDataType;
	//==============================================================================
	typedef struct {
		bool revertDataBytes;
		bool appendToFile;
	} Storage_Flags;
	//==============================================================================

#ifdef __cplusplus
}
#endif