#pragma once
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
} VTKHeaderLines;
//==============================================================================
typedef enum {
	BINARY_DATA = 1,
	ASCII_DATA
} loadDataType;
//==============================================================================