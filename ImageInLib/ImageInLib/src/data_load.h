#ifdef __cplusplus
extern "C" {
#endif

	/*
	* Author: Markjoe Olunna UBA
	* Purpose: ImageInLife project - 4D Image Segmentation Methods
	* Language:  C
	*/
#pragma once
#include <stdbool.h>
#include <stddef.h>
#include "common_functions.h"
#include "vtk_params.h"

	bool load3dDataArrayD(dataType ** imageDataPtr, const size_t imageLength, const size_t imageWidth,
		const size_t imageHeight, unsigned char * pathPtr);

	bool load3dDataArrayRAW(dataType ** imageDataPtr, const size_t imageLength, const size_t imageWidth,
		const size_t imageHeight, unsigned char * pathPtr, LoadDataType dType);

	bool load3dDataArrayVTK(unsigned char ** imageDataPtr, const size_t imageLength, const size_t imageWidth,
		const size_t imageHeight, unsigned char * pathPtr, VTK_Header_Lines * lines);

#ifdef __cplusplus
}
#endif