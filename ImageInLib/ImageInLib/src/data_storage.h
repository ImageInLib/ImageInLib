#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdbool.h>
#include <stddef.h>
#include "common_functions.h"
#include "vtk_params.h"

	bool store3dDataArrayUC(unsigned char ** array3DPtr, const size_t xDim, const size_t yDim,
		const size_t zDim, unsigned char * pathPtr, bool flagmode);

	bool store3dDataVtkUC(unsigned char ** array3DPtr, const size_t xDim, const size_t yDim,
		const size_t zDim, unsigned char * pathPtr, double h);

	//Stores data in 3D binary format.
	//xDim - x dimension
	//yDim - y dimension
	//zDim - z dimension
	//pathPtr - pointer to file path
	//appendToFile - flag defining whether data will be appended to file
	//revertDataBytes - flag defining whether output data bytes will be reversed
	bool store3dDataArrayD(dataType ** array3DPtr, const size_t xDim, const size_t yDim,
		const size_t zDim, unsigned char * pathPtr, Storage_Flags flags);

	bool store3dDataVtkD(dataType ** array3DPtr, const size_t xDim, const size_t yDim,
		const size_t zDim, unsigned char * pathPtr, double h, Storage_Flags flags);

	bool store3dDataArrayASCII(dataType ** array3DPtr, const size_t xDim, const size_t yDim,
		const size_t zDim, unsigned char * pathPtr, Storage_Flags flags);

	bool store3dRealDataVtkD(dataType ** array3DPtr, const size_t xDim, const size_t yDim,
		const size_t zDim, unsigned char * pathPtr, VTK_Header_Lines *lines, Storage_Flags flags);

	bool store3dRealDataVtkUC(unsigned char ** array3DPtr, const size_t imageLength, const size_t imageWidth,
		const size_t imageHeight, unsigned char * pathPtr, VTK_Header_Lines *lines);

	bool store2dPGM(dataType* imageDataPtr, const size_t xDim, const size_t yDim, const char* pathPtr, const bool writeRawData);

	bool store2dCSV(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const char* pathPtr);

	bool store2dRawData(dataType* array2DPtr, const size_t xDim, const size_t yDim, const char* pathPtr, Storage_Flags flags);
#ifdef __cplusplus
}
#endif