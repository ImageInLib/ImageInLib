#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"
#include "vtk_params.h"

	bool load3dDataArrayD(dataType ** imageDataPtr, const size_t imageLength, const size_t imageWidth,
		const size_t imageHeight, unsigned char * pathPtr);

	bool load3dDataArrayRAW(dataType ** imageDataPtr, const size_t imageLength, const size_t imageWidth,
		const size_t imageHeight, unsigned char * pathPtr, LoadDataType dType);

	bool load3dDataArrayVTK(unsigned char ** imageDataPtr, const size_t imageLength, const size_t imageWidth,
		const size_t imageHeight, unsigned char * pathPtr, VTK_Header_Lines * lines);

	//==================================
	//Load 2D .pgm (ascii) image
	bool load2dPGM(dataType* imageDataPtr, const size_t xDim, const size_t yDim, const char* pathPtr);

	//==================================

	bool load2dArrayRAW(dataType* imageDataPtr, const size_t length, const size_t width, const char* pathPtr, LoadDataType dType);

	//==================================

	typedef enum {
		REAL = 1,
		IMAGE
	} CoordinateSystem;

	/// <summary>
	/// Function to load a list of 3D points in real or image coordinate system
	/// </summary>
	/// <param name="image">Stuture holding the image</param>
	/// <param name="pCurve">curve to hold the 3D points</param>
	/// <param name="filePath">path to the file</param>
	/// <param name="cSystem">coordinate system</param>
	/// <returns>return true the loading succed, false otherwise</returns>
	bool loadListof3dPoints(Image_Data image, Curve3D* pCurve, const char* filePath, CoordinateSystem cSystem);

#ifdef __cplusplus
}
#endif