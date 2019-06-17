#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdio.h> // Standard lib for input and output functions
#include "common_functions.h"
#include "../src/vtk_params.h"
#include <stdbool.h>

	typedef enum {
		LOAD_DATA_VTK = 1,
		LOAD_DATA_RAW,
		STORE_DATA_VTK,
		STORE_DATA_RAW
	} OperationType;

	bool manageFile(void  ** imageDataPtr, const size_t length, const size_t width,
		const size_t height, unsigned char * pathPtr, OperationType operation, LoadDataType dType, Storage_Flags flags);

	/*Converts to doulbe*/
	void convertTodataType(unsigned char ** dataPtrUC, dataType ** dataPtrD, const size_t dimXY, const size_t height);

#ifdef __cplusplus
}
#endif