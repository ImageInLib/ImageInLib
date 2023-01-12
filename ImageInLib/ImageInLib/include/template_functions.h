#pragma once

#include<iostream>
#include "common_functions.h"
#include "../src/data_load.h"
#include "../src/endianity_bl.h"
#include <stdio.h>
#include <string.h>

template <typename T>
bool store3dRawData(T** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, const char* pathPtr);

template<typename T>
inline bool store3dRawData(T** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, const char* pathPtr)
{
	size_t i, j, k;
	FILE* file;

	//checks if the memory was allocated
	if (imageDataPtr == NULL)
		return false;

	//writing binary data to file
	if (fopen_s(&file, pathPtr, "wb") != 0)
		return false;

	for (k = 0; k < zDim; k++)
	{
		for (i = 0; i < xDim; i++)
		{
			for (j = 0; j < yDim; j++)
			{
				fwrite(&imageDataPtr[k][x_new(i, j, xDim)], sizeof(T), 1, file);
			}
		}
	}

	fclose(file);

	return true;
}

template <typename T>
bool load3dArrayRAW(T** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, const char* pathPtr);

template <typename T>
inline bool load3dArrayRAW(T** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, const char* pathPtr)
{
	dataType firstCpuTime, secondCpuTime;
	firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
	//checks if the memory was allocated
	if (imageDataPtr == NULL)
		return false;

	if (pathPtr == NULL)
		return false;

	size_t i, j, k, xd;
	T value;
	FILE* file;
	char rmode[4] = "rb";

	// Reading data file
	if (fopen_s(&file, pathPtr, rmode) != 0) {
		fprintf(stderr, "Error: Unable to open file %s\n\n", pathPtr);
		return false;
	}
	else {
		for (k = 0; k < zDim; k++){
			fread(imageDataPtr[k], sizeof(T), xDim*yDim, file);
		}
	}
	//----
	//Endianness check

	//if (*(char*)&imageDataPtr[0][0] == imageDataPtr[0][0]) {
	//	printf("little endian before swaping \n");
	//}
	//else {
	//	printf("big endian before swaping \n");
	//}

	//revert byte
	for (k = 0; k < zDim; k++) {
		for (i = 0; i < xDim * yDim; i++) {
			revertBytes(&imageDataPtr[k][i], sizeof(T));
		}
	}

	////Endianness check
	//if (*(char*)&imageDataPtr[0][0] == imageDataPtr[0][0]) {
	//	printf("little endian after swaping \n");
	//}
	//else {
	//	printf("big endian after swaping \n");
	//}

	secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
	printf("CPU time: %e secs\n", secondCpuTime - firstCpuTime);

	fclose(file);
	return true;
}
