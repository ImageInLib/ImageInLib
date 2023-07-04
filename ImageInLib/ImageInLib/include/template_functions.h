/*
* Author: Konan ALLALY
* Purpose: INFLANET project - Image Processing in Nuclear Medicine (2D/3D)
* Language:  C and C++
*/
#pragma once

#include<iostream>
#include "common_functions.h"
#include "../src/data_load.h"
#include "../src/endianity_bl.h"
#include <stdio.h>
#include <string.h>

typedef enum {
	LOAD_DATA = 1,
	STORE_DATA
} Operation;

//3D Image

template <typename T>
bool load3dArrayRAW(T** imageDataPtr, const size_t length, const size_t width, const size_t height, const char* pathPtr, bool revert);

template <typename T>
inline bool load3dArrayRAW(T** imageDataPtr, const size_t length, const size_t width, const size_t height, const char* pathPtr, bool revert)
{
	//checks if the memory was allocated
	if (imageDataPtr == NULL || pathPtr == NULL) {
		printf("Not initialized data \n");
		return false;
	}

	size_t i, k;
	FILE* file;
	char rmode[4] = "rb";
	const size_t pointsInSlice = width * length;

	// Reading data file
	if (fopen_s(&file, pathPtr, rmode) != 0) {
		fprintf(stderr, "Error: Unable to open file %s\n\n", pathPtr);
		return false;
	}
	else {
		for (k = 0; k < height; k++) {
			fread(imageDataPtr[k], sizeof(T), pointsInSlice, file);
		}
	}

	const size_t dataSize = sizeof(T);
	char* swapBuf = (char*)malloc(dataSize);

	//revert byte
	if (revert == true) {
		for (k = 0; k < height; k++) {
			for (i = 0; i < pointsInSlice; i++) {
				revertBytesEx(&imageDataPtr[k][i], dataSize, swapBuf);
			}
		}
	}

	free(swapBuf);

	fclose(file);
	return true;
}

template <typename T>
bool store3dRawData(T** imageDataPtr, const size_t length, const size_t width, const size_t height, const char* pathPtr);

template<typename T>
inline bool store3dRawData(T** imageDataPtr, const size_t length, const size_t width, const size_t height, const char* pathPtr)
{
	size_t k;
	FILE* file;

	//checks if the memory was allocated
	if (imageDataPtr == NULL)
		return false;

	//writing binary data to file
	if (fopen_s(&file, pathPtr, "wb") != 0)
		return false;

	const size_t pointsInSlice = length * width;

	for (k = 0; k < height; k++)
	{
		fwrite(imageDataPtr[k], sizeof(T), pointsInSlice, file);
	}

	fclose(file);

	return true;
}

template<typename T>
bool manageRAWFile3D(T** imageDataPtr, const size_t length, const size_t width, const size_t height, const char* pathPtr, Operation operation, bool revert);

template <typename T>
inline bool manageRAWFile3D(T** imageDataPtr, const size_t length, const size_t width, const size_t height, const char* pathPtr, Operation operation, bool revert) {

	bool status = false; // Initial Status, only changed to true if the operation is successful

	switch (operation)
	{
	case LOAD_DATA:
		status = load3dArrayRAW<T>(imageDataPtr, length, width, height, pathPtr, revert);
		break;
	case STORE_DATA:
		status = store3dRawData<T>(imageDataPtr, length, width, height, pathPtr);
		break;
	default:
		break;
	}

	return status;
}

//===========================================================================

//2D image

template <typename T>
bool load2dArrayRAW(T* imageDataPtr, const size_t length, const size_t width, const char* pathPtr, bool revert);

template <typename T>
inline bool load2dArrayRAW(T* imageDataPtr, const size_t length, const size_t width, const char* pathPtr, bool revert)
{
	//checks if the memory was allocated
	if (imageDataPtr == NULL || pathPtr == NULL) {
		printf("Not initialized data \n");
		return false;
	}

	size_t i;
	FILE* file;
	char rmode[4] = "rb";
	const size_t pointsInSlice = length * width;

	// Reading data file
	if (fopen_s(&file, pathPtr, rmode) != 0) {
		fprintf(stderr, "Error: Unable to open file %s\n\n", pathPtr);
		return false;
	}
	else {
		fread(imageDataPtr, sizeof(T), pointsInSlice, file);
	}

	const size_t dataSize = sizeof(T);
	char* swapBuf = (char*)malloc(dataSize);

	//revert byte
	if (revert == true) {
		for (i = 0; i < pointsInSlice; i++) {
			revertBytesEx(&imageDataPtr[i], dataSize, swapBuf);
		}
	}
	
	free(swapBuf);

	fclose(file);
	return true;
}

template <typename T>
bool store2dRawData(T* imageDataPtr, const size_t length, const size_t width, const char* pathPtr);

template<typename T>
inline bool store2dRawData(T* imageDataPtr, const size_t length, const size_t width, const char* pathPtr)
{
	FILE* file;

	//checks if the memory was allocated
	if (imageDataPtr == NULL)
		return false;

	//writing binary data to file
	if (fopen_s(&file, pathPtr, "wb") != 0)
		return false;

	const size_t pointsInSlice = length * width;

	fwrite(imageDataPtr, sizeof(T), pointsInSlice, file);

	fclose(file);

	return true;
}


template<typename T>
bool manageRAWFile2D(T* imageDataPtr, const size_t length, const size_t width, const char* pathPtr, Operation operation, bool revert);

template <typename T>
inline bool manageRAWFile2D(T* imageDataPtr, const size_t length, const size_t width, const char* pathPtr, Operation operation, bool revert) {

	bool status = false; // Initial Status, only changed to true if the operation is successful

	switch (operation)
	{
	case LOAD_DATA:
		status = load2dArrayRAW<T>(imageDataPtr, length, width, pathPtr, revert);
		break;
	case STORE_DATA:
		status = store2dRawData<T>(imageDataPtr, length, width, pathPtr);
		break;

	default:
		break;
	}

	return status;
}
