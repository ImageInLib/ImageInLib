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
			for (i = 0; i < xDim; i++){
				for (j = 0; j < yDim; j++){
					// 2D to 1D representation for i, j
					xd = x_new(i, j, xDim);
					fread(&value, sizeof(T), 1, file);
					imageDataPtr[k][xd] = (T)value;
				}
			}
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

	fclose(file);
	return true;
}

//template<typename T>
//bool load3dDataArrayVTK(T** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, const char* pathPtr);

//template <typename T>
//inline bool load3dDataArrayVTK(T** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, const char* pathPtr)
//{
//	//checks if the memory was allocated
//	if (imageDataPtr == NULL)
//		return false;
//
//	if (pathPtr == NULL)
//		return false;
//
//	char line1[100], line2[100], line3[100], line4[100], line5[100], line6[100], line7[100], line8[100], line9[100], line10[100];
//	size_t i, j, k;
//	size_t x; //x = x_new(i, j, xDim);
//	T value;
//	FILE* file;
//
//	//Reading data from file
//	if (fopen_s(&file, pathPtr, "rb") != 0) {
//		return false;
//	}
//	else {
//		//Read header
//		fgets(line1, 100, file);
//		fgets(line2, 100, file);
//		fgets(line3, 100, file);
//		fgets(line4, 100, file);
//		fgets(line5, 100, file);
//		fgets(line6, 100, file);
//		fgets(line7, 100, file);
//		fgets(line8, 100, file);
//		fgets(line9, 100, file);
//		fgets(line10, 100, file);
//
//		//Read other data
//		for (k = 0; k < zDim; k++)
//		{
//			for (i = 0; i < xDim; i++)
//			{
//				for (j = 0; j < yDim; j++)
//				{
//					// 2D to 1D representation for i, j
//					x = x_new(i, j, xDim);
//
//					//fscanf(file, "%d", &value);
//					value = getc(file);
//					imageDataPtr[k][x] = value;
//				}
//			}
//		}
//	}
//	fclose(file);
//	return true;
//}
