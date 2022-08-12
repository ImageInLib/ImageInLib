#pragma once

#include<iostream>

template <typename T>
bool store3dRawData(T* imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, const char* pathPtr);

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
