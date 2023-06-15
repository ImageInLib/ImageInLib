/*
* Author: Jozef Urbán
* Purpose: Business logic definitions for data endianity conversion
* Language:  C
*/

#include <stdlib.h>
#include <memory.h>
#include "endianity_bl.h"

// Byte swap
short swap2B(short val)
{
	return (val << 8) | ((val >> 8) & 0xFF);
}

// Byte couples swap
int swap4BToMidBigEndian(int val)
{
	return (val << 16) | ((val >> 16) & 0xFFFF);
}

// Reverts bytes of pointed by dataPtr
void revertBytes(void * dataPtr, const size_t dataSize)
{
	char * output8B = malloc(dataSize);
	revertBytesEx(dataPtr, dataSize, output8B);
	free(output8B);
}

void revertBytesEx(void* dataPtr, const size_t dataSize, char* swapBuff) {
	
	char* inputToConvert = (char*)dataPtr;
	//swap the bytes into a temporary buffer

	for (size_t i = 0; i < dataSize; i++)
	{
		swapBuff[i] = inputToConvert[dataSize - i - 1];
	}

	memcpy(inputToConvert, swapBuff, dataSize);
}