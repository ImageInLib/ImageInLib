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
	char * inputToConvert = (char *)dataPtr;
	char * output8B = malloc(dataSize);

	//swap the bytes into a temporary buffer

	for (size_t i = 0; i < dataSize; i++)
	{
		output8B[i] = inputToConvert[dataSize - i - 1];
	}

	memcpy(inputToConvert, output8B, dataSize);

	free(output8B);
}