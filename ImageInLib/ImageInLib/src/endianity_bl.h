#ifdef __cplusplus
extern "C" {
#endif

	/*
	* Author: Jozef Urbán
	* Purpose: Business logic declarations for data endianity conversion
	* Language:  C
	*/

#pragma once
#include <uchar.h>

	short swap2B(short val);// Byte swap
	int swap4BToMidBigEndian(int val);// Byte couples swap

	// Reverts bytes pointed by dataPtr
	// dataPtr - pointer to data which data will be reverted
	// dataSize - size of data type
	void revertBytes(void * dataPtr, const size_t dataSize);
	/*
	* Author: Jozef Urbán
	* Purpose: Business logic declarations for data endianity conversion
	* Language:  C
	*/

#pragma once
#include <uchar.h>

	short swap2B(short val);// Byte swap
	int swap4BToMidBigEndian(int val);// Byte couples swap

									  // Reverts bytes pointed by dataPtr
									  // dataPtr - pointer to data which data will be reverted
									  // dataSize - size of data type
	void revertBytes(void * dataPtr, const size_t dataSize);

#ifdef __cplusplus
}
#endif