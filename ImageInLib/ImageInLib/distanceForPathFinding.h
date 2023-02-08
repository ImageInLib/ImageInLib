#ifdef __cplusplus
extern "C" {
#endif

#pragma once

#include<iostream>
#include "common_functions.h"
#include "../src/data_load.h"
#include "../src/endianity_bl.h"
#include <stdio.h>
#include <string.h>

	typedef struct {
		size_t x, y;
	} Point2D;

	bool fastMarching2d(dataType* inputImagePtr, dataType* outputImage, const size_t imageHeight, const size_t imageWidth, Point2D* seedPoints);



#ifdef __cplusplus
}
#endif