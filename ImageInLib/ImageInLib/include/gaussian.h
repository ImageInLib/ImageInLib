#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"

	bool generateGaussianKernel(dataType* kernel, const dataType sigma);
	
	bool gaussianSmoothing2D(dataType* imageDataPtr, dataType* smoothImageDataPtr, const size_t length, const size_t width, const dataType sigma);

#ifdef __cplusplus
}
#endif
