#ifdef __cplusplus
extern "C" {
#endif

	/*
	* Authors: Polycarp Omondi Okock and Markjoe Olunna UBA
	* Purpose: ImageInLife project - 4D Image Segmentation Methods
	* Language:  C
	*/
#pragma once
#include <stdbool.h>
#include "common_functions.h"

	bool additive3dNoise_UC(unsigned char ** array3DPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType C);

	bool additive2dNoise_UC(unsigned char * array2DPtr, const size_t xDim, const size_t yDim, dataType C, bool flag);

	bool saltAndPepper3dNoise_UC(unsigned char ** array3DPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType K);

	bool saltAndPepper2dNoise_UC(unsigned char * array2DPtr, const size_t xDim, const size_t yDim, dataType K, bool flag);

	bool additive3dNoise_D(dataType** array3DPtr, const size_t xDim, const size_t yDim, const size_t zDim, int C, dataType fgMin, dataType bgMax);

	bool additive2dNoise_D(dataType * array2DPtr, const size_t xDim, const size_t yDim, dataType C, bool flag);

	bool saltAndPepper3dNoise_D(dataType ** array3DPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType density, const dataType pepper);

	bool saltAndPepper2dNoise_D(dataType* array2DPtr, const size_t xDim, const size_t yDim, dataType density, const dataType pepper);

	/*
	* Multiplicative noise adds noise to imageDataPtr
	* imageHeight represent the Z coordinate, imageWidth is the X*Y coordinate - 2D representation of 3D
	* Variance for a Continuous Uniform Distribution
	* b - upper value, a - lower boundary (assume to be Zero)
	* Variance is used to evaluate b, then be and a used to generate random number n
	* N is used to add multiplicative noise to the original image
	*/
	/*void addMultiplicativeNoise(dataType** imageDataPtr, const size_t imageHeight, const size_t imageWidth, dataType variance);*/
	/*
	* Structural noise generate adds noise to imageDataPtr
	* imageDataPtr contains points that we will add noise to
	* imageHeight represent the Z coordinate, imageWidth is the X*Y coordinate - 2D representation of 3D
	* Uses a sinusoid function to generate regular noise
	*/
	/*void addStructuralNoise(dataType** imageDataPtr, const size_t imageHeight, const size_t imageWidth);*/

#ifdef __cplusplus
}
#endif
