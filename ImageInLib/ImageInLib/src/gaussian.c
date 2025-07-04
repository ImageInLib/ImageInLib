#include <stdlib.h> //malloc, free
#include "math.h" // exp, pow, sqrt
#include "gaussian.h"

bool generateGaussianKernel(dataType* kernel, const dataType sigma) {

	if (kernel == NULL || sigma <= 0) {
		return false;
	}

	dataType sum = 0.0;
	size_t kernelSize = (size_t)(2 * sigma + 1);
	for (size_t i = 0; i < kernelSize; i++) 
	{
		for(size_t j = 0; j < kernelSize; j++) 
		{
			size_t xd = x_new(i, j, kernelSize);
			kernel[xd] = (dataType)(exp(-( pow(i, 2) + pow(j, 2) )) / (2.0 * pow(sigma, 2)));
			sum += kernel[xd];
		}
	}

	//The empirical normalization is considered to avoid accumulated errors
	// Analitical normalization ---> 1.0 / pow(sqrt(2.0 * M_PI * sigma * sigma), d);  d is the dimension
	for (size_t i = 0; i < kernelSize * kernelSize; i++) {
		kernel[i] /= sum;
	}
	return true;
}

bool gaussianSmoothing2D(dataType* imageDataPtr, dataType* smoothImageDataPtr, const size_t length, const size_t width, const dataType sigma) {

	if (imageDataPtr == NULL || smoothImageDataPtr == NULL || sigma <= 0) {
		return false;
	}

	size_t i, j, k, xd;
	const size_t kernelSize = (size_t)(2 * sigma + 1);
	const size_t p = kernelSize / 2;
	if (kernelSize < 3) {
		return false;
	}

	dataType* kernel = (dataType*)malloc(kernelSize * kernelSize * sizeof(dataType));
	if (kernel == NULL) {
		return false;
	}

	//generate the Gaussian Kernel
	generateGaussianKernel(kernel, sigma);

	const size_t pLength = length + 2 * p;
	const size_t pWidth = width + 2 * p;
	dataType* imageDataExt = (dataType*)malloc(pLength * pWidth * sizeof(dataType));
	for(i = 0; i < pLength * pWidth; i++) {
		imageDataExt[i] = 0.0;
	}

	//Copy to extended area
	size_t i_ext, j_ext, k_ext, xd_ext;
	for (i = 0, i_ext = p; i < length; i++, i_ext++) {
		for (j = 0, j_ext = p; j < width; j++, j_ext++) {
			xd = x_new(i, j, length);
			xd_ext = x_new(i_ext, j_ext, pLength);
			imageDataExt[xd_ext] = imageDataPtr[xd];
		}
	}

	reflection2DB(imageDataExt, length, width, p);

	//Convolution
	size_t im, jm, km, xdm, xdn;
	for (i = 0, i_ext = p; i < length; i++, i_ext++) {
		for (j = 0, j_ext = p; j < width; j++, j_ext++) {
			xd = x_new(i, j, length);
			xd_ext = x_new(i_ext, j_ext, pLength);
			for (size_t ki = 0; ki < kernelSize; ki++) {
				for (size_t kj = 0; kj < kernelSize; kj++) {
					im = i_ext + ki - p;
					jm = j_ext + kj - p;
					xdm = x_new(im, jm, pLength);
					xdn = x_new(ki, kj, kernelSize);
					smoothImageDataPtr[xd] += kernel[xdn] * imageDataExt[xdm];
				}
			}
		}
	}

	free(imageDataExt);
	free(kernel);

	return true;
}
