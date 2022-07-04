#include "noiseGeneration.h"
#include "common_functions.h"
#include "noise_generator.h"

// length == xDim, width == yDim, height == zDim
void addNoiseToImage(dataType ** array3DPtr, const size_t xDim, const size_t yDim, const size_t zDim,
	dataType density, const NoiseType method)
{
	switch (method)
	{
	case SALT_AND_PEPPER:
		saltAndPepper3dNoise_D(array3DPtr, xDim, yDim, zDim, density);
		break;
	case ADDITIVE_NOISE:
		additive3dNoise_D(array3DPtr, xDim, yDim, zDim, density);
		break;
	case MULTIPLICATIVE_NOISE:
		addMultiplicativeNoise(array3DPtr, zDim, yDim, density);
		break;
	default:
		break;
	}
}