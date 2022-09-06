#include "noiseGeneration.h"
#include "common_functions.h"
#include "noise_generator.h"

// length == xDim, width == yDim, height == zDim
void addNoiseToImage(dataType** array3DPtr, const size_t xDim, const size_t yDim, const size_t zDim,
	NoiseParameters parameters, const NoiseType method)
{
	switch (method)
	{
	case SALT_AND_PEPPER:
		saltAndPepper3dNoise_D(array3DPtr, xDim, yDim, zDim, parameters.salt_pepper_density, parameters.bgMax);
		break;
	case ADDITIVE_NOISE:
		additive3dNoise_D(array3DPtr, xDim, yDim, zDim, parameters.additive_value, parameters.fgMin, parameters.bgMax);
		break;
	case MULTIPLICATIVE_NOISE:
		addMultiplicativeNoise(array3DPtr, zDim, xDim, yDim, parameters.multiplicative_variance, parameters.fgMin, parameters.bgMax);
		break;
	case STRUCTURAL_NOISE:
		addStructuralNoise(array3DPtr, zDim, xDim, yDim, parameters.fgMin, parameters.bgMax);
		break;
	default:
		break;
	}
}