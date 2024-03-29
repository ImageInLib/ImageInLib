#include "common_functions.h"
#include "distanceMaps.h"
#include "fast_marching.h"
#include "distance_function.h"

void computeDistanceMap(dataType **outputDistanceData, dataType **inputImageData, const size_t imageLength, const size_t imageWidth, const size_t imageHeight, Distance_Map_Params distParams, DistanceMapMethod method)
{
	switch (method)
	{
	case FAST_SWEEP:
		fastSweepingFunction_3D(outputDistanceData, inputImageData, imageLength, imageWidth, imageHeight, distParams.h, distParams.initValue, distParams.objectPixel);
		break;
	case FAST_MARCH:
		fastMarching(outputDistanceData, inputImageData, imageHeight, imageLength, imageWidth, distParams.objectPixel);
		break;
	case ROUY_TOURIN:
		rouyTourinFunction_3D(outputDistanceData, inputImageData, distParams.tolerance, imageLength, imageWidth, imageHeight, distParams.tau, distParams.h);
		break;
	case BRUTE_FORCE:
		bruteForceFunction_3D(outputDistanceData, inputImageData, imageLength, imageWidth, imageHeight, distParams.initValue, distParams.objectPixel);
		break;
	default:
		break;
	}
}