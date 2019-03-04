#include "registration.h"
//==============================================================================
#include "shape_registration.h"
//==============================================================================

void imageRegistration(dataType **movingImage, dataType **fixedImage, dataType **registeredImage, const size_t imageHeight, const size_t imageLength, const size_t imageWidth, registrationParams params, optimizationMethod optimization)
{
	switch (optimization)
	{
	case GRADIENT_DESCENT:
		run_registration(movingImage, fixedImage, registeredImage, imageHeight, imageLength, imageWidth, params, 1);
		break;
	case STOCHASTIC_DESCENT:
		run_registration(movingImage, fixedImage, registeredImage, imageHeight, imageLength, imageWidth, params, 2);
		break;
	default:
		break;
	}
}