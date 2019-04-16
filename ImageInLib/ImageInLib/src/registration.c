#include "registration.h"
//==============================================================================
#include "shape_registration.h"
//==============================================================================

void imageRegistration(dataType **movingImage, dataType **fixedImage, dataType **registeredImage, const size_t imageHeight, const size_t imageLength, const size_t imageWidth, Registration_Params params, Optimization_Method optimization)
{
	run_registration(fixedImage, movingImage, registeredImage, imageHeight, imageLength, imageWidth, params, optimization);
}