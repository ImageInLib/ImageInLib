#pragma once
#include "transform_params.h"

typedef enum
{
	FORWARD_TRANSFORM = 1,
	INVERSE_TRANSFORM
} transformationMethod;

void transformImage(ImageData image, void **, affineTransformParameters affineTransformParam, transformationMethod transform);