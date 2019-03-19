#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "transform_params.h"

	typedef enum
	{
		FORWARD_TRANSFORM = 1,
		INVERSE_TRANSFORM
	} transformationMethod;

	void transformImage(ImageData image, dataType **, affineTransformParameters affineTransformParam, transformationMethod transform);

#ifdef __cplusplus
}
#endif