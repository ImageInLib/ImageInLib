#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "../src/transform_params.h"

	typedef enum
	{
		FORWARD_TRANSFORM = 1,
		INVERSE_TRANSFORM
	} TransformationMethod;

	void transformImage(Image_Data image, dataType **, Affine_Transform_Parameters affineTransformParam, TransformationMethod transform);

#ifdef __cplusplus
}
#endif