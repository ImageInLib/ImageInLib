#include "transformation.h"
#include "transformations.h"

void transformImage(Image_Data image, dataType **outputTransformedImage, Affine_Transform_Parameters affineTransformParam, TransformationMethod transform)
{
	switch (transform)
	{
	case FORWARD_TRANSFORM:
		transform3DImage(image.imageDataPtr, outputTransformedImage, affineTransformParam.translation, affineTransformParam.scaling, affineTransformParam.rotation, image.height, image.length, image.width, affineTransformParam.imageBackground, affineTransformParam.centroid);
		break;
	case INVERSE_TRANSFORM:
		transformInverse3DImage(image.imageDataPtr, outputTransformedImage, affineTransformParam.translation, affineTransformParam.scaling, affineTransformParam.rotation, image.height, image.length, image.width, affineTransformParam.imageBackground, affineTransformParam.centroid);
		break;
	default:
		break;
	}
}