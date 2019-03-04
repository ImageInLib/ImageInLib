#include "transformation.h"
#include "transformations.h"

void transformImage(ImageData image, void **outputTransformedImage, affineTransformParameters affineTransformParam, transformationMethod transform)
{
	switch (transform)
	{
	case FORWARD_TRANSFORM:
		transform3DImage(image.imageDataPtr, (dataType **) outputTransformedImage, affineTransformParam.translation, affineTransformParam.scaling, affineTransformParam.rotation, image.height, image.length, image.width, affineTransformParam.imageBackground, affineTransformParam.centroid);
		break;
	case INVERSE_TRANSFORM:
		transformInverse3DImage(image.imageDataPtr, (dataType **) outputTransformedImage, affineTransformParam.translation, affineTransformParam.scaling, affineTransformParam.rotation, image.height, image.length, image.width, affineTransformParam.imageBackground, affineTransformParam.centroid);
		break;
	default:
		break;
	}
}