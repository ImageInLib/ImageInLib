
#pragma once

#include "front_propagation.h"

bool partialFrontPropagation2D(Image_Data2D imageData, dataType* actionMapPtr, dataType* potentialPtr, Point2D* endPoints);

bool partialFrontPropagation(Image_Data inputImageData, dataType** actionMapPtr, dataType** potentialPtr, Point3D* endPoints);
