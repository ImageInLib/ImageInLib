#pragma once
#include "common_functions.h"

void doubleFrontPropagation2D(Image_Data2D imageData, dataType* actionFirstFront, dataType* actionSecondFront, dataType* potentialPtr, Point2D* endPoints);
