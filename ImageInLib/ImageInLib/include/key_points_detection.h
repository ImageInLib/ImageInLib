#pragma once
#include <iostream>
#include <vector>
#include "common_functions.h"
#include "front_propagation.h"

bool frontPropagationWithKeyPointDetection(Image_Data2D inputImageData, dataType* actionPtr, dataType* potentialPtr, Point2D* endPoints, const dataType LengthKeyPoints, std::vector<Point2D>& key_points);
