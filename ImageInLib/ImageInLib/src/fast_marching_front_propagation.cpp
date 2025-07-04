#include "fast_marching_front_propagation.h"
#include "partial_front_propagation.h"
#include "front_propagation.h"
#include "double_front_propagation.h"
#include "key_points_detection.h"

void fastMarchingfrontPropagation2D(Image_Data2D inputImageData, dataType* firstActionMapPtr, dataType* secondActionMapPtr, dataType* potentialPtr, Point2D* endPoints, const dataType LengthKeyPoints, std::vector<Point2D>& keyPoints, const PropagationType pType)
{
	switch (pType)
	{
    case PARTIAL_FRONT_PROPAGATION:
		partialFrontPropagation2D(inputImageData, firstActionMapPtr, potentialPtr, endPoints);
		break;
	case FRONT_PROPAGATION:
		frontPropagation2D(inputImageData, firstActionMapPtr, potentialPtr, endPoints);
		break;
	case DOUBLE_FRONT_PROPAGATION:
		doubleFrontPropagation2D(inputImageData, firstActionMapPtr, secondActionMapPtr, potentialPtr, endPoints);
		break;
	case KEY_POINT_DETECTION:
		frontPropagationWithKeyPointDetection(inputImageData, firstActionMapPtr, potentialPtr, endPoints, LengthKeyPoints, keyPoints);
		break;
	default:
		break;
	}
}

void fastMarchingfrontPropagation(Image_Data inputImageData, dataType** actionMapPtr, dataType** potentialPtr, Point3D* endPoints, const PropagationType pType)
{
	switch (pType)
	{
		case PARTIAL_FRONT_PROPAGATION:
			partialFrontPropagation(inputImageData, actionMapPtr, potentialPtr, endPoints);
			break;
			/*
		case FRONT_PROPAGATION:
			frontPropagation(inputImageData, potential, endPoints);
			break;
		case DOUBLE_FRONT_PROPAGATION:
			doubleFrontPropagation(inputImageData, potential, endPoints);
			break;
		case KEY_POINT_DETECTION:
			frontPropagationsWithKeyPointDetection(inputImageData, potential, endPoints);
			break;
		case DISTANCE_MAP:
			distanceMapFastMarching(inputImageData, potential, endPoints);
			break*/;
	default:
		break;
	}
}

void computePotential2D(dataType* imageDataPtr, dataType* potentialPtr, const size_t length, const size_t width, Point2D* endPoints, const dataType epsilon) {
	if(imageDataPtr == NULL || potentialPtr == NULL || endPoints == NULL) {
		return; // Handle null pointers
	}

	size_t x1 = (size_t)endPoints[0].x;
	size_t y1 = (size_t)endPoints[0].y;
	size_t x2 = (size_t)endPoints[1].x;
	size_t y2 = (size_t)endPoints[1].y;
	dataType seedValue = ( imageDataPtr[x_new(x1, y1, length)] + imageDataPtr[x_new(x2, y2, length)] ) / 2.0;

	size_t dim2D = length * width;
	for(size_t i = 0; i < dim2D; i++) {
		potentialPtr[i] = epsilon + fabs(seedValue - imageDataPtr[i]);
	}
}