
#pragma once
#include "front_propagation.h"

	typedef enum
	{
		PARTIAL_FRONT_PROPAGATION = 1,
		FRONT_PROPAGATION = 2,
		DOUBLE_FRONT_PROPAGATION = 3,
		KEY_POINT_DETECTION = 4,
		DISTANCE_MAP = 5
	} PropagationType;

	void fastMarchingfrontPropagation2D(Image_Data2D inputImageData, dataType* firstActionMapPtr, dataType* secondActionMapPtr, dataType* potentialPtr, Point2D* endPoints, const dataType LengthKeyPoints, std::vector<Point2D>& keyPoints, const PropagationType pType);

	void fastMarchingfrontPropagation(Image_Data inputImageData, dataType** actionMapPtr, dataType** potentialPtr, Point3D* endPoints, const PropagationType pType);

	void computePotential2D(dataType* imageDataPtr, dataType* potentialPtr, const size_t length, const size_t width, Point2D* endPoints, const dataType epsilon);
