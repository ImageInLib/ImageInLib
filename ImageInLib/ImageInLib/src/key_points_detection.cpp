#include <iostream>
#include <vector>
#include "key_points_detection.h"
#include "imageInterpolation.h"

bool frontPropagationWithKeyPointDetection(Image_Data2D inputImageData, dataType* actionPtr, dataType* potentialPtr, Point2D* endPoints, const dataType LengthKeyPoints, std::vector<Point2D>& key_points) {

	if (inputImageData.imageDataPtr == NULL || actionPtr == NULL || potentialPtr == NULL || endPoints == NULL) {
		return false;
	}

	const size_t length = inputImageData.height;
	const size_t width = inputImageData.width;

	PixelSpacing spacing = inputImageData.spacing;
	if (endPoints[0].x < 0 || endPoints[0].x > length || endPoints[0].y < 0 || endPoints[0].y > width) {
		std::cout << "Error in the input seed point" << std::endl;
		return false;
	}
	if (endPoints[1].x < 0 || endPoints[1].x > length || endPoints[1].y < 0 || endPoints[1].y > width) {
		std::cout << "Error in the input seed point" << std::endl;
		return false;
	}

	std::vector <pointFastMarching2D> inProcess;
	size_t i = 0, j = 0, k = 0, dim2D = length * width;

	short* labelArray = new short[dim2D];

	//Initialization
	//All the points are notProcessed ---> label = 3
	for (k = 0; k < dim2D; k++) {
		actionPtr[k] = INFINITY;
		labelArray[k] = 3;
	}

	pointFastMarching2D current;
	i = (size_t)endPoints[0].x;
	j = (size_t)endPoints[0].y;
	size_t currentIndx = x_new(i, j, length);
	actionPtr[currentIndx] = 0.0;
	labelArray[currentIndx] = 1;

	//find the neighbours of the initial point add add them to inProcess
	size_t length_minus = length - 1, width_minus = width - 1;

	//West
	if (i > 0 && i < length && j >= 0 && j < width) {
		size_t iminus = i - 1;
		size_t indxWest = x_new(iminus, j, length);
		dataType x = selectX(actionPtr, length, width, iminus, j);
		dataType y = selectY(actionPtr, length, width, iminus, j);
		dataType coefSpeed = potentialPtr[indxWest];
		dataType dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D WestNeighbor = { iminus, j, dWest };
		actionPtr[indxWest] = dWest;
		inProcess.push_back(WestNeighbor);
		labelArray[indxWest] = 2;
	}

	//East
	if (i < length_minus && i >= 0 && j >= 0 && j < width) {
		size_t iplus = i + 1;
		size_t indxEast = x_new(iplus, j, length);
		dataType x = selectX(actionPtr, length, width, iplus, j);
		dataType y = selectY(actionPtr, length, width, iplus, j);
		dataType coefSpeed = potentialPtr[indxEast];
		dataType dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D EastNeighbor = { iplus, j, dEast };
		actionPtr[indxEast] = dEast;
		inProcess.push_back(EastNeighbor);
		labelArray[indxEast] = 2;
	}

	//North
	if (j > 0 && j < width && i >= 0 && i < length) {
		size_t jminus = j - 1;
		size_t indxNorth = x_new(i, jminus, length);
		dataType x = selectX(actionPtr, length, width, i, jminus);
		dataType y = selectY(actionPtr, length, width, i, jminus);
		dataType coefSpeed = potentialPtr[indxNorth];
		dataType dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D NorthNeighbor = { i, jminus, dNorth };
		actionPtr[indxNorth] = dNorth;
		inProcess.push_back(NorthNeighbor);
		labelArray[indxNorth] = 2;
	}

	//South
	if (j >= 0 && j < width_minus && i >= 0 && i < length) {
		size_t jplus = j + 1;
		size_t indxSouth = x_new(i, jplus, length);
		dataType x = selectX(actionPtr, length, width, i, jplus);
		dataType y = selectY(actionPtr, length, width, i, jplus);
		dataType coefSpeed = potentialPtr[indxSouth];
		dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D SouthNeighbor = { i, jplus, dSouth };
		actionPtr[indxSouth] = dSouth;
		inProcess.push_back(SouthNeighbor);
		labelArray[indxSouth] = 2;
	}

	heapifyVector2D(inProcess);
	size_t label = 0;

	int l = 0, m = 0;
	size_t nbSourcePoint = 0;
	dataType max_weighted_distance = 0.0;

	size_t iEnd = (size_t)endPoints[1].x;
	size_t jEnd = (size_t)endPoints[1].y;

	//Visualize the front propagation
	size_t id_keyPoint = 1;
	std::string storing_path;
	std::vector<Point2D> savingList;

	double distanceToCurrentSourcePoint = 0.0;
	//Set the starting point as initial source point
	Point2D currentSourcePoint = { (dataType)i, (dataType)j };
	key_points.push_back(currentSourcePoint);

	while (inProcess.size() > 0) {

		//processed the point with minimum distance
		current = inProcess[0];
		i = current.x;
		j = current.y;
		currentIndx = x_new(i, j, length);
		labelArray[currentIndx] = 1;

		//Exit of the while loop if the end point is reached
		if (labelArray[x_new(iEnd, jEnd, length)] == 1) {
			break;
		}

		//check if the current point is a key point
		Point2D pSource = { i, j };
		Point2D pSourceReal = getRealCoordFromImageCoord2D(pSource, inputImageData.origin, inputImageData.spacing, inputImageData.orientation);
		Point2D pCurrent = getRealCoordFromImageCoord2D(currentSourcePoint, inputImageData.origin, inputImageData.spacing, inputImageData.orientation);
		distanceToCurrentSourcePoint = getPoint2DDistance(pCurrent, pSourceReal);
		
		if (distanceToCurrentSourcePoint >= LengthKeyPoints) {

			//If the condition is true ---> new key point is found so we need to initilize it neighbors
			currentSourcePoint = pSource;
			key_points.push_back(currentSourcePoint);
			actionPtr[currentIndx] = 0;

			//Initialize all the points inside the narrow band as not processed ---> label = 3
			for (size_t it = 0; it < inProcess.size(); it++) {
				labelArray[x_new(inProcess[it].x, inProcess[it].y, length)] = 1;
			}

			//West
			if (i > 0 && i < length && j >= 0 && j < width)
			{
				size_t iminus = i - 1;
				size_t indxWest = x_new(iminus, j, length);
				actionPtr[indxWest] = INFINITY;
				labelArray[indxWest] = 3;
			}

			//East
			if (i < length_minus && i >= 0 && j >= 0 && j < width)
			{
				size_t iplus = i + 1;
				size_t indxEast = x_new(iplus, j, length);
				actionPtr[indxEast] = INFINITY;
				labelArray[indxEast] = 3;
			}

			//North
			if (j > 0 && j < width && i >= 0 && i < length)
			{
				size_t jminus = j - 1;
				size_t indxNorth = x_new(i, jminus, length);
				actionPtr[indxNorth] = INFINITY;
				labelArray[indxNorth] = 3;
			}

			//South
			if (j >= 0 && j < width_minus && i >= 0 && i < length)
			{
				size_t jplus = j + 1;
				size_t indxSouth = x_new(i, jplus, length);
				actionPtr[indxSouth] = INFINITY;
				labelArray[indxSouth] = 3;
			}

			inProcess.clear();
			id_keyPoint++;

		}
		else {
			//actionMapStr.imageDataPtr[currentIndx] = current.arrival;
			deleteRootHeap2D(inProcess);
		}

		//====================

		//processed neighbors of the minimum in the narrow band

		//West
		if (i > 0 && i < length && j >= 0 && j < width) {
			size_t iminus = i - 1;
			size_t indxWest = x_new(iminus, j, length);
			size_t label = labelArray[indxWest];
			if (label != 1) {
				dataType x = selectX(actionPtr, length, width, iminus, j);
				dataType y = selectY(actionPtr, length, width, iminus, j);
				dataType coefSpeed = potentialPtr[indxWest];
				dataType dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
				if (label == 3) {
					labelArray[indxWest] = 2;
					actionPtr[indxWest] = dWest;
					pointFastMarching2D WestNeighbor = { iminus, j, dWest };
					addPointHeap2D(inProcess, WestNeighbor);
				}
				else {
					if (dWest < actionPtr[indxWest]) {
						actionPtr[indxWest] = dWest;
						size_t pt_pos = getIndexFromHeap2D(inProcess, iminus, j);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dWest;
							heapifyUp2D(inProcess, pt_pos);
						}
					}
				}
			}
		}

		//East
		if (i < length_minus && i >= 0 && j >= 0 && j < width) {
			size_t iplus = i + 1;
			size_t indxEast = x_new(iplus, j, length);
			size_t label = labelArray[indxEast];
			if (label != 1) {
				dataType x = selectX(actionPtr, length, width, iplus, j);
				dataType y = selectY(actionPtr, length, width, iplus, j);
				dataType coefSpeed = potentialPtr[indxEast];
				dataType dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
				if (label == 3) {
					labelArray[indxEast] = 2;
					actionPtr[indxEast] = dEast;
					pointFastMarching2D EastNeighbor = { iplus, j, dEast };
					addPointHeap2D(inProcess, EastNeighbor);
				}
				else {
					if (dEast < actionPtr[indxEast]) {
						actionPtr[indxEast] = dEast;
						size_t pt_pos = getIndexFromHeap2D(inProcess, iplus, j);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dEast;
							heapifyUp2D(inProcess, pt_pos);
						}
					}
				}
			}
		}

		//North
		if (j > 0 && j < width && i >= 0 && i < length) {
			size_t jminus = j - 1;
			size_t indxNorth = x_new(i, jminus, length);
			size_t label = labelArray[indxNorth];
			if (label != 1) {
				dataType x = selectX(actionPtr, length, width, i, jminus);
				dataType y = selectY(actionPtr, length, width, i, jminus);
				dataType coefSpeed = potentialPtr[indxNorth];
				dataType dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
				if (label == 3) {
					labelArray[indxNorth] = 2;
					actionPtr[indxNorth] = dNorth;
					pointFastMarching2D NorthNeighbor = { i, jminus, dNorth };
					addPointHeap2D(inProcess, NorthNeighbor);
				}
				else {
					if (dNorth < actionPtr[indxNorth]) {
						actionPtr[indxNorth] = dNorth;
						size_t pt_pos = getIndexFromHeap2D(inProcess, i, jminus);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dNorth;
							heapifyUp2D(inProcess, pt_pos);
						}
					}
				}
			}
		}

		//South
		if (j >= 0 && j < width_minus && i >= 0 && i < length) {
			size_t jplus = j + 1;
			size_t indxSouth = x_new(i, jplus, length);
			size_t label = labelArray[indxSouth];
			if (label != 1) {
				dataType x = selectX(actionPtr, length, width, i, jplus);
				dataType y = selectY(actionPtr, length, width, i, jplus);
				dataType coefSpeed = potentialPtr[indxSouth];
				dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
				if (label == 3) {
					labelArray[indxSouth] = 2;
					actionPtr[indxSouth] = dSouth;
					pointFastMarching2D NorthNeighbor = { i, jplus, dSouth };
					addPointHeap2D(inProcess, NorthNeighbor);
				}
				else {
					if (dSouth < actionPtr[indxSouth]) {
						actionPtr[indxSouth] = dSouth;
						size_t pt_pos = getIndexFromHeap2D(inProcess, i, jplus);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dSouth;
							heapifyUp2D(inProcess, pt_pos);
						}
					}
				}
			}
		}

	}

	delete[] labelArray;

	return true;
}