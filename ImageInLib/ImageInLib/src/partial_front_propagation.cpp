#include<iostream>
#include<vector>

#include "partial_front_propagation.h"

bool partialFrontPropagation2D(Image_Data2D imageData, dataType* actionMapPtr, dataType* potentialPtr, Point2D* endPoints) {

	const size_t length = imageData.height;
	const size_t width = imageData.width;
	PixelSpacing spacing = imageData.spacing;

	size_t dim2D = length * width;
	size_t length_minus = length - 1;
	size_t width_minus = width - 1;

	short* labelArray = new short[dim2D];

	if (imageData.imageDataPtr == NULL || actionMapPtr == NULL || potentialPtr == NULL || endPoints == NULL || labelArray == NULL) {
		return false;
	}

	std::vector<pointFastMarching2D> narrowBand;

	size_t i = (size_t)endPoints[0].x;
	size_t j = (size_t)endPoints[0].y;
	size_t currentIndx = x_new(i, j, length);

	//STEP 1
	//In labelAray we have : 1 ---> already processed, 2 ---> narrow band and 3 ---> not processed
	for (size_t k = 0; k < dim2D; k++) {
		actionMapPtr[k] = INFINITY;
		labelArray[k] = 3;
	}

	actionMapPtr[currentIndx] = 0.0;
	labelArray[currentIndx] = 1;

	//East
	if (j < width - 1 && i >= 0 && i < length) {
		size_t jplus = j + 1;
		dataType x = selectX(actionMapPtr, length, width, i, jplus);
		dataType y = selectY(actionMapPtr, length, width, i, jplus);
		size_t indxEast = x_new(i, jplus, length);
		dataType coefSpeed = potentialPtr[indxEast];
		dataType dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D EastNeighbor = { i, jplus, dEast };
		actionMapPtr[indxEast] = dEast;
		narrowBand.push_back(EastNeighbor);
		labelArray[indxEast] = 2;
	}

	//West
	if (j > 0 && i >= 0 && i < length) {
		size_t jminus = j - 1;
		dataType x = selectX(actionMapPtr, length, width, i, jminus);
		dataType y = selectY(actionMapPtr, length, width, i, jminus);
		size_t indxWest = x_new(i, jminus, length);
		dataType coefSpeed = potentialPtr[indxWest];
		dataType dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D WestNeighbor = { i, jminus, dWest };
		actionMapPtr[indxWest] = dWest;
		narrowBand.push_back(WestNeighbor);
		labelArray[indxWest] = 2;
	}

	//North
	if (j >= 0 && j < width && i > 0) {
		size_t iminus = i - 1;
		dataType x = selectX(actionMapPtr, length, width, iminus, j);
		dataType y = selectY(actionMapPtr, length, width, iminus, j);
		size_t indxNorth = x_new(iminus, j, length);
		dataType coefSpeed = potentialPtr[indxNorth];
		dataType dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D NorthNeighbor = { iminus, j, dNorth };
		actionMapPtr[indxNorth] = dNorth;
		narrowBand.push_back(NorthNeighbor);
		labelArray[indxNorth] = 2;
	}

	//South
	if (j >= 0 && j < width && i < length - 1) {
		size_t iplus = i + 1;
		dataType x = selectX(actionMapPtr, length, width, iplus, j);
		dataType y = selectY(actionMapPtr, length, width, iplus, j);
		size_t indxSouth = x_new(iplus, j, length);
		dataType coefSpeed = potentialPtr[indxSouth];
		dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D SouthNeighbor = { iplus, j, dSouth };
		actionMapPtr[indxSouth] = dSouth;
		narrowBand.push_back(SouthNeighbor);
		labelArray[indxSouth] = 2;
	}

	//heapify 2D vector
	heapifyVector2D(narrowBand);

	//Save points for visualization
	size_t x_final_point = (size_t)endPoints[1].x;
	size_t y_final_point = (size_t)endPoints[1].y;

	dataType max_save_action = 0.0;
	
	while (narrowBand.size() > 0) {

		pointFastMarching2D current = narrowBand[0];
		size_t i = current.x;
		size_t j = current.y;
		size_t currentIndx = x_new(i, j, length);
		labelArray[currentIndx] = 1;

		// Exit when the final point is reached by the front
		if (i == x_final_point && j == y_final_point) {
			labelArray[currentIndx] = 1; // Mark the seed point as processed
			break; 
		}

		if (actionMapPtr[currentIndx] > max_save_action) {
			max_save_action = actionMapPtr[currentIndx];
		}

		deleteRootHeap2D(narrowBand);

		//West
		if (i > 0 && i < length && j >= 0 && j < width) {
			size_t iminus = i - 1;
			size_t indxWest = x_new(iminus, j, length);
			short label = labelArray[indxWest];
			if (label != 1) {
				dataType x = selectX(actionMapPtr, length, width, iminus, j);
				dataType y = selectY(actionMapPtr, length, width, iminus, j);
				dataType coefSpeed = potentialPtr[indxWest];
				dataType dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D WestNeighbor = { iminus, j, dWest };
				if (label == 3) {
					actionMapPtr[indxWest] = dWest;
					labelArray[indxWest] = 2;
					addPointHeap2D(narrowBand, WestNeighbor);
				}
				else {
					if (dWest < actionMapPtr[indxWest]) {
						actionMapPtr[indxWest] = dWest;
						size_t pIndex = getIndexFromHeap2D(narrowBand, iminus, j);
						if (pIndex != -1) {
							heapifyUp2D(narrowBand, pIndex);
						}
					}
				}
			}
		}

		//East
		if (i >= 0 && i < length_minus && j >= 0 && j < width) {
			size_t iplus = i + 1;
			size_t indxEast = x_new(iplus, j, length);
			short label = labelArray[indxEast];
			if (label != 1) {
				dataType x = selectX(actionMapPtr, length, width, iplus, j);
				dataType y = selectY(actionMapPtr, length, width, iplus, j);
				dataType coefSpeed = potentialPtr[indxEast];
				dataType dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D EastNeighbor = { iplus, j, dEast };
				if (label == 3) {
					actionMapPtr[indxEast] = dEast;
					labelArray[indxEast] = 2;
					addPointHeap2D(narrowBand, EastNeighbor);
				}
				else {
					if (dEast < actionMapPtr[indxEast]) {
						actionMapPtr[indxEast] = dEast;
						size_t pIndex = getIndexFromHeap2D(narrowBand, iplus, j);
						if (pIndex != -1) {
							heapifyUp2D(narrowBand, pIndex);
						}
					}
				}
			}
		}

		//North
		if (j > 0 && j < width && i >= 0 && i < length) {
			size_t jminus = j - 1;
			size_t indxNorth = x_new(i, jminus, length);
			short label = labelArray[indxNorth];
			if (label != 1) {
				dataType x = selectX(actionMapPtr, length, width, i, jminus);
				dataType y = selectY(actionMapPtr, length, width, i, jminus);
				dataType coefSpeed = potentialPtr[indxNorth];
				dataType dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D NorthNeighbor = { i, jminus, dNorth };
				if (label == 3) {
					actionMapPtr[indxNorth] = dNorth;
					labelArray[indxNorth] = 2;
					addPointHeap2D(narrowBand, NorthNeighbor);
				}
				else {
					if (dNorth < actionMapPtr[indxNorth]) {
						actionMapPtr[indxNorth] = dNorth;
						size_t pIndex = getIndexFromHeap2D(narrowBand, i, jminus);
						if (pIndex != -1) {
							heapifyUp2D(narrowBand, pIndex);
						}
					}
				}
			}
		}

		//South
		if (j >= 0 && j < width_minus && i >= 0 && i < length) {
			size_t jplus = j + 1;
			size_t indxSouth = x_new(i, jplus, length);
			short label = labelArray[indxSouth];
			if (label != 1) {
				dataType x = selectX(actionMapPtr, length, width, i, jplus);
				dataType y = selectY(actionMapPtr, length, width, i, jplus);
				dataType coefSpeed = potentialPtr[indxSouth];
				dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D SouthNeighbor = { i, jplus, dSouth };
				if (label == 3) {
					actionMapPtr[indxSouth] = dSouth;
					labelArray[indxSouth] = 2;
					addPointHeap2D(narrowBand, SouthNeighbor);
				}
				else {
					if (dSouth < actionMapPtr[indxSouth]) {
						actionMapPtr[indxSouth] = dSouth;
						size_t pIndex = getIndexFromHeap2D(narrowBand, i, jplus);
						if (pIndex != -1) {
							heapifyUp2D(narrowBand, pIndex);
						}
					}
				}
			}
		}
	}

	//Set the action of the non-processed points to maximum action value + 1
	for (size_t k = 0; k < dim2D; k++) {
		if (labelArray[k] == 3) {
			actionMapPtr[k] = max_save_action + 1;
		}
	}
	narrowBand.clear();

	delete[] labelArray;
}

bool partialFrontPropagation(Image_Data inputImageData, dataType** actionMapPtr, dataType** potentialPtr, Point3D* endPoints) {
	//TO DO
	return true;
}