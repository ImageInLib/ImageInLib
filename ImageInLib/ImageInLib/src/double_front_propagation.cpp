#include <iostream>
#include <vector>
#include "front_propagation.h"
#include "double_front_propagation.h"

void doubleFrontPropagation2D(Image_Data2D imageData, dataType* actionFirstFront, dataType* actionSecondFront, dataType* potentialPtr, Point2D* endPoints) {

	if (imageData.imageDataPtr == NULL || actionFirstFront == NULL || actionSecondFront == NULL || potentialPtr == NULL || endPoints == NULL) {
		return;// false;
	}
	const size_t length = imageData.height;
	const size_t width = imageData.width;
	PixelSpacing spacing = imageData.spacing;

	size_t dim2D = length * width;
	size_t length_minus = length - 1;
	size_t width_minus = width - 1;

	short* firstLabelArray = new short[dim2D];
	short* secondLabelArray = new short[dim2D];
	if (firstLabelArray == NULL || secondLabelArray == NULL || actionFirstFront == NULL || actionSecondFront == NULL) {
		return;// false;
	}

	//Initialize action
	for (size_t k = 0; k < dim2D; k++) {
		actionFirstFront[k] = INFINITY;
		actionSecondFront[k] = INFINITY;
		firstLabelArray[k] = 3; //3 ---> not processed
		secondLabelArray[k] = 3; //3 ---> not processed
	}

	size_t x1 = (size_t)endPoints[0].x;
	size_t y1 = (size_t)endPoints[0].y;
	if (x1 >= length || y1 >= width) {
		delete[] firstLabelArray;
		delete[] secondLabelArray;
		return;// false; //Invalid end point
	}
	firstLabelArray[x_new(x1, y1, length)] = 1; //1 ---> already processed
	actionFirstFront[x_new(x1, y1, length)] = 0.0;

	size_t x2 = (size_t)endPoints[1].x;
	size_t y2 = (size_t)endPoints[1].y;
	if (x2 >= length || y2 >= width) {
		delete[] secondLabelArray;
		delete[] firstLabelArray;
		return;// false; //Invalid end point
	}
	secondLabelArray[x_new(x2, y2, length)] = 1; //1 ---> already processed
	actionSecondFront[x_new(x2, y2, length)] = 0.0;

	std::vector<pointFastMarching2D> narrowBandFirstFront, narrowBandSecondFront;

	//Initialize neighbors for the first front

	if (x1 > 0 && y1 >= 0 && y1 < width) {
		size_t xminus = x1 - 1;
		dataType x = selectX(actionFirstFront, length, width, xminus, y1);
		dataType y = selectY(actionFirstFront, length, width, xminus, y1);
		size_t indxWest = x_new(xminus, y1, length);
		dataType coefSpeed = potentialPtr[indxWest];
		dataType dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D WestNeighbor = { xminus, y1, dWest };
		actionFirstFront[indxWest] = dWest;
		narrowBandFirstFront.push_back(WestNeighbor);
		firstLabelArray[indxWest] = 2; //2 ---> in process
	}

	if (x1 < length_minus && y1 >= 0 && y1 < width) {
		size_t xplus = x1 + 1;
		dataType x = selectX(actionFirstFront, length, width, xplus, y1);
		dataType y = selectY(actionFirstFront, length, width, xplus, y1);
		size_t indxEast = x_new(xplus, y1, length);
		dataType coefSpeed = potentialPtr[indxEast];
		dataType dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D EastNeighbor = { xplus, y1, dEast };
		actionFirstFront[indxEast] = dEast;
		narrowBandFirstFront.push_back(EastNeighbor);
		firstLabelArray[indxEast] = 2; //2 ---> in process
	}

	if (y1 > 0 && x1 >= 0 && x1 < length) {
		size_t yminus = y1 - 1;
		dataType x = selectX(actionFirstFront, length, width, x1, yminus);
		dataType y = selectY(actionFirstFront, length, width, x1, yminus);
		size_t indxNorth = x_new(x1, yminus, length);
		dataType coefSpeed = potentialPtr[indxNorth];
		dataType dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D WestNeighbor = { x1, yminus, dNorth };
		actionFirstFront[indxNorth] = dNorth;
		narrowBandFirstFront.push_back(WestNeighbor);
		firstLabelArray[indxNorth] = 2; //2 ---> in process
	}

	if (y1 < width_minus && x1 >= 0 && x1 < length) {
		size_t yplus = y1 + 1;
		dataType x = selectX(actionFirstFront, length, width, x1, yplus);
		dataType y = selectY(actionFirstFront, length, width, x1, yplus);
		size_t indxSouth = x_new(x1, yplus, length);
		dataType coefSpeed = potentialPtr[indxSouth];
		dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D EastNeighbor = { x1, yplus, dSouth };
		actionFirstFront[indxSouth] = dSouth;
		narrowBandFirstFront.push_back(EastNeighbor);
		firstLabelArray[indxSouth] = 2; //2 ---> in process
	}

	//Initialize neighbors for the second front

	if (x2 > 0 && y2 >= 0 && y2 < width) {
		size_t xminus = x2 - 1;
		dataType x = selectX(actionSecondFront, length, width, xminus, y2);
		dataType y = selectY(actionSecondFront, length, width, xminus, y2);
		size_t indxWest = x_new(xminus, y2, length);
		dataType coefSpeed = potentialPtr[indxWest];
		dataType dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D WestNeighbor = { xminus, y2, dWest };
		actionSecondFront[indxWest] = dWest;
		narrowBandSecondFront.push_back(WestNeighbor);
		secondLabelArray[indxWest] = 2; //2 ---> in process
	}

	if (x2 < length_minus && y2 >= 0 && y2 < width) {
		size_t xplus = x2 + 1;
		dataType x = selectX(actionSecondFront, length, width, xplus, y2);
		dataType y = selectY(actionSecondFront, length, width, xplus, y2);
		size_t indxEast = x_new(xplus, y2, length);
		dataType coefSpeed = potentialPtr[indxEast];
		dataType dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D EastNeighbor = { xplus, y2, dEast };
		actionSecondFront[indxEast] = dEast;
		narrowBandSecondFront.push_back(EastNeighbor);
		secondLabelArray[indxEast] = 2; //2 ---> in process
	}

	if (y2 > 0 && x2 >= 0 && x2 < length) {
		size_t yminus = y2 - 1;
		dataType x = selectX(actionSecondFront, length, width, x2, yminus);
		dataType y = selectY(actionSecondFront, length, width, x2, yminus);
		size_t indxNorth = x_new(x2, yminus, length);
		dataType coefSpeed = potentialPtr[indxNorth];
		dataType dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D NorthNeighbor = { x2, yminus, dNorth };
		actionSecondFront[indxNorth] = dNorth;
		narrowBandSecondFront.push_back(NorthNeighbor);
		secondLabelArray[indxNorth] = 2; //2 ---> in process
	}

	if (y2 < width_minus && x2 >= 0 && x2 < length) {
		size_t yplus = y2 + 1;
		dataType x = selectX(actionSecondFront, length, width, x2, yplus);
		dataType y = selectY(actionSecondFront, length, width, x2, yplus);
		size_t indxSouth = x_new(x2, yplus, length);
		dataType coefSpeed = potentialPtr[indxSouth];
		dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D SouthNeighbor = { x2, yplus, dSouth };
		actionSecondFront[indxSouth] = dSouth;
		narrowBandSecondFront.push_back(SouthNeighbor);
		secondLabelArray[indxSouth] = 2; //2 ---> in process
	}

	//heapify 2D vector
	heapifyVector2D(narrowBandFirstFront);
	heapifyVector2D(narrowBandSecondFront);

	dataType max_save_action = 0.0;
	size_t nb_computed_points = 0;

	while (narrowBandFirstFront.size() > 0 && narrowBandSecondFront.size() > 0) {

		pointFastMarching2D firstFrontPoint = narrowBandFirstFront[0];
		x1 = firstFrontPoint.x;
		y1 = firstFrontPoint.y;
		size_t indexFirst = x_new(x1, y1, length);
		//Exit the while loop when the fronts have met
		if (secondLabelArray[indexFirst] == 1) {
			//The fronts have met
			endPoints[2].x = x1;
			endPoints[2].y = y1;
			break;
		}
		else {
			firstLabelArray[indexFirst] = 1;
		}

		if (firstFrontPoint.arrival > max_save_action) {
			max_save_action = firstFrontPoint.arrival;
		}
		deleteRootHeap2D(narrowBandFirstFront);
		nb_computed_points++;

		pointFastMarching2D secondFrontPoint = narrowBandSecondFront[0];
		x2 = secondFrontPoint.x;
		y2 = secondFrontPoint.y;
		size_t indexSecond = x_new(x2, y2, length);
		if (firstLabelArray[indexSecond] == 1) {
			//The fronts have met
			endPoints[2].x = x2;
			endPoints[2].y = y2;
			break;
		}
		else {
			secondLabelArray[indexSecond] = 1;
		}

		if (secondFrontPoint.arrival > max_save_action) {
			max_save_action = secondFrontPoint.arrival;
		}
		deleteRootHeap2D(narrowBandSecondFront);
		nb_computed_points++;

		//West first front
		if (x1 > 0 && x1 < length && y1 >= 0 && y1 < width) {
			size_t xminus = x1 - 1;
			size_t indxWest = x_new(xminus, y1, length);
			short label = firstLabelArray[indxWest];
			if (label != 1) {
				dataType x = selectX(actionFirstFront, length, width, xminus, y1);
				dataType y = selectY(actionFirstFront, length, width, xminus, y1);
				dataType coefSpeed = potentialPtr[indxWest];
				dataType dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D WestNeighbor = { xminus, y1, dWest };
				if (label == 3) {
					actionFirstFront[indxWest] = dWest;
					firstLabelArray[indxWest] = 2;
					addPointHeap2D(narrowBandFirstFront, WestNeighbor);
				}
				else {
					if (dWest < actionFirstFront[indxWest]) {
						actionFirstFront[indxWest] = dWest;
						size_t pIndex = getIndexFromHeap2D(narrowBandFirstFront, xminus, y1);
						if (pIndex != -1) {
							heapifyUp2D(narrowBandFirstFront, pIndex);
						}
					}
				}
			}
		}

		//West second front
		if (x2 > 0 && x2 < length && y2 >= 0 && y2 < width) {
			size_t xminus = x2 - 1;
			size_t indxWest = x_new(xminus, y2, length);
			short label = secondLabelArray[indxWest];
			if (label != 1) {
				dataType x = selectX(actionSecondFront, length, width, xminus, y2);
				dataType y = selectY(actionSecondFront, length, width, xminus, y2);
				dataType coefSpeed = potentialPtr[indxWest];
				dataType dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D WestNeighbor = { xminus, y2, dWest };
				if (label == 3) {
					actionSecondFront[indxWest] = dWest;
					secondLabelArray[indxWest] = 2;
					addPointHeap2D(narrowBandSecondFront, WestNeighbor);
				}
				else {
					if (dWest < actionSecondFront[indxWest]) {
						actionSecondFront[indxWest] = dWest;
						size_t pIndex = getIndexFromHeap2D(narrowBandSecondFront, xminus, y2);
						if (pIndex != -1) {
							heapifyUp2D(narrowBandSecondFront, pIndex);
						}
					}
				}
			}
		}

		//East first front
		if (x1 < length_minus && x1 >= 0 && y1 >= 0 && y1 < width) {
			size_t xplus = x1 + 1;
			size_t indxEast = x_new(xplus, y1, length);
			short label = firstLabelArray[indxEast];
			if (label != 1) {
				dataType x = selectX(actionFirstFront, length, width, xplus, y1);
				dataType y = selectY(actionFirstFront, length, width, xplus, y1);
				dataType coefSpeed = potentialPtr[indxEast];
				dataType dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D EastNeighbor = { xplus, y1, dEast };
				if (label == 3) {
					actionFirstFront[indxEast] = dEast;
					firstLabelArray[indxEast] = 2;
					addPointHeap2D(narrowBandFirstFront, EastNeighbor);
				}
				else {
					if (dEast < actionFirstFront[indxEast]) {
						actionFirstFront[indxEast] = dEast;
						size_t pIndex = getIndexFromHeap2D(narrowBandFirstFront, xplus, y1);
						if (pIndex != -1) {
							heapifyUp2D(narrowBandFirstFront, pIndex);
						}
					}
				}
			}
		}

		//East second front
		if (x2 < length_minus && x2 >= 0 && y2 >= 0 && y2 < width) {
			size_t xplus = x2 + 1;
			size_t indxEast = x_new(xplus, y2, length);
			short label = secondLabelArray[indxEast];
			if (label != 1) {
				dataType x = selectX(actionSecondFront, length, width, xplus, y2);
				dataType y = selectY(actionSecondFront, length, width, xplus, y2);
				dataType coefSpeed = potentialPtr[indxEast];
				dataType dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D EastNeighbor = { xplus, y2, dEast };
				if (label == 3) {
					actionSecondFront[indxEast] = dEast;
					secondLabelArray[indxEast] = 2;
					addPointHeap2D(narrowBandSecondFront, EastNeighbor);
				}
				else {
					if (dEast < actionSecondFront[indxEast]) {
						actionSecondFront[indxEast] = dEast;
						size_t pIndex = getIndexFromHeap2D(narrowBandSecondFront, xplus, y2);
						if (pIndex != -1) {
							heapifyUp2D(narrowBandSecondFront, pIndex);
						}
					}
				}
			}
		}

		//North first front
		if (y1 > 0 && y1 < width && x1 >= 0 && x1 < length) {
			size_t yminus = y1 - 1;
			size_t indxNorth = x_new(x1, yminus, length);
			short label = firstLabelArray[indxNorth];
			if (label != 1) {
				dataType x = selectX(actionFirstFront, length, width, x1, yminus);
				dataType y = selectY(actionFirstFront, length, width, x1, yminus);
				dataType coefSpeed = potentialPtr[indxNorth];
				dataType dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D NorthNeighbor = { x1, yminus, dNorth };
				if (label == 3) {
					actionFirstFront[indxNorth] = dNorth;
					firstLabelArray[indxNorth] = 2;
					addPointHeap2D(narrowBandFirstFront, NorthNeighbor);
				}
				else {
					if (dNorth < actionFirstFront[indxNorth]) {
						actionFirstFront[indxNorth] = dNorth;
						size_t pIndex = getIndexFromHeap2D(narrowBandFirstFront, x1, yminus);
						if (pIndex != -1) {
							heapifyUp2D(narrowBandFirstFront, pIndex);
						}
					}
				}
			}
		}

		//North second front
		if (y2 > 0 && y2 < width && x2 >= 0 && x2 < length) {
			size_t yminus = y2 - 1;
			size_t indxNorth = x_new(x2, yminus, length);
			short label = secondLabelArray[indxNorth];
			if (label != 1) {
				dataType x = selectX(actionSecondFront, length, width, x2, yminus);
				dataType y = selectY(actionSecondFront, length, width, x2, yminus);
				dataType coefSpeed = potentialPtr[indxNorth];
				dataType dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D NorthNeighbor = { x2, yminus, dNorth };
				if (label == 3) {
					actionSecondFront[indxNorth] = dNorth;
					secondLabelArray[indxNorth] = 2;
					addPointHeap2D(narrowBandSecondFront, NorthNeighbor);
				}
				else {
					if (dNorth < actionSecondFront[indxNorth]) {
						actionSecondFront[indxNorth] = dNorth;
						size_t pIndex = getIndexFromHeap2D(narrowBandSecondFront, x2, yminus);
						if (pIndex != -1) {
							heapifyUp2D(narrowBandSecondFront, pIndex);
						}
					}
				}
			}
		}

		//South first front
		if (y1 >= 0 && y1 < width_minus && x1 >= 0 && x1 < length) {
			size_t yplus = y1 + 1;
			size_t indxSouth = x_new(x1, yplus, length);
			short label = firstLabelArray[indxSouth];
			if (label != 1) {
				dataType x = selectX(actionFirstFront, length, width, x1, yplus);
				dataType y = selectY(actionFirstFront, length, width, x1, yplus);
				dataType coefSpeed = potentialPtr[indxSouth];
				dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D SouthNeighbor = { x1, yplus, dSouth };
				if (label == 3) {
					actionFirstFront[indxSouth] = dSouth;
					firstLabelArray[indxSouth] = 2;
					addPointHeap2D(narrowBandFirstFront, SouthNeighbor);
				}
				else {
					if (dSouth < actionFirstFront[indxSouth]) {
						actionFirstFront[indxSouth] = dSouth;
						size_t pIndex = getIndexFromHeap2D(narrowBandFirstFront, x1, yplus);
						if (pIndex != -1) {
							heapifyUp2D(narrowBandFirstFront, pIndex);
						}
					}
				}
			}
		}

		//South second front
		if (y2 >= 0 && y2 < width_minus && x2 >= 0 && x2 < length) {
			size_t yplus = y2 + 1;
			size_t indxSouth = x_new(x2, yplus, length);
			short label = secondLabelArray[indxSouth];
			if (label != 1) {
				dataType x = selectX(actionSecondFront, length, width, x2, yplus);
				dataType y = selectY(actionSecondFront, length, width, x2, yplus);
				dataType coefSpeed = potentialPtr[indxSouth];
				dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D SouthNeighbor = { x2, yplus, dSouth };
				if (label == 3) {
					actionSecondFront[indxSouth] = dSouth;
					secondLabelArray[indxSouth] = 2;
					addPointHeap2D(narrowBandSecondFront, SouthNeighbor);
				}
				else {
					if (dSouth < actionSecondFront[indxSouth]) {
						actionSecondFront[indxSouth] = dSouth;
						size_t pIndex = getIndexFromHeap2D(narrowBandSecondFront, x2, yplus);
						if (pIndex != -1) {
							heapifyUp2D(narrowBandSecondFront, pIndex);
						}
					}
				}
			}
		}

	}

	for (size_t k = 0; k < dim2D; k++) {
		if (actionFirstFront[k] == INFINITY) {
			actionFirstFront[k] = max_save_action + 1;
		}
		if (actionSecondFront[k] == INFINITY) {
			actionSecondFront[k] = max_save_action + 1;
		}
	}
	narrowBandFirstFront.clear();
	narrowBandSecondFront.clear();

	delete[] firstLabelArray;
	delete[] secondLabelArray;

	//return true;
}