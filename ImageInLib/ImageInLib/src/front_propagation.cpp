#include<iostream>
#include<vector>
#include<cmath>     // sqrt()
#include<algorithm> // min() , max()

#include "front_propagation.h"

void swap2dPoints(pointFastMarching2D* a, pointFastMarching2D* b) {
	pointFastMarching2D temp = *a;
	*a = *b;
	*b = temp;
}

void heapifyDown2D(std::vector<pointFastMarching2D>& in_Process, int pos) {
	int length_array = in_Process.size();
	int current = pos;
	int left_child = 2 * pos + 1;
	int right_child = 2 * pos + 2;
	dataType val_current = 0.0, val_left = 0.0, val_right = 0.0;
	if (current >= 0 && current < length_array) {
		val_current = in_Process[current].arrival;
	}
	if (left_child < length_array) {
		val_left = in_Process[left_child].arrival;
		if (val_left < val_current) {
			current = left_child;
			val_current = in_Process[current].arrival;
		}
	}
	if (right_child < length_array) {
		val_right = in_Process[right_child].arrival;
		if (val_right < val_current) {
			current = right_child;
		}
	}
	if (current != pos) {
		swap2dPoints(&in_Process[pos], &in_Process[current]);
		heapifyDown2D(in_Process, current);
	}
}

void heapifyUp2D(std::vector<pointFastMarching2D>& in_Process, int i) {
	int current = i;
	if (i > 0) {
		int parent = (i - 1) / 2;
		dataType val_current = in_Process[current].arrival;
		dataType val_parent = in_Process[parent].arrival;
		if (val_current < val_parent) {
			current = parent;
		}
	}
	if (current != i) {
		swap2dPoints(&in_Process[current], &in_Process[i]);
		heapifyUp2D(in_Process, current);
	}

}

void heapifyVector2D(std::vector<pointFastMarching2D>& in_Process) {
	int length_array = in_Process.size();
	int indx, start = length_array / 2 - 1;
	for (indx = start; indx >= 0; indx--) {
		heapifyDown2D(in_Process, indx);
	}
}

void deleteRootHeap2D(std::vector<pointFastMarching2D>& in_Process) {
	int l = in_Process.size();
	swap2dPoints(&in_Process[0], &in_Process[l - 1]);
	in_Process.pop_back();
	heapifyDown2D(in_Process, 0);
}

void addPointHeap2D(std::vector<pointFastMarching2D>& in_Process, pointFastMarching2D point) {
	in_Process.push_back(point);
	int l = in_Process.size();
	heapifyUp2D(in_Process, l - 1);
}

int getIndexFromHeap2D(std::vector<pointFastMarching2D>& in_Process, size_t i, size_t j) {
	for (int ind = 0; ind < in_Process.size(); ind++) {
		if (in_Process[ind].x == i && in_Process[ind].y == j) {
			return ind;
		}
	}
	return -1; //not found
}

dataType selectX(dataType* actionMapPtr, const size_t length, const size_t width, const size_t ind_x, const size_t ind_y) {
	dataType x_minus, x_plus;
	if (ind_x == 0) {
		x_minus = INFINITY;
	}
	else {
		x_minus = actionMapPtr[x_new(ind_x - 1, ind_y, length)];
	}
	if (ind_x == length - 1) {
		x_plus = INFINITY;
	}
	else {
		x_plus = actionMapPtr[x_new(ind_x + 1, ind_y, length)];
	}
	return std::min(x_minus, x_plus);
}

dataType selectY(dataType* actionMapPtr, const size_t length, const size_t width, const size_t ind_x, const size_t ind_y) {
	dataType y_minus, y_plus;
	if (ind_y == 0) {
		y_minus = INFINITY;
	}
	else {
		y_minus = actionMapPtr[x_new(ind_x, ind_y - 1, length)];
	}
	if (ind_y == width - 1) {
		y_plus = INFINITY;
	}
	else {
		y_plus = actionMapPtr[x_new(ind_x, ind_y + 1, length)];
	}
	return std::min(y_minus, y_plus);
}

dataType solve2dQuadratic(dataType X, dataType Y, dataType P, PixelSpacing h) {

	dataType solution = 0.0, a = 0.0, b = 0.0, c = 0.0, delta = 0.0;
	dataType P_2 = P * P;

	dataType hx = h.sx;
	dataType hx_2 = h.sx * h.sx;

	dataType hy = h.sy;
	dataType hy_2 = h.sy * h.sy;

	if (P <= 0.0) {
		std::cerr << "Error: P must be positive." << std::endl;
		return solution; // Return 0 if P is not positive
	}

	if (X == INFINITY && Y != INFINITY)
	{
		a = 1.0;
		b = -2 * Y;
		c = (dataType)(Y * Y - hy_2 * P_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			solution = (dataType)((-b + std::sqrt(delta)) / (2 * a));
			if (solution >= Y) {
				return solution;
			}
			else {
				return (dataType)(Y + hy * P);
			}
		}
		else {
			return (dataType)(Y + hy * P);
		}
	}

	if (X != INFINITY && Y == INFINITY)
	{
		a = 1.0;
		b = -2 * X;
		c = (dataType)(X * X - hx_2 * P_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0)
		{
			solution = (dataType)((-b + std::sqrt(delta)) / (2 * a));
			if (solution >= X)
			{
				return solution;
			}
			else {
				return (dataType)(X + hx * P);
			}
		}
		else {
			return (dataType)(X + hx * P);
		}
	}

	if (X != INFINITY && Y != INFINITY)
	{
		a = hx_2 + hy_2;
		b = -2 * (hy_2 * X + hx_2 * Y);
		c = (dataType)(hy_2 * X * X + hx_2 * Y * Y - hx_2 * hy_2 * P_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0)
		{
			solution = (dataType)((-b + std::sqrt(delta)) / (2 * a));
			if (solution >= std::max(X, Y)) {
				return solution;
			}
			else {
				return (dataType)(std::min(X + hx * P, Y + hy * P));
			}
		}
		else {
			return (dataType)(std::min(X + hx * P, Y + hy * P));
		}
	}

}

bool frontPropagation2D(Image_Data2D imageData, dataType* actionMapPtr, dataType* potentialPtr, Point2D* seedPoints) {

	const size_t length = imageData.height;
	const size_t width = imageData.width;
	PixelSpacing spacing = imageData.spacing;
	size_t length_minus = length - 1;
	size_t width_minus = width - 1;
	size_t dim2D = length * width;

	std::vector<pointFastMarching2D> inProcess;
	short* labelArray = new short[dim2D] { 0 };
	if (imageData.imageDataPtr == NULL || actionMapPtr == NULL || potentialPtr == NULL || seedPoints == NULL || labelArray == NULL) {
		return false;
	}

	//1 ---> already processed, 
	//2 ---> in process and 
	//3 ---> not processed
	for (size_t k = 0; k < dim2D; k++) {
		actionMapPtr[k] = INFINITY;
		labelArray[k] = 3;
	}

	size_t i = (size_t)seedPoints[0].x;
	size_t j = (size_t)seedPoints[0].y;
	size_t currentIndx = x_new(i, j, length);
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
		inProcess.push_back(EastNeighbor);
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
		inProcess.push_back(WestNeighbor);
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
		inProcess.push_back(NorthNeighbor);
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
		inProcess.push_back(SouthNeighbor);
		labelArray[indxSouth] = 2;
	}

	//heapify 2D vector
	heapifyVector2D(inProcess);

	pointFastMarching2D current;
	short label = 0;

	while (inProcess.size() > 0) {

		current = inProcess[0];
		i = current.x;
		j = current.y;
		currentIndx = x_new(i, j, length);
		labelArray[currentIndx] = 1;
		deleteRootHeap2D(inProcess);

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
					addPointHeap2D(inProcess, WestNeighbor);
				}
				else {
					if (dWest < actionMapPtr[indxWest]) {
						actionMapPtr[indxWest] = dWest;
						size_t pIndex = getIndexFromHeap2D(inProcess, iminus, j);
						if (pIndex != -1) {
							heapifyUp2D(inProcess, pIndex);
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
					addPointHeap2D(inProcess, EastNeighbor);
				}
				else {
					if (dEast < actionMapPtr[indxEast]) {
						actionMapPtr[indxEast] = dEast;
						size_t pIndex = getIndexFromHeap2D(inProcess, iplus, j);
						if (pIndex != -1) {
							heapifyUp2D(inProcess, pIndex);
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
					addPointHeap2D(inProcess, NorthNeighbor);
				}
				else {
					if (dNorth < actionMapPtr[indxNorth]) {
						actionMapPtr[indxNorth] = dNorth;
						size_t pIndex = getIndexFromHeap2D(inProcess, i, jminus);
						if (pIndex != -1) {
							heapifyUp2D(inProcess, pIndex);
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
					addPointHeap2D(inProcess, SouthNeighbor);
				}
				else {
					if (dSouth < actionMapPtr[indxSouth]) {
						actionMapPtr[indxSouth] = dSouth;
						size_t pIndex = getIndexFromHeap2D(inProcess, i, jplus);
						if (pIndex != -1) {
							heapifyUp2D(inProcess, pIndex);
						}
					}
				}
			}
		}

	}

	inProcess.clear();
	delete[] labelArray;
}
