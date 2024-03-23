#include <memory.h>
#include "common_functions.h"
#include <math.h>
#include <stdlib.h>


//==============================================================================
// Local Function Prototype
//==============================================================================
void reflection3DB(dataType** toReflectImage, size_t imageHeight, size_t imageLength, size_t imageWidth, size_t p)
{
	size_t k, i, j;
	size_t height = imageHeight, length = imageLength, width = imageWidth; // Actual Z, X, Y dimensions
	size_t rowLength = length + 2 * p;
	// Modified Dimension less 1
	size_t mheight = height - 1, mlength = length - 1, mwidth = width - 1;
	// Z Reflection
	for (i = p; i <= (mlength)+p; i++)
	{
		for (j = p; j <= (mwidth)+p; j++)
		{
			// If
			toReflectImage[p - 1][x_new(i, j, rowLength)] = toReflectImage[p + 1][x_new(i, j, rowLength)];
			toReflectImage[(mheight)+p + 1][x_new(i, j, rowLength)] = toReflectImage[(mheight)+p - 1][x_new(i, j, rowLength)];
		}
	}
	// Y reflection
	for (j = p; j <= (mwidth)+p; j++)
	{
		for (k = 0; k <= (mheight)+2 * p; k++)
		{
			toReflectImage[k][x_new(p - 1, j, rowLength)] = toReflectImage[k][x_new(p + 1, j, rowLength)];
			toReflectImage[k][x_new((mlength)+p + 1, j, rowLength)] = toReflectImage[k][x_new((mlength)+p - 1, j, rowLength)];
		}
	}
	// X Direction
	for (i = 0; i <= (mlength)+2 * p; i++)
	{
		for (k = 0; k <= (mheight)+2 * p; k++)
		{
			toReflectImage[k][x_new(i, p - 1, rowLength)] = toReflectImage[k][x_new(i, p + 1, rowLength)];
			toReflectImage[k][x_new(i, (mwidth)+p + 1, rowLength)] = toReflectImage[k][x_new(i, (mwidth)+p - 1, rowLength)];
		}
	}
}
//==============================================================================
//imageHeight and imageLength is considered as data extension together with boundary
void reflection3D(dataType** toReflectImage, size_t imageHeight, size_t imageLength, size_t imageWidth)
{
	size_t k, i, j;
	size_t length = imageLength, width = imageWidth; // Actual X, Y dimensions

	const size_t heightMin = imageHeight - 1;
	const size_t widthMin = width - 1;
	const size_t lengthMin = length - 1;
	const size_t sliceSize = imageLength * imageWidth;

	size_t x;

	// Y reflection
	for (k = 1; k <= heightMin; k++)
	{
		for (i = 0; i < length; i++)
		{
			toReflectImage[k][x_new(i, 0, length)] = toReflectImage[k][x_new(i, 1, length)];
			toReflectImage[k][x_new(i, widthMin, length)] = toReflectImage[k][x_new(i, widthMin - 1, length)];
		}
	}

	// X Direction
	for (k = 1; k <= heightMin; k++)
	{
		for (j = 0; j < width; j++)
		{
			x = x_new(0, j, length);
			toReflectImage[k][x] = toReflectImage[k][x + 1];

			x = x_new(lengthMin, j, length);
			toReflectImage[k][x] = toReflectImage[k][x - 1];
		}
	}

	// Z Reflection
	for (i = 0; i < sliceSize; i++)
	{
		toReflectImage[0][i] = toReflectImage[1][i];
		toReflectImage[heightMin][i] = toReflectImage[heightMin - 1][i];
	}
}
//==============================================================================
/*
* Evaluates the diffusivity value
*/
dataType edgeDetector(dataType value, dataType coef)
{
	return (dataType)(1.0 / (1 + coef * value));
}

dataType similarIntensityDetector(dataType currValue, dataType refValue, dataType coef)
{
	return (dataType)(1.0 / (1 + coef * pow(currValue - refValue,2)));
}


//==============================================================================
size_t x_new(const size_t rowIndex, const size_t columnIndex, const size_t rowLength)
{
	return rowIndex + columnIndex * rowLength; // x + y*DimX
}
//==============================================================================
void copyDataToExtendedArea(dataType** originalDataPtr, dataType** extendedDataPtr, const size_t originalHeight, const size_t originalLength, const size_t originalWidth)
{
	const size_t height_ext = originalHeight + 2;
	const size_t length_ext = originalLength + 2;
	const size_t width_ext = originalWidth + 2;

	size_t sliceBound = (length_ext - 1) * width_ext;
	size_t i, k, k_ext, i_d = 0;

	for (k = 0, k_ext = 1; k < originalHeight; k++, k_ext++)
	{
		i_d = 0;
		for (i = length_ext + 1; i < sliceBound; i += length_ext)
		{
			memcpy(&(extendedDataPtr[k_ext][i]), &(originalDataPtr[k][i_d]), originalLength * sizeof(dataType));
			i_d += originalLength;
		}
	}
}
//==============================================================================
void copyDataToReducedArea(dataType** originalDataPtr, const dataType** extendedDataPtr, const size_t originalHeight, const size_t originalLength, const size_t originalWidth)
{
	const size_t height_ext = originalHeight + 2;
	const size_t length_ext = originalLength + 2;
	const size_t width_ext = originalWidth + 2;

	size_t sliceBound = (length_ext - 1) * width_ext;
	size_t i, k, k_ext, i_d = 0;

	for (k = 0, k_ext = 1; k < originalHeight; k++, k_ext++)
	{
		i_d = 0;
		for (i = length_ext + 1; i < sliceBound; i += length_ext)
		{
			memcpy(&(originalDataPtr[k][i_d]), &(extendedDataPtr[k_ext][i]), originalLength * sizeof(dataType));
			i_d += originalLength;
		}
	}
}
//==============================================================================
void copyDataToAnotherArray(dataType** source, dataType** destination, size_t height, size_t length, size_t width)
{
	size_t k, i, j, xd;
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				// 2D flattening
				xd = x_new(i, j, length);
				destination[k][xd] = source[k][xd];
			}
		}
	}
}
//==============================================================================
size_t x_flat(const size_t rowIndex, const size_t columnIndex, const size_t heightIndex, const size_t rowLength, const size_t columnLength)
{
	return rowIndex + rowLength * (columnIndex + columnLength * heightIndex); // x + xlen * (y + ylen * z)
}
//==============================================================================
void rescaleNewRange(dataType** imageDataPtr, size_t imageLength, size_t imageWidth, size_t imageHeight, dataType minNew, dataType maxNew, dataType max_dta, dataType min_dta) {
	size_t k, i, j, xd;

	// Find the Min and Max Intensity
	for (k = 0; k < imageHeight; k++) {
		for (i = 0; i < imageLength; i++) {
			for (j = 0; j < imageWidth; j++) {

				// 1D Conversion of row and column
				xd = x_new(i, j, imageLength);
				if (imageDataPtr[k][xd] > max_dta) {
					max_dta = imageDataPtr[k][xd];
				}
				if (imageDataPtr[k][xd] < min_dta) {
					min_dta = imageDataPtr[k][xd];
				}
			}
		}
	}
	// Rescale from min_new to max_new
	dataType diffOld = max_dta - min_dta;
	dataType diffNew = maxNew - minNew;
	dataType scale_factor = (diffNew) / (diffOld);
	for (k = 0; k < imageHeight; k++) {
		for (i = 0; i < imageLength; i++) {
			for (j = 0; j < imageWidth; j++) {
				// 1D Conversion of row and column
				xd = x_new(i, j, imageLength);
				imageDataPtr[k][xd] = scale_factor * (imageDataPtr[k][xd] - max_dta) + maxNew;
				// Alternatively
				//dta[k][xd] = scale_factor * (dta[k][xd] - min_dta) + min_new; }
			}
		}
	}
}
//==============================================================================

//2D function
void copyDataToAnother2dArray(dataType* source, dataType* destination, size_t imageHeight, size_t imageWidth) {
	size_t i, j, xd;
	for (i = 0; i < imageHeight; i++) {
		for (j = 0; j < imageWidth; j++) {
			xd = x_new(i, j, imageHeight);
			destination[xd] = source[xd];
		}
	}
}
//==============================================================================
void copyDataTo2dExtendedArea(dataType* originalDataPtr, dataType* extendedDataPtr, const size_t originalHeight, const size_t originalWidth)
{
	const size_t height_ext = originalHeight + 2;
	const size_t width_ext = originalWidth + 2;

	size_t sliceBound = (height_ext - 1) * width_ext;
	size_t i, j, i_ext, j_ext, i_d = 0;

	for (i = 0, i_ext = 1; i < originalHeight; i++, i_ext++) {
		for (j = 0, j_ext = 1; j < originalWidth; j++, j_ext++) {
			extendedDataPtr[x_new(i_ext, j_ext, height_ext)] = originalDataPtr[x_new(i, j, originalHeight)];
		}
	}

}
//==============================================================================
void copyDataTo2dReducedArea(dataType* originalDataPtr, const dataType* extendedDataPtr, const size_t originalHeight, const size_t originalWidth)
{
	const size_t height_ext = originalHeight + 2;
	const size_t width_ext = originalWidth + 2;

	size_t sliceBound = (height_ext - 1) * width_ext;
	size_t i, j, i_ext, j_ext, i_d = 0;

	for (i = 0, i_ext = 1; i < originalHeight; i++, i_ext++) {
		for (j = 0, j_ext = 1; j < originalWidth; j++, j_ext++) {
			originalDataPtr[x_new(i, j, originalHeight)] = extendedDataPtr[x_new(i_ext, j_ext, height_ext)];
		}
	}

}
//==============================================================================
void reflection2D(dataType* toReflectImage, size_t imageHeight, size_t imageWidth)
{
	size_t i, j;
	size_t height = imageHeight, width = imageWidth;

	const size_t heightMin = height - 1;
	const size_t widthMin = width - 1;

	size_t x;

	// Y reflection
	for (i = 0; i < height; i++)
	{
		toReflectImage[x_new(i, 0, height)] = toReflectImage[x_new(i, 1, height)];
		toReflectImage[x_new(i, widthMin, height)] = toReflectImage[x_new(i, widthMin - 1, height)];
	}

	// X Direction
	for (j = 0; j < width; j++)
	{
		x = x_new(0, j, height);
		toReflectImage[x] = toReflectImage[x + 1];

		x = x_new(heightMin, j, height);
		toReflectImage[x] = toReflectImage[x - 1];
	}

}
//==============================================================================
double getPoint2DDistance(const Point2D a, const Point2D b)
{
	return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}
//==============================================================================
double getPoint3DDistance(const Point3D a, const Point3D b)
{
	return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}
//==============================================================================
Point3D getPointWithTheHighestValue(dataType** distanceMapPtr, const size_t length, const size_t width, const size_t height) {

	Point3D result = { 0.0, 0.0, 0.0 };
	dataType max_value = 0.0;

	for (size_t k = 0; k < height; k++) {
		for (size_t i = 0; i < length; i++) {
			for (size_t j = 0; j < width; j++) {
				size_t x = x_new(i, j, length);
				if (distanceMapPtr[k][x] > max_value) {
					max_value = distanceMapPtr[k][x];
					result.x = (dataType)i; 
					result.y = (dataType)j; 
					result.z = (dataType)k;
				}
			}
		}
	}
	return result;
}

bool copyCurve2DPointsToArray(const Curve2D* pCurve, dataType** pArray)
{
	if (pCurve == NULL || pArray == NULL) {
		return false;
	}

	for (int i = 0; i < pCurve->numPoints; i++) {
		pArray[i][0] = pCurve->pPoints[i].x;
		pArray[i][1] = pCurve->pPoints[i].y;
	}

	return true;
}

bool isCurveClosed(const Curve2D* pcurve)
{
	size_t num_end_points = 0;
	for (size_t i = 0; i < pcurve->numPoints; i++)
	{
		if (pcurve->pPoints[i].isEndPoint)
		{
			num_end_points++;
		}
	}

	if (num_end_points == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool isCurveOrientedPositively(const Curve2D* pcurve)
{
	if (pcurve == NULL) {
		return false;
	}

	double signed_area = 0;
	double x_i_plus = -1;
	double y_i_plus = -1;

	for (size_t i = 0; i < pcurve->numPoints; i++)
	{
		if (i + 1 == pcurve->numPoints)
		{
			x_i_plus = pcurve->pPoints[0].x;
			y_i_plus = pcurve->pPoints[0].y;
		}
		else
		{
			x_i_plus = pcurve->pPoints[i+1].x;
			y_i_plus = pcurve->pPoints[i + 1].y;
		}

		signed_area += (pcurve->pPoints[i].x * y_i_plus - x_i_plus * pcurve->pPoints[i].y);
	}

	signed_area *= 0.5;

	return signum(signed_area) > 0.0;
}

Point2D getCurveCentroid(const Curve2D* pcurve)
{
	if (pcurve == NULL)
	{
		return (Point2D){ 0.0, 0.0};
	}
	dataType x = 0, y = 0;

	for (size_t i = 0; i < pcurve->numPoints; i++) {
		x += pcurve->pPoints[i].x;
		y += pcurve->pPoints[i].y;
	}

	return (Point2D){x/ pcurve->numPoints, y/ pcurve->numPoints	};
}

bool getGradient2D(dataType* pbase_data, const size_t width, const size_t height, const size_t ind_x, const size_t ind_y, const FiniteVolumeSize2D sz, Point2D* grad)
{
	if (pbase_data == NULL || height < 2 || width < 2 ||
		sz.hx == 0 || sz.hy == 0) {
		return false;
	}

	size_t x = ind_x, y = ind_y;
	if (x >= width)
	{
		x = width - 1;
	}

	if (y >= height)
	{
		y = height - 1;
	}

	dataType dx = 0, dy = 0;
	dataType hx_c = 2 * sz.hx;
	dataType hy_c = 2 * sz.hy;

	if (x == 0)
	{
		dx = (pbase_data[x_new(x + 1, y, width)] - pbase_data[x_new(x, y, width)]) / sz.hx;
	}
	else if (x == width - 1)
	{
		dx = (pbase_data[x_new(x, y, width)] - pbase_data[x_new(x - 1, y, width)]) / sz.hx;
	}
	else
	{
		dx = (pbase_data[x_new(x + 1, y, width)] - pbase_data[x_new(x - 1, y, width)]) / hx_c;
	}

	if (y == 0)
	{
		dy = (pbase_data[x_new(x, y + 1, width)] - pbase_data[x_new(x, y, width)]) / sz.hy;
	}
	else if (y == height - 1)
	{
		dy = (pbase_data[x_new(x, y, width)] - pbase_data[x_new(x, y - 1, width)]) / sz.hy;
	}
	else
	{
		const size_t xtmp = x_new(x, y + 1, width);
		dy = (pbase_data[x_new(x, y + 1, width)] - pbase_data[x_new(x, y - 1, width)]) / hy_c;
	}


	grad->x = dx;
	grad->y = dy;

	return true;
}

dataType norm(const Point2D pt)
{
	return (dataType)sqrt(pt.x * pt.x + pt.y * pt.y);
}

double signum(const double value)
{
	if (value > 0.0)
		return 1.0;
	else if (value == 0.0)
		return 0.0;
	else
		return -1.0;
}

//LinkedCurve
LinkedCurve createLinkedCurve()
{
	LinkedCurve linkedCurve;
	linkedCurve.number_of_points = 0;
	//linkedCurve.first_point = CreateLinkedPoint(first_point_x, first_point_y);
	linkedCurve.length = 0;
	return linkedCurve;
}

LinkedPoint* createLinkedPoint(const double point_x, const double point_y)
{
	LinkedPoint* linked_point = (LinkedPoint*)malloc(sizeof(LinkedPoint));
	linked_point->x = point_x;
	linked_point->y = point_y;
	linked_point->next = NULL;
	linked_point->previous = NULL;
	linked_point->distance_to_next = 0;
	linked_point->id = getNextID();
	return linked_point;
}

bool initializeLinkedCurve(Curve2D * pcurve, LinkedCurve * plinked_curve, const bool reverse, const bool close_curve)
{
	if (pcurve == NULL || plinked_curve == NULL)
	{
		return false;
	}

	size_t from = 0, to = pcurve->numPoints - 1;
	size_t increment = 1;
	if (reverse)
	{
		from = pcurve->numPoints;
		to = 0;
		size_t increment = -1;
	}

	size_t i = from;
	LinkedPoint * pcurrent_point = createLinkedPoint(pcurve->pPoints[from].x, pcurve->pPoints[from].y);
	plinked_curve->first_point = pcurrent_point;
	plinked_curve->number_of_points = 1;

	while (i != to)
	{
		i += increment;

		pcurrent_point = pushAfterPoint(plinked_curve, pcurrent_point, pcurve->pPoints[i].x, pcurve->pPoints[i].y);
	}

	pcurrent_point->next = plinked_curve->first_point;
	plinked_curve->first_point->previous = pcurrent_point;
	updateDistanceToNext(plinked_curve, pcurrent_point);

	return true;
}

LinkedPoint* pushAfterPoint(LinkedCurve* linked_curve, LinkedPoint* linked_point, const double point_x, const double point_y)
{
	LinkedPoint* new_linked_point = createLinkedPoint(point_x, point_y);

	if (linked_point->next != NULL)
	{
		new_linked_point->next = linked_point->next;
		(linked_point->next)->previous = new_linked_point;
	}
	else if (linked_curve->first_point->previous == NULL)
	{
		linked_curve->first_point->previous = new_linked_point;
		new_linked_point->next = linked_curve->first_point;
	}

	linked_point->next = new_linked_point;
	new_linked_point->previous = linked_point;

	linked_curve->number_of_points++;

	updateDistanceToNext(linked_curve, linked_point);
	updateDistanceToNext(linked_curve, new_linked_point);

	return new_linked_point;
}

void releaseLinkedCurve(LinkedCurve* linked_curve)
{
	if (linked_curve == NULL)
	{
		return;
	}

	LinkedPoint* current_point = linked_curve->first_point;
	if (current_point != NULL &&
		current_point->previous != NULL)
	{
		(current_point->previous)->next = NULL;
	}

	LinkedPoint* next_point = NULL;

	while (current_point != NULL) {
		next_point = current_point->next;
		free(current_point);
		current_point = NULL;
		linked_curve->number_of_points--;

		current_point = next_point;
	}
}

double updateDistanceToNext(LinkedCurve* linked_curve, LinkedPoint* linked_point)
{
	if (linked_curve == NULL ||
		linked_point == NULL ||
		linked_point->next == NULL)
	{
		return -1;
	}

	const double old_distance = linked_point->distance_to_next;

	double hx = fabs(linked_point->x - linked_point->next->x);
	double hy = fabs(linked_point->y - linked_point->next->y);

	linked_point->distance_to_next = sqrt(hx * hx + hy * hy);
	linked_point->distance_to_next_x = hx;
	linked_point->distance_to_next_y = hy;

	double diff = linked_point->distance_to_next - old_distance;

	linked_curve->length = max(linked_curve->length + diff, 0);

	return linked_point->distance_to_next;
}

bool updateLinkedCurveLengths(LinkedCurve* linked_curve)
{
	if (linked_curve == NULL)
	{
		return false;
	}

	LinkedPoint* current_point = linked_curve->first_point;
	const unsigned long long first_id = current_point->id;
	double curve_length = 0;
	size_t number_of_points = 0;

	do {
		number_of_points++;
		curve_length += updateDistanceToNext(linked_curve, current_point);
		current_point = current_point->next;
	} while (current_point->id != first_id);

	linked_curve->number_of_points = number_of_points;
	linked_curve->length = curve_length;

	return true;
}

bool updatePoint(LinkedCurve* linked_curve, LinkedPoint* linked_point, const double x, const double y)
{
	if (linked_point == NULL)
	{
		return false;
	}

	if (linked_point->x == x &&
		linked_point->y == y) {
		return true;
	}

	linked_point->x = x;
	linked_point->y = y;

	if (linked_point->previous != NULL)
	{
		updateDistanceToNext(linked_curve, linked_point);
		updateDistanceToNext(linked_curve, linked_point->previous);
	}

	return true;
}

static unsigned long long id = 1;

unsigned long long getNextID()
{
	return id++;
}

void resetIDGenerator()
{
	id = 1;
}