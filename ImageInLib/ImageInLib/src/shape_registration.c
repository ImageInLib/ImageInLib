//==============================================================================
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
//==============================================================================
#include "shape_registration.h"
//==============================================================================
size_t bandPts = 0, clipBoxPts = 0;
//==============================================================================
// Stores generated random points
typedef struct {
	size_t k, i, j;
} Random3DPoints;
// stores calc. distances from generated points
typedef struct {
	dataType distDifference;
} distanceCalculated;
//
typedef struct {
	dataType pvl, qvl, hvl;
} xNbDistances;

typedef struct {
	dataType pvl, qvl, hvl;
} yNbDistances;

typedef struct {
	dataType pvl, qvl, hvl;
} zNbDistances;

typedef struct {
	dataType xFwd, yFwd, zFwd;
} Finite_Differences;
//==============================================================================
typedef struct { size_t k, i, j; } CoordPoints;
//==============================================================================
CoordPoints transformPoint(CoordPoints * inputPoints, Point3D translation, Point3D scaling, Point3D rotation, dataType centroid[3], size_t imageHeight, size_t imageLength, size_t imageWidth, int loc);
dataType getDistance(dataType ** binaryImage, size_t imageHeight, size_t imageLength, size_t dim2D, const size_t k1, const size_t x1, ClipBox bestfitBox, Point3D * surface_points, size_t ptsNum, dataType insideShapevalue, bool parallelize);
size_t surfacePoints(dataType ** binaryImage, size_t imageLength, const unsigned char fgroundValue, ClipBox bestfitBox);
//==============================================================================
void nbPointsX(dataType ** transformedBinaryData, dataType pixelSize, size_t x, size_t k, size_t i, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType * pValue, dataType * qValue, dataType *hh, ClipBox bestfitBox, Point3D * surface_points, size_t ptsNum, dataType insideShapevalue, bool parallelize);
void nbPointsY(dataType ** transformedBinaryData, dataType pixelSize, size_t k, size_t i, size_t j, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType * pValue, dataType * qValue, dataType *hh, ClipBox bestfitBox, Point3D * surface_points, size_t ptsNum, dataType insideShapevalue, bool parallelize);
void nbPointsZ(dataType ** transformedBinaryData, dataType pixelSize, size_t x, size_t k, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType * pValue, dataType * qValue, dataType *hh, ClipBox bestfitBox, Point3D * surface_points, size_t ptsNum, dataType insideShapevalue, bool parallelize);
//==============================================================================
dataType finiteDifAll(dataType pValue, dataType qValue, dataType hh);
//==============================================================================
void run_registration(dataType **fixedData, dataType **movingData, dataType **resultPtr, size_t zDim, size_t xDim, size_t yDim, Registration_Params params, Optimization_Method gdescentMethod)
{
	//==============================================================================
	// Variable definition
	dataType  step_size = params.step_size;
	dataType tol = params.tolerance;
	//==============================================================================
	dataType firstCpuTime, secondCpuTime;
	//==============================================================================
	// Centroid Parameters
	dataType fixedCentroid[3], movingCentroid[3];
	// Centroid Method Fixed/Destination Data
	centroidImage(fixedData, fixedCentroid, zDim, xDim, yDim, params.imageBackground);
	// Centroid Method Moving/Source Data
	centroidImage(movingData, movingCentroid, zDim, xDim, yDim, params.imageBackground);
	//==============================================================================
	// Set the Translation approximation
	Point3D translationTran;
	translationTran.x = fixedCentroid[0] - movingCentroid[0];
	translationTran.y = fixedCentroid[1] - movingCentroid[1];
	translationTran.z = fixedCentroid[2] - movingCentroid[2];
	//==============================================================================
	// Sets Transformation Parameters to be used in Registration
	Affine_Parameter finalResults;
	Point3D rotationTran = { 0.0, 0.0, 0.0 };
	Point3D scalingTran = { 1.0, 1.0, 1.0 };
	finalResults.rotation = rotationTran, finalResults.scaling = scalingTran, finalResults.translation = translationTran;
	//==============================================================================
	// Call the registration Function
	//==============================================================================
	// Begin Record Time
#ifdef MEASURE_TIME
	firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
	//==============================================================================
	if (gdescentMethod == GRADIENT_DESCENT)
	{
		finalResults = registration3D(fixedData, movingData, finalResults, step_size, tol, zDim, xDim, yDim, movingCentroid, params);
	}
	else if (gdescentMethod == STOCHASTIC_DESCENT)
	{
		finalResults = registrationStochastic3D(fixedData, movingData, finalResults, step_size, tol, zDim, xDim, yDim, movingCentroid, params);
	}
	else if (gdescentMethod == BLOCK_COORDINATE_DESCENT)
	{
		finalResults = registrationCoorDinateDescent3D(fixedData, movingData, finalResults, step_size, tol, zDim, xDim, yDim, movingCentroid, params);
	}
	//==============================================================================
#ifdef MEASURE_TIME
	secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
	//==============================================================================
	printf("Total Registration CPU time: %e secs\n", secondCpuTime - firstCpuTime);
	//==============================================================================
	// Apply transformation results to destination - Expect same as source~approximately
	// Begin Record Time
#ifdef MEASURE_TIME
	firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
	//==============================================================================
	transform3DImage(movingData, resultPtr, finalResults.translation, finalResults.scaling, finalResults.rotation, zDim, xDim, yDim, params.imageBackground, movingCentroid, params.parallelize);
	//==============================================================================
#ifdef MEASURE_TIME
	secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
	//==============================================================================
#endif
#ifdef CONSOLE_OUTPUT
	printf("Final Resulting Transformation CPU time: %e secs\n\n", secondCpuTime - firstCpuTime);
#endif
	//==============================================================================
	// Save the Resultant Transformation
	params.affineResults = finalResults;
	//==============================================================================
}
//==============================================================================
void fastMarching(dataType ** distancePtr, dataType ** dataSourcePtr, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType objPixel)
{
	size_t k, i, j, x;

	struct Node * band = NULL; // Holds all the Objects
							   // Sets the structure size, to hold all the calculated arrival times
	Obj_Structure ** objectNthD = (Obj_Structure **)malloc(sizeof(Obj_Structure*) * imageHeight);
	for (i = 0; i < imageHeight; i++)
	{
		objectNthD[i] = (Obj_Structure *)malloc(sizeof(Obj_Structure) * (imageLength * imageWidth));
	}
	// Initialize Object2D
	for (k = 0; k < imageHeight; k++)
	{
		for (i = 0; i < imageLength; i++)
		{
			for (j = 0; j < imageWidth; j++)
			{
				// 2D to 1D representation for i, j
				size_t x = x_new(i, j, imageLength);
				objectNthD[k][x].arrival = INFINITY;
				objectNthD[k][x].state = UNKNOWN;
				objectNthD[k][x].xpos = i;
				objectNthD[k][x].ypos = j;
				objectNthD[k][x].zpos = k;
				objectNthD[k][x].position = x_flat(i, j, k, imageLength, imageWidth);
			}
		}
	}
	Point3D *shapePoints = (Point3D *)malloc(sizeof(Point3D)*(imageHeight * imageLength * imageWidth));
	int loop = 0;
	// Derive the points

	for (k = 0; k < imageHeight; k++)
	{
		for (i = 0; i < imageLength; i++)
		{
			for (j = 0; j < imageWidth; j++)
			{
				// 1D representation
				x = x_new(i, j, imageLength);
				if (dataSourcePtr[k][x] == objPixel) // Fill value for block
				{
					//Save the dimension with those values
					shapePoints[loop].x = (dataType)i;
					shapePoints[loop].y = (dataType)j;
					shapePoints[loop].z = (dataType)k;
					loop++;
				}
			}
		}
	}

	// Arrival times
	Arrival_Time *shapeArrival = (Arrival_Time *)malloc(sizeof(Arrival_Time)*loop);
	for (i = 0; i < loop; i++)
	{
		shapeArrival[i].T = 0.0;
	}

	// Calls Fm3D
	fastMarching3D(band, objectNthD, shapePoints, shapeArrival, imageHeight, imageLength, imageWidth, loop);
	// Copy Fast marching modified to distancePtr
	for (k = 0; k < imageHeight; k++)
	{
		for (i = 0; i < (imageLength*imageWidth); i++)
		{
			distancePtr[k][i] = objectNthD[k][i].arrival;
		}
	}
	shapePoints = NULL;
	shapeArrival = NULL;
	for (i = 0; i < imageHeight; i++)
	{
		free(objectNthD[i]);
	}
	//deleteList(band);
}
//==============================================================================
// Prototypes
void centroidImage(dataType ** imageDataPtr, dataType *centroid, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType imageBackground)
{
	size_t k, i, j, counts = 0;
	dataType x = 0.0, y = 0.0, z = 0.0;
	for (k = 0; k < imageHeight; k++)
	{
		for (i = 0; i < imageLength; i++)
		{
			for (j = 0; j < imageWidth; j++)
			{
				// 2D to 1D representation for i, j
				size_t xD = x_new(i, j, imageLength);
				if (imageDataPtr[k][xD] != imageBackground)
				{
					x += i;
					y += j;
					z += k;
					counts++;
				}
			}
		}
	}
	// Check K incase no 0 was found - not a shape
	if (counts == 0)
	{
		counts = 1;
	}
	// Set the Centers of the shape
	centroid[0] = x / counts, centroid[1] = y / counts, centroid[2] = z / counts;
}
//==============================================================================
void centroidClipBox(dataType *centroid, ClipBox coord, dataType ** imageDataPtr, size_t imageLength, dataType imageBackground)
{
	size_t k, i, j, counts = 0;
	dataType x = 0.0, y = 0.0, z = 0.0;
	for (k = coord.k_min; k <= coord.k_max; k++)
	{
		for (i = coord.i_min; i <= coord.i_max; i++)
		{
			for (j = coord.j_min; j <= coord.j_max; j++)
			{
				// 2D to 1D representation for i, j
				size_t xD = x_new(i, j, imageLength);
				if (imageDataPtr[k][xD] != imageBackground)
				{
					x += i;
					y += j;
					z += k;
					counts++;
				}
			}
		}
	}
	// Check K incase no 0 was found - not a shape
	if (counts == 0)
	{
		counts = 1;
	}
	// Set the Centers of the shape
	centroid[0] = x / counts, centroid[1] = y / counts, centroid[2] = z / counts;
}
//==============================================================================
inline int NFunctionBinary(dataType v1, dataType v2, dataType delta)
{
	//==============================================================================
	if (abs((int)v1) == delta || abs((int)v2) == delta)
	{
		return 1;
	}
	else
	{
		return 0;
	}
	//==============================================================================
}
//==============================================================================
int NFunction(dataType val1, dataType val2, dataType delta)
{
	if (fabs(val1) < fabs(val2))
	{
		if (fabs(val1) > delta)
		{
			return 0;
		}
		else
		{
			return 1;
		}
	}
	else
	{
		if (fabs(val2) > delta)
		{
			return 0;
		}
		else
		{
			return 1;
		}
	}
}
inline int NFunctionOne(dataType v1, dataType delta)
{
#ifdef USE_NARROWBAND
	if (abs((int)v1) > delta)
	{
		return 0;
	}
	else
	{
		return 1;
	}
#else
	return 1;
#endif // USE_NARROWBAND
}
//==============================================================================
dataType energyFunction(dataType ** destination, dataType ** distTrans, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType h)
{
	size_t k, i, j, counter = 0;
	dataType energy = 0.0;
	// Calcaulate the error in grid functions
	for (k = 0; k < imageHeight; k++)
	{
		for (i = 0; i < imageLength; i++)
		{
			for (j = 0; j < imageWidth; j++)
			{
				// 2D to 1D representation for i, j
				size_t x = x_new(i, j, imageLength);
				if (NFunction(destination[k][x], distTrans[k][x], NDelta) == 1)
				{
					energy += (destination[k][x] - distTrans[k][x])   * (destination[k][x] - distTrans[k][x]);
					counter++;
				}
			}
		}
	}
	return ((h * energy) / (2 * counter));
}
//==============================================================================
dataType energyFunctionClip(dataType ** destination, dataType **distTrans, ClipBox coord, size_t imageLength)
{
	size_t k, i, j, counter = 0, clipPts = 0;
	dataType energy = 0.0;
	//==============================================================================
	// Energy calculation withing the clipped box and inside narrow band
	for (k = coord.k_min; k <= coord.k_max; k++)
	{
		for (i = coord.i_min; i <= coord.i_max; i++)
		{
			for (j = coord.j_min; j <= coord.j_max; j++)
			{
				// 2D to 1D representation for i, j
				size_t x = x_new(i, j, imageLength);
				if (NFunction(destination[k][x], distTrans[k][x], NDelta) == 1)
				{
					energy += (destination[k][x] - distTrans[k][x])   * (destination[k][x] - distTrans[k][x]);
					counter++; // Points within the narrow band
				}
				clipPts++; // Count points within the clipbox
			}
		}
	}
	//==============================================================================
	bandPts = counter;
	clipBoxPts = clipPts;
	return (energy / (2 * counter));
}
//==============================================================================
dataType energyFunctionClipBandArea(dataType ** destination, dataType ** distTrans, ClipBox coord, size_t imageLength, dataType ** fixedNBandPtr, dataType ** movingNBandPtr, dataType imageForeground)
{
	size_t k, i, j, counter = 0, clipPts = 0;
	dataType energy = 0.0;
	//==============================================================================
	// Energy calculation withing the clipped box and inside narrow band
	for (k = coord.k_min; k <= coord.k_max; k++)
	{
		for (i = coord.i_min; i <= coord.i_max; i++)
		{
			for (j = coord.j_min; j <= coord.j_max; j++)
			{
				// 2D to 1D representation for i, j
				size_t x = x_new(i, j, imageLength);
				if (NFunctionBinary(fixedNBandPtr[k][x], movingNBandPtr[k][x], imageForeground) == 1)
				{
					energy += (destination[k][x] - distTrans[k][x])   * (destination[k][x] - distTrans[k][x]);
					counter++; // Points within the narrow band
				}
				clipPts++; // Count points within the clipbox
			}
		}
	}
	//==============================================================================
	bandPts = counter;
	clipBoxPts = clipPts;
	return (energy / (2 * counter));
}
//==============================================================================
dataType finiteDifX(dataType ** distPtr, dataType h, size_t x, size_t k, size_t i, size_t imageLength)
{
	if (i == 0) // Apply Forward Difference
	{
		return (distPtr[k][x + 1] - distPtr[k][x]) / (h);
	}
	else if (i >= imageLength - 1) // Apply Backward Difference
	{
		return (distPtr[k][x] - distPtr[k][x - 1]) / (h);
	}
	else // Apply Central Difference
	{
		return (distPtr[k][x + 1] - distPtr[k][x - 1]) / (2 * h);
	}
}
//==============================================================================
dataType finiteDifY(dataType ** distPtr, dataType h, size_t k, size_t i, size_t j, size_t imageLength, size_t imageWidth)
{
	size_t x_n, x_p;
	if (j == 0) // Apply Forward Difference
	{
		// 2D to 1D representation for i, j
		x_n = x_new(i, j + 1, imageLength);
		x_p = x_new(i, j, imageLength);

		return (distPtr[k][x_n] - distPtr[k][x_p]) / (h);
	}
	else if (j >= imageWidth - 1) // Apply Backward Difference
	{
		// 2D to 1D representation for i, j
		x_n = x_new(i, j, imageLength);
		x_p = x_new(i, j - 1, imageLength);

		return (distPtr[k][x_n] - distPtr[k][x_p]) / (h);
	}
	else // Apply Central Difference
	{
		// 2D to 1D representation for i, j
		x_n = x_new(i, j + 1, imageLength);
		x_p = x_new(i, j - 1, imageLength);

		return (distPtr[k][x_n] - distPtr[k][x_p]) / (2 * h);
	}
}
//==============================================================================
dataType finiteDifZ(dataType ** distPtr, dataType h, size_t x, size_t k, size_t i, size_t imageLength, size_t imageHeight)
{
	if (k == 0) // Apply Forward Difference
	{
		return (distPtr[k + 1][x] - distPtr[k][x]) / (h);
	}
	else if (k >= imageHeight - 1) // Apply Backward Difference
	{
		return (distPtr[k][x] - distPtr[k - 1][x]) / (h);
	}
	else // Apply Central Difference
	{
		return (distPtr[k + 1][x] - distPtr[k - 1][x]) / (2 * h);
	}
}
//==============================================================================
void nbPointsX(dataType ** transformedBinaryData, dataType pixelSize, size_t x, size_t k, size_t i, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType * pValue, dataType * qValue, dataType *hh, ClipBox bestfitBox, Point3D * surface_points, size_t ptsNum, dataType insideShapevalue, bool parallelize)
{
	////==============================================================================
	size_t  dim2D = imageLength * imageWidth;
	//==============================================================================
	// neighbourPoints nb;
	dataType pv, qv;
	//==============================================================================
	if (i == 0) // Beginings
	{
		pv = getDistance(transformedBinaryData, imageHeight, imageLength, dim2D, k, x + 1, bestfitBox, surface_points, ptsNum, insideShapevalue, parallelize);
		qv = getDistance(transformedBinaryData, imageHeight, imageLength, dim2D, k, x, bestfitBox, surface_points, ptsNum, insideShapevalue, parallelize);
		*hh = pixelSize * 1;
	}
	else if (i >= imageLength - 1) // End
	{
		pv = getDistance(transformedBinaryData, imageHeight, imageLength, dim2D, k, x, bestfitBox, surface_points, ptsNum, insideShapevalue, parallelize);
		qv = getDistance(transformedBinaryData, imageHeight, imageLength, dim2D, k, x - 1, bestfitBox, surface_points, ptsNum, insideShapevalue, parallelize);
		*hh = pixelSize * 1;
	}
	else // Central values
	{
		pv = getDistance(transformedBinaryData, imageHeight, imageLength, dim2D, k, x + 1, bestfitBox, surface_points, ptsNum, insideShapevalue, parallelize);
		qv = getDistance(transformedBinaryData, imageHeight, imageLength, dim2D, k, x - 1, bestfitBox, surface_points, ptsNum, insideShapevalue, parallelize);
		//*qValue = distanceResultsPtr[k][x - 1];
		*hh = pixelSize * 2;
	}
	// *pValue = nb.pv; *qValue = nb.qv;
	*pValue = pv; *qValue = qv;
}
//==============================================================================
void nbPointsY(dataType ** transformedBinaryData, dataType pixelSize, size_t k, size_t i, size_t j, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType * pValue, dataType * qValue, dataType *hh, ClipBox bestfitBox, Point3D * surface_points, size_t ptsNum, dataType insideShapevalue, bool parallelize)
{
	////==============================================================================
	size_t  dim2D = imageLength * imageWidth;
	//==============================================================================
	// neighbourPoints nb;
	dataType pv, qv;
	//==============================================================================
	size_t x_n, x_p;
	if (j == 0) // Apply Forward Difference
	{
		// 2D to 1D representation for i, j
		x_n = x_new(i, j + 1, imageLength);
		x_p = x_new(i, j, imageLength);

		pv = getDistance(transformedBinaryData, imageHeight, imageLength, dim2D, k, x_n, bestfitBox, surface_points, ptsNum, insideShapevalue, parallelize);
		qv = getDistance(transformedBinaryData, imageHeight, imageLength, dim2D, k, x_p, bestfitBox, surface_points, ptsNum, insideShapevalue, parallelize);
		*hh = pixelSize * 1;
	}
	else if (j >= imageWidth - 1) // Apply Backward Difference
	{
		// 2D to 1D representation for i, j
		x_n = x_new(i, j, imageLength);
		x_p = x_new(i, j - 1, imageLength);

		pv = getDistance(transformedBinaryData, imageHeight, imageLength, dim2D, k, x_n, bestfitBox, surface_points, ptsNum, insideShapevalue, parallelize);
		qv = getDistance(transformedBinaryData, imageHeight, imageLength, dim2D, k, x_p, bestfitBox, surface_points, ptsNum, insideShapevalue, parallelize);
		*hh = pixelSize * 1;
	}
	else // Apply Central Difference
	{
		// 2D to 1D representation for i, j
		x_n = x_new(i, j + 1, imageLength);
		x_p = x_new(i, j - 1, imageLength);

		pv = getDistance(transformedBinaryData, imageHeight, imageLength, dim2D, k, x_n, bestfitBox, surface_points, ptsNum, insideShapevalue, parallelize);
		qv = getDistance(transformedBinaryData, imageHeight, imageLength, dim2D, k, x_p, bestfitBox, surface_points, ptsNum, insideShapevalue, parallelize);
		*hh = pixelSize * 2;
	}
	// *pValue = nb.pv; *qValue = nb.qv;
	*pValue = pv; *qValue = qv;
}
//==============================================================================
void nbPointsZ(dataType ** transformedBinaryData, dataType pixelSize, size_t x, size_t k, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType * pValue, dataType * qValue, dataType *hh, ClipBox bestfitBox, Point3D * surface_points, size_t ptsNum, dataType insideShapevalue, bool parallelize)
{
	////==============================================================================
	size_t  dim2D = imageLength * imageWidth;
	//==============================================================================
	// neighbourPoints nb;
	dataType pv, qv;
	//==============================================================================
	if (k == 0) // Apply Forward Difference
	{
		pv = getDistance(transformedBinaryData, imageHeight, imageLength, dim2D, k + 1, x, bestfitBox, surface_points, ptsNum, insideShapevalue, parallelize);
		qv = getDistance(transformedBinaryData, imageHeight, imageLength, dim2D, k, x, bestfitBox, surface_points, ptsNum, insideShapevalue, parallelize);
		*hh = pixelSize * 1;
		//==============================================================================
	}
	else if (k >= imageHeight - 1) // Apply Backward Difference
	{
		pv = getDistance(transformedBinaryData, imageHeight, imageLength, dim2D, k, x, bestfitBox, surface_points, ptsNum, insideShapevalue, parallelize);
		qv = getDistance(transformedBinaryData, imageHeight, imageLength, dim2D, k - 1, x, bestfitBox, surface_points, ptsNum, insideShapevalue, parallelize);
		//*qValue = distanceResultsPtr[k - 1][x];
		*hh = pixelSize * 1;
		//==============================================================================
	}
	else // Apply Central Difference
	{
		pv = getDistance(transformedBinaryData, imageHeight, imageLength, dim2D, k + 1, x, bestfitBox, surface_points, ptsNum, insideShapevalue, parallelize);
		qv = getDistance(transformedBinaryData, imageHeight, imageLength, dim2D, k - 1, x, bestfitBox, surface_points, ptsNum, insideShapevalue, parallelize);
		*hh = pixelSize * 2;
		//==============================================================================
	}
	// *pValue = nb.pv; *qValue = nb.qv;
	*pValue = pv; *qValue = qv;
}
//==============================================================================
dataType finiteDifAll(dataType pValue, dataType qValue, dataType hh)
{
	return (pValue - qValue) / hh;
}
//==============================================================================
Affine_Parameter gradientComponents(dataType ** destPtr, dataType ** distTrans, dataType h, Affine_Parameter * params, size_t imageHeight, size_t imageLength, size_t imageWidth)
{
	Affine_Parameter results;
	// Initialize the parameters
	size_t k, i, j, x, counter = 0;
	// Initialize the results
	results.rotation.x = 0.0, results.rotation.y = 0.0, results.rotation.z = 0.0;
	results.scaling.x = 0.0, results.scaling.y = 0.0, results.scaling.z = 0.0;
	results.translation.x = 0.0, results.translation.y = 0.0, results.translation.z = 0.0;
	// Forward difference parameters
	dataType xFwd, yFwd, zFwd;
	// Derivative component
	dataType componentX, componentY, componentZ;
	// Stores the difference between two distance pointers
	dataType distDifference;
	// Shorter Transformation names
	dataType phi = params->rotation.x, theta = params->rotation.y, psi = params->rotation.z;
	dataType sx = params->scaling.x, sy = params->scaling.y, sz = params->scaling.z;
	dataType tx = params->translation.x, ty = params->translation.y, tz = params->translation.z;
	// Trigonometric functions as const
	// Rotation
	const double neg_sin_phi_sin_psi = -sin(phi)*sin(psi), _cos_phi_cos_psi_sin_theta = cos(phi)*cos(psi)*sin(theta), neg_cos_psi_sin_phi = -cos(psi)*sin(phi), _cos_phi_sin_psi_sin_theta = cos(phi)*sin(psi)*sin(theta), _cos_phi_cos_theta = cos(phi)*cos(theta);
	const double _cos_phi_sin_psi = cos(phi)*sin(psi), neg_cos_phi_sin_psi = -cos(phi)*sin(psi), _cos_psi_sin_phi_sin_theta = cos(psi)*sin(phi)*sin(theta), _cos_phi_cos_psi = cos(phi)*cos(psi), _sin_phi_sin_psi_sin_theta = sin(phi)*sin(psi)*sin(theta), _cos_theta_sin_phi = cos(theta) * sin(phi);
	const double _cos_psi_sin_theta = cos(psi)*sin(theta), _sin_psi_sin_theta = sin(psi)*sin(theta), _cos_theta = cos(theta), _cos_psi_cos_theta_sin_phi = cos(psi)*cos(theta)*sin(phi), _cos_theta_sin_phi_sin_psi = cos(theta)*sin(phi)*sin(psi), _sin_phi_sin_theta = sin(phi)*sin(theta);
	const double _cos_phi_cos_psi_cos_theta = cos(phi)*cos(psi)*cos(theta), _cos_phi_cos_theta_sin_psi_ = cos(phi)*cos(theta)*sin(psi), _cos_phi_sin_theta = cos(phi)*sin(theta);
	// Scaling
	const double _cos_psi_cos_theta = cos(psi)*cos(theta), _cos_theta_sin_psi = cos(theta)*sin(psi), _sin_theta = sin(theta);
	const double _sin_phi_sin_psi = sin(phi)*sin(psi);
	const double _cos_psi_sin_phi = cos(psi)*sin(phi);
	// Scales
	const dataType _sx_sx = sx * sx, _sy_sy = sy * sy, _sz_sz = sz * sz;
	// Begin Evaluation
	for (k = 0; k < imageHeight; k++)
	{
		for (j = 0; j < imageWidth; j++)
		{
			x = x_new(0, j, imageLength);

			for (i = 0; i < imageLength; i++)
			{
				// 2D to 1D representation for i, j
				/*x = x_new(i, j, imageLength);*/
				if (NFunction(destPtr[k][x], distTrans[k][x], NDelta) == 1)
				{
					counter++;
					// Store the distance function difference
					distDifference = (dataType)((destPtr[k][x] - distTrans[k][x]) * 2.0);
					// Directional component vector derivatives - i, j, k
					dataType tmpI = i / (dataType)imageLength, tmpJ = j / (dataType)imageWidth, tmpK = k / (dataType)imageHeight;
					// Apply Forward Differences to the distTrans pointer
					xFwd = finiteDifX(distTrans, h, x, k, i, imageLength);
					yFwd = finiteDifY(distTrans, h, k, i, j, imageLength, imageWidth);
					zFwd = finiteDifZ(distTrans, h, x, k, i, imageLength, imageHeight);
					// Evaluate Individual Gradient Components
#ifdef DIRECTIONAL
					// Rotation Components - Directionnal
					componentX = (dataType)(yFwd * (((tmpI)*((neg_sin_phi_sin_psi + _cos_phi_cos_psi_sin_theta) / sy)) + ((tmpJ)*((neg_cos_psi_sin_phi - _cos_phi_sin_psi_sin_theta) / sy)) + ((-tmpK)*((_cos_phi_cos_theta) / sy))) +
						zFwd * (((tmpI)*((_cos_phi_sin_psi + _cos_psi_sin_phi_sin_theta) / sz)) + ((tmpJ)*((_cos_phi_cos_psi - _sin_phi_sin_psi_sin_theta) / sz)) + ((-tmpK)*((_cos_theta_sin_phi) / sz))));
					results.rotation.x += (componentX)*(distDifference);
					componentY = (dataType)(xFwd * (((-tmpI)*((_cos_psi_sin_theta) / sx)) + ((tmpJ)*((_sin_psi_sin_theta) / sx)) + ((tmpK)*((_cos_theta) / sx))) +
						yFwd * (((tmpI)*((_cos_psi_cos_theta_sin_phi) / sy)) + ((-tmpJ)*((_cos_theta_sin_phi_sin_psi) / sy)) + ((tmpK)*((_sin_phi_sin_theta) / sy))) +
						zFwd * (((-tmpI)*((_cos_phi_cos_psi_cos_theta) / sz)) + ((tmpJ)*((_cos_phi_cos_theta_sin_psi_) / sz)) + ((-tmpK)*(_cos_phi_sin_theta / sz))));
					results.rotation.y += (componentY)*(distDifference);
					componentZ = (dataType)(xFwd * (((-tmpI)*((_cos_theta_sin_psi) / sx)) + ((-tmpJ)*((_cos_psi_cos_theta) / sx))) +
						yFwd * (((tmpI)*((_cos_phi_cos_psi - _sin_phi_sin_psi_sin_theta) / sy)) + ((tmpJ)*((neg_cos_phi_sin_psi - _cos_psi_sin_phi_sin_theta) / sy))) +
						zFwd * (((tmpI)*((_cos_psi_sin_phi + _cos_phi_sin_psi_sin_theta) / sz)) + ((tmpJ)*((neg_sin_phi_sin_psi + _cos_phi_cos_psi_sin_theta) / sz))));					
					results.rotation.z += (componentZ)*(distDifference);
					// Directional Scale Components
					componentX = (dataType)(xFwd * ((-tmpI)*((_cos_psi_cos_theta) / (_sx_sx)) + ((tmpJ)*((_cos_theta_sin_psi) / (_sx_sx))) + ((-tmpK)*((_sin_theta) / (_sx_sx)))));
					results.scaling.x += (componentX)*(distDifference);
					componentY = (dataType)(yFwd * (((-tmpI)*((_cos_phi_sin_psi + _cos_psi_sin_phi_sin_theta) / (_sy_sy))) + ((-tmpJ)*((_cos_phi_cos_psi - _sin_phi_sin_psi_sin_theta) / (_sy_sy))) + ((tmpK)*((_cos_theta_sin_phi) / (_sy_sy)))));
					results.scaling.y += (componentY)*(distDifference);
					componentZ = (dataType)(zFwd * (((-tmpI)*((_sin_phi_sin_psi - _cos_phi_cos_psi_sin_theta) / (_sz_sz))) + ((-tmpJ)*((_cos_psi_sin_phi + _cos_phi_sin_psi_sin_theta) / (_sz_sz))) + ((-tmpK)*((_cos_phi_cos_theta) / (_sz_sz)))));
					results.scaling.z += (componentZ)*(distDifference);			
#endif // DIRECTIONAL
					// Translation Parameters - Always directional
					// Tx
					results.translation.x += (-xFwd)*(distDifference);
					// Ty
					results.translation.y += (-yFwd)*(distDifference);
					// Tz
					results.translation.z += (-zFwd)*(distDifference);
				}
				x++;
			}
		}
	}
	// Normalize results
	if (counter == 0)
	{
		counter = 1;
	}
	// NORMALIZATION
	results.scaling.x = results.scaling.x / counter;
	results.scaling.y = results.scaling.y / counter;
	results.scaling.z = results.scaling.z / counter;

	results.rotation.x = results.rotation.x / counter;
	results.rotation.y = results.rotation.y / counter;
	results.rotation.z = results.rotation.z / counter;

	results.translation.x = results.translation.x / counter;
	// Ty
	results.translation.y = results.translation.y / counter;
	// Tz
	results.translation.z = results.translation.z / counter;
	return results;
}
//==============================================================================
Affine_Parameter gradientCoorDinateDescentComp(dataType ** destPtr, dataType ** distTrans, dataType h, Affine_Parameter * params, size_t imageHeight, size_t imageLength, size_t imageWidth, size_t updateComponent)
{
	Affine_Parameter results;
	// Initialize the parameters
	size_t k, i, j, x, counter = 0;
	// Initialize the results
	results.rotation.x = 0.0, results.rotation.y = 0.0, results.rotation.z = 0.0;
	results.scaling.x = 0.0, results.scaling.y = 0.0, results.scaling.z = 0.0;
	results.translation.x = 0.0, results.translation.y = 0.0, results.translation.z = 0.0;
	// Forward difference parameters
	dataType xFwd, yFwd, zFwd;
	// Derivative component
	dataType componentX, componentY, componentZ;
	// Stores the difference between two distance pointers
	dataType distDifference;
	// Shorter Transformation names
	dataType phi = params->rotation.x, theta = params->rotation.y, psi = params->rotation.z;
	dataType sx = params->scaling.x, sy = params->scaling.y, sz = params->scaling.z;
	dataType tx = params->translation.x, ty = params->translation.y, tz = params->translation.z;
	// Trigonometric functions as const
	// Rotation
	const double neg_sin_phi_sin_psi = -sin(phi)*sin(psi), _cos_phi_cos_psi_sin_theta = cos(phi)*cos(psi)*sin(theta), neg_cos_psi_sin_phi = -cos(psi)*sin(phi), _cos_phi_sin_psi_sin_theta = cos(phi)*sin(psi)*sin(theta), _cos_phi_cos_theta = cos(phi)*cos(theta);
	const double _cos_phi_sin_psi = cos(phi)*sin(psi), neg_cos_phi_sin_psi = -cos(phi)*sin(psi), _cos_psi_sin_phi_sin_theta = cos(psi)*sin(phi)*sin(theta), _cos_phi_cos_psi = cos(phi)*cos(psi), _sin_phi_sin_psi_sin_theta = sin(phi)*sin(psi)*sin(theta), _cos_theta_sin_phi = cos(theta) * sin(phi);
	const double _cos_psi_sin_theta = cos(psi)*sin(theta), _sin_psi_sin_theta = sin(psi)*sin(theta), _cos_theta = cos(theta), _cos_psi_cos_theta_sin_phi = cos(psi)*cos(theta)*sin(phi), _cos_theta_sin_phi_sin_psi = cos(theta)*sin(phi)*sin(psi), _sin_phi_sin_theta = sin(phi)*sin(theta);
	const double _cos_phi_cos_psi_cos_theta = cos(phi)*cos(psi)*cos(theta), _cos_phi_cos_theta_sin_psi_ = cos(phi)*cos(theta)*sin(psi), _cos_phi_sin_theta = cos(phi)*sin(theta);
	// Scaling
	const double _cos_psi_cos_theta = cos(psi)*cos(theta), _cos_theta_sin_psi = cos(theta)*sin(psi), _sin_theta = sin(theta);
	const double _sin_phi_sin_psi = sin(phi)*sin(psi);
	const double _cos_psi_sin_phi = cos(psi)*sin(phi);
	// Scales
	const dataType _sx_sx = sx * sx, _sy_sy = sy * sy, _sz_sz = sz * sz;
	// Begin Evaluation
	for (k = 0; k < imageHeight; k++)
	{
		for (j = 0; j < imageWidth; j++)
		{
			x = x_new(0, j, imageLength);

			for (i = 0; i < imageLength; i++)
			{
				// 2D to 1D representation for i, j
				/*x = x_new(i, j, imageLength);*/
				if (NFunction(destPtr[k][x], distTrans[k][x], NDelta) == 1)
				{
					counter++;
					// Store the distance function difference
					distDifference = (dataType)((destPtr[k][x] - distTrans[k][x]) * 2.0);

					// Directional component vector derivatives - i, j, k
					dataType tmpI = i / (dataType)imageLength, tmpJ = j / (dataType)imageWidth, tmpK = k / (dataType)imageHeight;
					// Trignometry functions inside the component evaluation equation
					//T a = cos(phi), b = sin(phi), ab = cos(phi)*sin(phi), aa = cos(phi)*cos(phi), bb = sin(phi)*sin(phi), aaa = a * aa, bbb = b * bb, aab = aa * b, abb = a * bb;
					//==============================================================================
					// Apply Forward Differences to the distTrans pointer
					xFwd = finiteDifX(distTrans, h, x, k, i, imageLength);
					yFwd = finiteDifY(distTrans, h, k, i, j, imageLength, imageWidth);
					zFwd = finiteDifZ(distTrans, h, x, k, i, imageLength, imageHeight);
					//==============================================================================
#ifdef DIRECTIONAL
					// Begin update for Passed component
					//==============================================================================
					if (updateComponent == 1) // Rotation Component
					{
						componentX = (dataType)(yFwd * (((tmpI)*((neg_sin_phi_sin_psi + _cos_phi_cos_psi_sin_theta) / sy)) + ((tmpJ)*((neg_cos_psi_sin_phi - _cos_phi_sin_psi_sin_theta) / sy)) + ((-tmpK)*((_cos_phi_cos_theta) / sy))) +
							zFwd * (((tmpI)*((_cos_phi_sin_psi + _cos_psi_sin_phi_sin_theta) / sz)) + ((tmpJ)*((_cos_phi_cos_psi - _sin_phi_sin_psi_sin_theta) / sz)) + ((-tmpK)*((_cos_theta_sin_phi) / sz))));
						results.rotation.x += (componentX)*(distDifference);
						componentY = (dataType)(xFwd * (((-tmpI)*((_cos_psi_sin_theta) / sx)) + ((tmpJ)*((_sin_psi_sin_theta) / sx)) + ((tmpK)*((_cos_theta) / sx))) +
							yFwd * (((tmpI)*((_cos_psi_cos_theta_sin_phi) / sy)) + ((-tmpJ)*((_cos_theta_sin_phi_sin_psi) / sy)) + ((tmpK)*((_sin_phi_sin_theta) / sy))) +
							zFwd * (((-tmpI)*((_cos_phi_cos_psi_cos_theta) / sz)) + ((tmpJ)*((_cos_phi_cos_theta_sin_psi_) / sz)) + ((-tmpK)*(_cos_phi_sin_theta / sz))));
						results.rotation.y += (componentY)*(distDifference);
						componentZ = (dataType)(xFwd * (((-tmpI)*((_cos_theta_sin_psi) / sx)) + ((-tmpJ)*((_cos_psi_cos_theta) / sx))) +
							yFwd * (((tmpI)*((_cos_phi_cos_psi - _sin_phi_sin_psi_sin_theta) / sy)) + ((tmpJ)*((neg_cos_phi_sin_psi - _cos_psi_sin_phi_sin_theta) / sy))) +
							zFwd * (((tmpI)*((_cos_psi_sin_phi + _cos_phi_sin_psi_sin_theta) / sz)) + ((tmpJ)*((neg_sin_phi_sin_psi + _cos_phi_cos_psi_sin_theta) / sz))));
						results.rotation.z += (componentZ)*(distDifference);
					}
					else if (updateComponent == 2) // Scaling Component
					{
						componentX = (dataType)(xFwd * ((-tmpI)*((_cos_psi_cos_theta) / (_sx_sx)) + ((tmpJ)*((_cos_theta_sin_psi) / (_sx_sx))) + ((-tmpK)*((_sin_theta) / (_sx_sx)))));
						results.scaling.x += (componentX)*(distDifference);
						componentY = (dataType)(yFwd * (((-tmpI)*((_cos_phi_sin_psi + _cos_psi_sin_phi_sin_theta) / (_sy_sy))) + ((-tmpJ)*((_cos_phi_cos_psi - _sin_phi_sin_psi_sin_theta) / (_sy_sy))) + ((tmpK)*((_cos_theta_sin_phi) / (_sy_sy)))));
						results.scaling.y += (componentY)*(distDifference);
						componentZ = (dataType)(zFwd * (((-tmpI)*((_sin_phi_sin_psi - _cos_phi_cos_psi_sin_theta) / (_sz_sz))) + ((-tmpJ)*((_cos_psi_sin_phi + _cos_phi_sin_psi_sin_theta) / (_sz_sz))) + ((-tmpK)*((_cos_phi_cos_theta) / (_sz_sz)))));
						results.scaling.z += (componentZ)*(distDifference);
					}
					else if (updateComponent == 3) // Translation Component
					{
						// Tx
						results.translation.x += (-xFwd)*(distDifference);
						// Ty
						results.translation.y += (-yFwd)*(distDifference);
						// Tz
						results.translation.z += (-zFwd)*(distDifference);
					}
					//==============================================================================
#endif // DIRECTIONAL
				}
				x++;
			}
		}
	}
	// Normalize results
	if (counter == 0)
	{
		counter = 1;
	}
	// NORMALIZATION
	if (updateComponent == 1) // Rotation
	{
		results.rotation.x = results.rotation.x / counter;
		results.rotation.y = results.rotation.y / counter;
		results.rotation.z = results.rotation.z / counter;
	}
	else if (updateComponent == 2) // Scaling
	{
		results.scaling.x = results.scaling.x / counter;
		results.scaling.y = results.scaling.y / counter;
		results.scaling.z = results.scaling.z / counter;
	}
	else if (updateComponent == 3) // Translation
	{
		// Tx
		results.translation.x = results.translation.x / counter;
		// Ty
		results.translation.y = results.translation.y / counter;
		// Tz
		results.translation.z = results.translation.z / counter;
	}
	return results;
}
//==============================================================================
Affine_Parameter gradientComponentsClip(dataType ** destPtr, dataType ** distTrans, dataType h, Affine_Parameter * params, size_t imageHeight, size_t imageLength, size_t imageWidth, ClipBox bestFit)
{
	Affine_Parameter results;
	// Initialize the parameters
	size_t k, i, j, x, counter = 0;

	// Initialize the results
	results.rotation.x = 0.0, results.rotation.y = 0.0, results.rotation.z = 0.0;
	results.scaling.x = 0.0, results.scaling.y = 0.0, results.scaling.z = 0.0;
	results.translation.x = 0.0, results.translation.y = 0.0, results.translation.z = 0.0;

	// Forward difference parameters
	dataType xFwd, yFwd, zFwd;

	// Derivative component
	dataType componentX, componentY, componentZ;

	// Stores the difference between two distance pointers
	dataType distDifference;

	// Shorter Transformation names
	dataType phi = params->rotation.x, theta = params->rotation.y, psi = params->rotation.z;
	dataType sx = params->scaling.x, sy = params->scaling.y, sz = params->scaling.z;
	dataType tx = params->translation.x, ty = params->translation.y, tz = params->translation.z;
	// Trigonometric functions as const
	// Rotation
	const double neg_sin_phi_sin_psi = -sin(phi)*sin(psi), _cos_phi_cos_psi_sin_theta = cos(phi)*cos(psi)*sin(theta), neg_cos_psi_sin_phi = -cos(psi)*sin(phi), _cos_phi_sin_psi_sin_theta = cos(phi)*sin(psi)*sin(theta), _cos_phi_cos_theta = cos(phi)*cos(theta);
	const double _cos_phi_sin_psi = cos(phi)*sin(psi), neg_cos_phi_sin_psi = -cos(phi)*sin(psi), _cos_psi_sin_phi_sin_theta = cos(psi)*sin(phi)*sin(theta), _cos_phi_cos_psi = cos(phi)*cos(psi), _sin_phi_sin_psi_sin_theta = sin(phi)*sin(psi)*sin(theta), _cos_theta_sin_phi = cos(theta) * sin(phi);
	const double _cos_psi_sin_theta = cos(psi)*sin(theta), _sin_psi_sin_theta = sin(psi)*sin(theta), _cos_theta = cos(theta), _cos_psi_cos_theta_sin_phi = cos(psi)*cos(theta)*sin(phi), _cos_theta_sin_phi_sin_psi = cos(theta)*sin(phi)*sin(psi), _sin_phi_sin_theta = sin(phi)*sin(theta);
	const double _cos_phi_cos_psi_cos_theta = cos(phi)*cos(psi)*cos(theta), _cos_phi_cos_theta_sin_psi_ = cos(phi)*cos(theta)*sin(psi), _cos_phi_sin_theta = cos(phi)*sin(theta);
	// Scaling
	const double _cos_psi_cos_theta = cos(psi)*cos(theta), _cos_theta_sin_psi = cos(theta)*sin(psi), _sin_theta = sin(theta);
	const double _sin_phi_sin_psi = sin(phi)*sin(psi);
	const double _cos_psi_sin_phi = cos(psi)*sin(phi);
	// Scales
	const dataType _sx_sx = sx * sx, _sy_sy = sy * sy, _sz_sz = sz * sz;

	// Begin Evaluation
	for (k = bestFit.k_min; k < bestFit.k_max + 1; k++)
	{
		for (j = bestFit.j_min; j < bestFit.j_max + 1; j++)
		{
			// 2D to 1D representation for i, j
			x = x_new(0, j, imageLength);
			for (i = bestFit.i_min; i < bestFit.i_max + 1; i++)
			{
				// 2D to 1D representation for i, j
				//x = x_new(i, j, imageLength);
				if (NFunction(destPtr[k][x], distTrans[k][x], NDelta) == 1)
				{
					counter++;
					// Store the distance function difference
					distDifference = (dataType)((destPtr[k][x] - distTrans[k][x]) * 2.0);
					// Directional component vector derivatives - i, j, k
					dataType tmpI = i / (dataType)imageLength, tmpJ = j / (dataType)imageWidth, tmpK = k / (dataType)imageHeight;
					// Apply Forward Differences to the distTrans pointer
					xFwd = finiteDifX(distTrans, h, x, k, i, imageLength);
					yFwd = finiteDifY(distTrans, h, k, i, j, imageLength, imageWidth);
					zFwd = finiteDifZ(distTrans, h, x, k, i, imageLength, imageHeight);
					// Evaluate Individual Gradient Components
#ifdef DIRECTIONAL
					// Rotation Components - Directionnal
					componentX = (dataType)(yFwd * (((tmpI)*((neg_sin_phi_sin_psi + _cos_phi_cos_psi_sin_theta) / sy)) + ((tmpJ)*((neg_cos_psi_sin_phi - _cos_phi_sin_psi_sin_theta) / sy)) + ((-tmpK)*((_cos_phi_cos_theta) / sy))) +
						zFwd * (((tmpI)*((_cos_phi_sin_psi + _cos_psi_sin_phi_sin_theta) / sz)) + ((tmpJ)*((_cos_phi_cos_psi - _sin_phi_sin_psi_sin_theta) / sz)) + ((-tmpK)*((_cos_theta_sin_phi) / sz))));
					results.rotation.x += (componentX)*(distDifference);
					componentY = (dataType)(xFwd * (((-tmpI)*((_cos_psi_sin_theta) / sx)) + ((tmpJ)*((_sin_psi_sin_theta) / sx)) + ((tmpK)*((_cos_theta) / sx))) +
						yFwd * (((tmpI)*((_cos_psi_cos_theta_sin_phi) / sy)) + ((-tmpJ)*((_cos_theta_sin_phi_sin_psi) / sy)) + ((tmpK)*((_sin_phi_sin_theta) / sy))) +
						zFwd * (((-tmpI)*((_cos_phi_cos_psi_cos_theta) / sz)) + ((tmpJ)*((_cos_phi_cos_theta_sin_psi_) / sz)) + ((-tmpK)*(_cos_phi_sin_theta / sz))));
					results.rotation.y += (componentY)*(distDifference);
					componentZ = (dataType)(xFwd * (((-tmpI)*((_cos_theta_sin_psi) / sx)) + ((-tmpJ)*((_cos_psi_cos_theta) / sx))) +
						yFwd * (((tmpI)*((_cos_phi_cos_psi - _sin_phi_sin_psi_sin_theta) / sy)) + ((tmpJ)*((neg_cos_phi_sin_psi - _cos_psi_sin_phi_sin_theta) / sy))) +
						zFwd * (((tmpI)*((_cos_psi_sin_phi + _cos_phi_sin_psi_sin_theta) / sz)) + ((tmpJ)*((neg_sin_phi_sin_psi + _cos_phi_cos_psi_sin_theta) / sz))));
					results.rotation.z += (componentZ)*(distDifference);
					// Directional Scale Components
					componentX = (dataType)(xFwd * ((-tmpI)*((_cos_psi_cos_theta) / (_sx_sx)) + ((tmpJ)*((_cos_theta_sin_psi) / (_sx_sx))) + ((-tmpK)*((_sin_theta) / (_sx_sx)))));
					results.scaling.x += (componentX)*(distDifference);
					componentY = (dataType)(yFwd * (((-tmpI)*((_cos_phi_sin_psi + _cos_psi_sin_phi_sin_theta) / (_sy_sy))) + ((-tmpJ)*((_cos_phi_cos_psi - _sin_phi_sin_psi_sin_theta) / (_sy_sy))) + ((tmpK)*((_cos_theta_sin_phi) / (_sy_sy)))));
					results.scaling.y += (componentY)*(distDifference);
					componentZ = (dataType)(zFwd * (((-tmpI)*((_sin_phi_sin_psi - _cos_phi_cos_psi_sin_theta) / (_sz_sz))) + ((-tmpJ)*((_cos_psi_sin_phi + _cos_phi_sin_psi_sin_theta) / (_sz_sz))) + ((-tmpK)*((_cos_phi_cos_theta) / (_sz_sz)))));
					results.scaling.z += (componentZ)*(distDifference);
#endif // DIRECTIONAL
					// Translation Parameters - Always directional
					// Tx
					results.translation.x += (-xFwd)*(distDifference);
					// Ty
					results.translation.y += (-yFwd)*(distDifference);
					// Tz
					results.translation.z += (-zFwd)*(distDifference);
				}
				x++;
			}
		}
	}
	// Normalize results
	if (counter == 0)
	{
		counter = 1;
	}
	// NORMALIZATION
	results.scaling.x = results.scaling.x / counter;
	results.scaling.y = results.scaling.y / counter;
	results.scaling.z = results.scaling.z / counter;

	results.rotation.x = results.rotation.x / counter;
	results.rotation.y = results.rotation.y / counter;
	results.rotation.z = results.rotation.z / counter;

	results.translation.x = results.translation.x / counter;
	// Ty
	results.translation.y = results.translation.y / counter;
	// Tz
	results.translation.z = results.translation.z / counter;
	return results;
}
//==============================================================================
Affine_Parameter registration3D(dataType ** fixedData, dataType ** movingData, Affine_Parameter initTransform, dataType step_size, dataType tol, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType centroid[3], Registration_Params params)
{
	//==============================================================================
	size_t k, i, dim2D = imageLength * imageWidth;
	int iteration = 0;
	dataType firstCpuTime, secondCpuTime, regStartCpuTime, regStopCpuTime, regTotalCpuTimen = 0.;
	dataType energyTotalCpuTime = 0., distanceTotalCpuTime = 0., gradientTotalCpuTime = 0., transformationTotalCpuTime = 0.;
	//==============================================================================
	// Affine Parameters
	Affine_Parameter affineResult, affineTmp;
	//==============================================================================
	// Create a new shape Pointers to be used
	dataType ** destPtr = (dataType **)malloc(sizeof(dataType *) * imageHeight); // distances for destination
	//==============================================================================
	// ClipBox Variable
	dataType ** movInitPtr = (dataType **)malloc(sizeof(dataType *) * imageHeight); // distances for Moving
	//==============================================================================
	ClipBox coordFixed, coordMoving, bestFit, coordMovingTmp;
	//==============================================================================
	const size_t mem_alloc_2D_block = sizeof(dataType) * dim2D;
	//==============================================================================
	const dataType large_value = 50000;
	//==============================================================================
	if (params.use_clipbox)
	{
		for (i = 0; i < imageHeight; i++)
		{
			movInitPtr[i] = (dataType *)malloc(mem_alloc_2D_block);
		}
	}
	//==============================================================================
	// Initializations of Pointers
	for (i = 0; i < imageHeight; i++)
	{
		destPtr[i] = (dataType *)malloc(mem_alloc_2D_block);
	}
	//==============================================================================
	// Instantiate Affine Parameters
	affineResult.rotation = initTransform.rotation;
	affineResult.scaling = initTransform.scaling;
	affineResult.translation = initTransform.translation;
	//==============================================================================
	// Energy tmp optimal, stop boolean
	dataType energyTmp;
	bool stopCond = false;
	//==============================================================================
	// Apply distance function between transPtr and distTrans
	// Begin Record Time
#ifdef MEASURE_TIME
	firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
	//==============================================================================
	fastSweepingFunction_3D(destPtr, fixedData, imageLength, imageWidth, imageHeight, params.h, large_value, params.imageForeground);
	//==============================================================================
	if (params.use_clipbox)
	{
		// Initial dist. fn for moving image before adding any transformation
		fastSweepingFunction_3D(movInitPtr, movingData, imageLength, imageWidth, imageHeight, params.h, large_value, params.imageForeground);
	}
	//==============================================================================
#ifdef MEASURE_TIME
	secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
	// Store the time
	distanceTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
	//==============================================================================
	if (params.use_clipbox)
	{
		//==============================================================================
		// ClipBoxes from the calculated distances
		// Finding the clip box points for the fixed image
		coordFixed = findClipBoxSingle(destPtr, imageHeight, imageLength, imageWidth);
		//==============================================================================
		coordMoving = findClipBoxSingle(movInitPtr, imageHeight, imageLength, imageWidth);
		//==============================================================================
		// Free after
		for (k = 0; k < imageHeight; k++)
		{
			free(movInitPtr[k]);
		}
		free(movInitPtr);
		movInitPtr = NULL;
		//==============================================================================
	}
	//==============================================================================
#ifdef CONSOLE_OUTPUT
	printf("Distance calc before Registration CPU time: %e secs\n\n", secondCpuTime - firstCpuTime);
#endif
	//==============================================================================
	// Begin Registration of Distances between shapes
	// Start Timing the Registration Process
#ifdef MEASURE_TIME
	regStartCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
	//==============================================================================
	// Create a new shape Pointers to be used
	dataType ** transPtr = (dataType **)malloc(sizeof(dataType*) * imageHeight); // Transformed Ptr
	dataType ** distTransPtr = (dataType **)malloc(sizeof(dataType*) * imageHeight); // distances for Transformed Ptr
	for (i = 0; i < imageHeight; i++)
	{
		transPtr[i] = (dataType *)malloc(mem_alloc_2D_block);
		distTransPtr[i] = (dataType *)malloc(mem_alloc_2D_block);
	}
	//==============================================================================
	while (!stopCond)
	{
		//==============================================================================
		// Timing The Transformation Inside the registration function
#ifdef MEASURE_TIME
		firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
		transform3DImage(movingData, transPtr, affineResult.translation, affineResult.scaling, affineResult.rotation, imageHeight, imageLength, imageWidth, params.imageBackground, centroid, params.parallelize);
		//==============================================================================
#ifdef MEASURE_TIME
		secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
		// Store the time
		transformationTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
		//==============================================================================
		if (params.use_clipbox)
		{
			//==============================================================================
			// Transform the coordMoving clip box using calc. transform component results
			// Copy to coordMovingTmp
			coordMovingTmp = coordMoving;
			transformClip(&coordMovingTmp, affineResult.translation, affineResult.scaling, affineResult.rotation, centroid, imageHeight, imageLength, imageWidth);
			// Find the bestFit from transformed clip
			bestFit.k_min = min(coordFixed.k_min, coordMovingTmp.k_min);
			bestFit.i_min = min(coordFixed.i_min, coordMovingTmp.i_min);
			bestFit.j_min = min(coordFixed.j_min, coordMovingTmp.j_min);

			bestFit.k_max = max(coordFixed.k_max, coordMovingTmp.k_max);
			bestFit.i_max = max(coordFixed.i_max, coordMovingTmp.i_max);
			bestFit.j_max = max(coordFixed.j_max, coordMovingTmp.j_max);
			//==============================================================================
		}
		//==============================================================================
#ifdef CONSOLE_OUTPUT
		printf("Registration Transformation calc. CPU time at iteration %4d: %e secs\n", iteration, secondCpuTime - firstCpuTime);
#endif
		//==============================================================================
		// Apply distance function between transPtr and distTrans
		// Begin Record Time
#ifdef MEASURE_TIME
		firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
		//==============================================================================
		if (params.use_clipbox)
		{
			fSweeping3D(distTransPtr, transPtr, imageLength, imageWidth, imageHeight, params.h, large_value, params.imageForeground, bestFit);
		}
		else
		{
			fastSweepingFunction_3D(distTransPtr, transPtr, imageLength, imageWidth, imageHeight, params.h, large_value, params.imageForeground);
		}
		//==============================================================================
#ifdef MEASURE_TIME
		secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
		// Store the time
		distanceTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
		//==============================================================================
#ifdef CONSOLE_OUTPUT
		printf("Distance calc during Registration CPU time at iteration %4d: %e secs\n", iteration, secondCpuTime - firstCpuTime);
#endif
		//==============================================================================
		// Evaluate Energy Function
		// Begin Record Time
#ifdef MEASURE_TIME
		firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
		//==============================================================================
		if (params.use_clipbox)
		{
			// Evaluate Energy Function - L2 Norm Between the two calc. distances within the band and clipbox
			energyTmp = energyFunctionClip(destPtr, distTransPtr, bestFit, imageLength);
		}
		else
		{
			// Evaluate Energy Function - L2 Norm Between the two calc. distances
			energyTmp = energyFunction(destPtr, distTransPtr, imageHeight, imageLength, imageWidth, params.h);
		}
		//==============================================================================
#ifdef MEASURE_TIME
		secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
		// Store the time
		energyTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
		//==============================================================================
#ifdef CONSOLE_OUTPUT
		printf("Energy Function calc. CPU time at iteration %4d: %e secs\n\n", iteration, secondCpuTime - firstCpuTime);
#endif
		//==============================================================================
		// Print Pre-evaluate affine values
#ifdef CONSOLE_OUTPUT
		printf("Energy = %5.5lf, iteration %4d, Phi = %3.5lf, Theta = %3.5lf, Psi = %3.5lf, Sx = %2.5lf, Sy = %2.5lf, Sz = %2.5lf, Tx = %2.5lf, Ty = %2.5lf, Tz = %2.5lf\n",
			energyTmp, iteration, affineResult.rotation.x, affineResult.rotation.y, affineResult.rotation.z, affineResult.scaling.x, affineResult.scaling.y,
			affineResult.scaling.z, affineResult.translation.x, affineResult.translation.y, affineResult.translation.z);
#endif
		//==============================================================================
		// Check Stoping condition with tolerance and number of ierations
		if (energyTmp < tol || iteration == params.max_iterations)
		{
			stopCond = true;
			for (k = 0; k < imageHeight; k++)
			{
				free(destPtr[k]);
				free(transPtr[k]);
				free(distTransPtr[k]);
			}
			free(destPtr);
			destPtr = NULL;
			//
			free(transPtr);
			transPtr = NULL;
			//
			free(distTransPtr);
			distTransPtr = NULL;
			//==============================================================================
			printf("Total distance Function calc. CPU Time is: %e secs\n", distanceTotalCpuTime);
			printf("Total energy Function calc. CPU Time is: %e secs\n", energyTotalCpuTime);
			printf("Total gradient Function calc. CPU Time is: %e secs\n", gradientTotalCpuTime);
			printf("Total transformation Function calc. CPU Time is: %e secs\n", transformationTotalCpuTime);
			//==============================================================================
			// Print the Calculated Transformation Parameters At the End of Registration
			printf("Energy = %5.5lf, iteration %4d, Phi = %3.5lf, Theta = %3.5lf, Psi = %3.5lf, Sx = %2.5lf, Sy = %2.5lf, Sz = %2.5lf, Tx = %2.5lf, Ty = %2.5lf, Tz = %2.5lf\n",
				energyTmp, iteration, affineResult.rotation.x, affineResult.rotation.y, affineResult.rotation.z, affineResult.scaling.x, affineResult.scaling.y,
				affineResult.scaling.z, affineResult.translation.x, affineResult.translation.y, affineResult.translation.z);
			//==============================================================================
		}
		else
		{
			//==============================================================================
			// Begin Record Time
#ifdef MEASURE_TIME
			firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
			//==============================================================================
			if (params.use_clipbox)
			{
				affineTmp = gradientComponentsClip(destPtr, distTransPtr, params.h, &affineResult, imageHeight, imageLength, imageWidth, bestFit);
			}
			else
			{
				affineTmp = gradientComponents(destPtr, distTransPtr, params.h, &affineResult, imageHeight, imageLength, imageWidth);
			}
			//==============================================================================
#ifdef MEASURE_TIME
			secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
			// Store the time
			gradientTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
			//==============================================================================
#ifdef CONSOLE_OUTPUT
			printf("Gradient Components calc. CPU time at iteration %4d: %e secs\n\n", iteration, secondCpuTime - firstCpuTime);
#endif
			//==============================================================================
			// Set new values for affine temp results
			// Rotation
			affineResult.rotation.x += params.rotation_weight * step_size*affineTmp.rotation.x;
			affineResult.rotation.y += params.rotation_weight * step_size*affineTmp.rotation.y;
			affineResult.rotation.z += params.rotation_weight * step_size*affineTmp.rotation.z;
			// Scaling
			affineResult.scaling.x += params.scaling_weight * step_size*affineTmp.scaling.x;
			affineResult.scaling.y += params.scaling_weight * step_size*affineTmp.scaling.y;
			affineResult.scaling.z += params.scaling_weight * step_size*affineTmp.scaling.z;
			//Translation
			affineResult.translation.x += params.translation_weight * step_size*affineTmp.translation.x;
			affineResult.translation.y += params.translation_weight * step_size*affineTmp.translation.y;
			affineResult.translation.z += params.translation_weight * step_size*affineTmp.translation.z;
			//==============================================================================
			// Increase Iteration
			iteration++;
			//==============================================================================
		}
	}
	//==============================================================================
	// Stop Timing the Registration Process
#ifdef MEASURE_TIME
	regStopCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
	// Store Time For Each Registration run
	regTotalCpuTimen = regStopCpuTime - regStartCpuTime;
#endif
	//==============================================================================
#ifdef CONSOLE_OUTPUT
	printf("Total Registration Function calc. CPU Time is: %e secs\n\n", regTotalCpuTimen);
#endif
	//==============================================================================
	return affineResult;
	//==============================================================================
}
//==============================================================================
Affine_Parameter registrationCoorDinateDescent3D(dataType ** fixedData, dataType ** movingData, Affine_Parameter initTransform, dataType step_size, dataType tol, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType centroid[3], Registration_Params params)
{
	dataType rotation_weight = 1.0;
	// Scaling
	dataType scaling_weight = 1.0;
	// Translation
	dataType translation_weight = 1.0;
	//==============================================================================
	int components = 1, switchcomponent;
	//==============================================================================
	dataType stepsize = step_size;
	//==============================================================================
	size_t k, i, dim2D = imageLength * imageWidth;
	int iteration = 0, max_ter = 1000;
	int count_rejected = 0, state_accept, count_steps_reset = 0, max_resets = 9;
	//==============================================================================
	dataType firstCpuTime, secondCpuTime, regStartCpuTime, regStopCpuTime, regTotalCpuTimen = 0.;
	dataType energyTotalCpuTime = 0., distanceTotalCpuTime = 0., gradientTotalCpuTime = 0., transformationTotalCpuTime = 0.;
	// Affine tmp prev init
	Point3D rotationTran = { 0.0, 0.0, 0.0 };
	Point3D scalingTran = { 1.0, 1.0, 1.0 };
	Point3D translationTran = { 0.0, 0.0, 0.0 };
	Affine_Parameter affineResult, affineTmp, affineTmp_prev, affinePrev;
	affineTmp_prev.rotation = rotationTran, affineTmp_prev.scaling = scalingTran, affineTmp_prev.translation = translationTran;
	// Create a new shape Pointers to be used
	dataType ** destPtr = (dataType **)malloc(sizeof(dataType*) * imageHeight); // distances for destination
	//==============================================================================
	const size_t mem_alloc_2D_block = sizeof(dataType)*dim2D;
	//==============================================================================
	const dataType large_value = 50000;
	//==============================================================================
	// USE_CLIP
	ClipBox coordFixed, coordMoving, bestFit, coordMovingTmp; // Clipbox for bestFit of both fixed and moving images, Moving image clipbox
	dataType ** movInitPtr = (dataType **)malloc(sizeof(dataType *) * imageHeight); // distances for Moving
	if (params.use_clipbox)
	{
		for (i = 0; i < imageHeight; i++)
		{
			movInitPtr[i] = (dataType*)malloc(mem_alloc_2D_block);
		}
	}
	//==============================================================================
	// Initializations of Pointers
	for (i = 0; i < imageHeight; i++)
	{
		destPtr[i] = (dataType*)malloc(mem_alloc_2D_block);
	}
	// Initialize to same background - default value is 255, 0, 0
	//initialize3dArrayD(destPtr, imageLength, imageWidth, imageHeight, foregound);
	//==============================================================================
	// Instantiate Affine Parameters
	affineResult.rotation = initTransform.rotation;
	affineResult.scaling = initTransform.scaling;
	affineResult.translation = initTransform.translation;
	//==============================================================================
	Point3D init_trans = { 0,0,0 };
	Point3D init_rot = { 0,0,0 };
	Point3D init_scale = { 1,1,1 };
	affinePrev = affineResult; // Start sames as affine results
	//==============================================================================
	// Energy tmp optimal, stop boolean
	double energyTmp, prev_energy = DBL_MAX;
	bool stopCond = false;
	// Apply distance function between transPtr and distTrans
	// Begin Record Time
#ifdef MEASURE_TIME
	firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
	fastSweepingFunction_3D(destPtr, fixedData, imageLength, imageWidth, imageHeight, params.h, large_value, params.imageForeground);
	//==============================================================================
	// USE_CLIP
	if (params.use_clipbox)
	{
		// Initial dist. fn for moving image before adding any transformation
		fastSweepingFunction_3D(movInitPtr, movingData, imageLength, imageWidth, imageHeight, params.h, large_value, params.imageForeground);
	}
	//==============================================================================
#ifdef MEASURE_TIME
	secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
	// Store the time
	distanceTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
	//==============================================================================
	// USE_CLIP
	if (params.use_clipbox)
	{
		//==============================================================================
		// Finding the clip box points for the fixed image
		coordFixed = findClipBoxSingle(destPtr, imageHeight, imageLength, imageWidth);
		//==============================================================================
		coordMoving = findClipBoxSingle(movInitPtr, imageHeight, imageLength, imageWidth);
		//==============================================================================
		// Free after
		for (k = 0; k < imageHeight; k++)
		{
			free(movInitPtr[k]);
		}
		free(movInitPtr);
		movInitPtr = NULL;
		//==============================================================================
	}
	//==============================================================================
#ifdef CONSOLE_OUTPUT
	printf("Distance calc before Registration CPU time: %e secs\n\n", secondCpuTime - firstCpuTime);
#endif
	// Begin Registration of Distances between shapes
	// Start Timing the Registration Process
#ifdef MEASURE_TIME
	regStartCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
	//==============================================================================
	// Create a new shape Pointers to be used
	dataType ** transPtr = (dataType**)malloc(sizeof(dataType*) * imageHeight); // Transformed Ptr
	dataType ** distTransPtr = (dataType**)malloc(sizeof(dataType*) * imageHeight); // distances for Transformed Ptr
	for (i = 0; i < imageHeight; i++)
	{
		transPtr[i] = (dataType*)malloc(mem_alloc_2D_block);
		distTransPtr[i] = (dataType*)malloc(mem_alloc_2D_block);
	}
	//==============================================================================
	while (!stopCond)
	{
		// Timing The Transformation Inside the registration function
#ifdef MEASURE_TIME
		firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
		transform3DImage(movingData, transPtr, affineResult.translation, affineResult.scaling, affineResult.rotation, imageHeight, imageLength, imageWidth, params.imageBackground, centroid, params.parallelize);
#ifdef MEASURE_TIME
		secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
		// Store the time
		transformationTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
		//==============================================================================
	// USE_CLIP
		if (params.use_clipbox)
		{
			//==============================================================================
			// Transform the coordMoving clip box using calc. transform component results
			// Copy to coordMovingTmp
			coordMovingTmp = coordMoving;
			transformClip(&coordMovingTmp, affineResult.translation, affineResult.scaling, affineResult.rotation, centroid, imageHeight, imageLength, imageWidth);
			// Find the bestFit from transformed clip
			bestFit.k_min = min(coordFixed.k_min, coordMovingTmp.k_min);
			bestFit.i_min = min(coordFixed.i_min, coordMovingTmp.i_min);
			bestFit.j_min = min(coordFixed.j_min, coordMovingTmp.j_min);

			bestFit.k_max = max(coordFixed.k_max, coordMovingTmp.k_max);
			bestFit.i_max = max(coordFixed.i_max, coordMovingTmp.i_max);
			bestFit.j_max = max(coordFixed.j_max, coordMovingTmp.j_max);
			//==============================================================================
		}
		//==============================================================================
#ifdef CONSOLE_OUTPUT
		printf("Registration Transformation calc. CPU time at iteration %4d: %e secs\n", iteration, secondCpuTime - firstCpuTime);
#endif
		// Apply distance function between transPtr and distTrans
		// Begin Record Time
#ifdef MEASURE_TIME
		firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
		//==============================================================================
// USE_CLIP
		if (params.use_clipbox)
		{
			// fSweeping3D(distTransPtr, transPtr, imageLength, imageWidth, imageHeight, 1, large_value, foregound, bestFit);
			fastSweepingFunction_3D(distTransPtr, transPtr, imageLength, imageWidth, imageHeight, params.h, large_value, params.imageForeground);
		}
		else
		{
			fastSweepingFunction_3D(distTransPtr, transPtr, imageLength, imageWidth, imageHeight, params.h, large_value, params.imageForeground);
		}
		//==============================================================================
#ifdef MEASURE_TIME
		secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
		// Store the time
		distanceTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
#ifdef CONSOLE_OUTPUT
		printf("Distance calc during Registration CPU time at iteration %4d: %e secs\n", iteration, secondCpuTime - firstCpuTime);
#endif
		// Evaluate Energy Function - L2 Norm Between the two calc. distances
		// Begin Record Time
#ifdef MEASURE_TIME
		firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
		//==============================================================================
		// Evaluate Energy Function - L2 Norm Between the two calc. distances within the band and clipbox
// USE_CLIP
		if (params.use_clipbox)
		{
			energyTmp = energyFunctionClip(destPtr, distTransPtr, bestFit, imageLength);
			//energyTmp = energyFunction(destPtr, distTransPtr, imageHeight, imageLength, imageWidth, h);
		}
		else {
			energyTmp = energyFunction(destPtr, distTransPtr, imageHeight, imageLength, imageWidth, params.h);
		}
		//==============================================================================
#ifdef MEASURE_TIME
		secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
		// Store the time
		energyTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
		//==============================================================================
		dataType diff_abs = (dataType)fabs(energyTmp - prev_energy), accept_diff = (dataType) 1e-04;
		// Check if current energy has reduced from previous
		if (energyTmp <= prev_energy)
		{
			prev_energy = energyTmp;
			//==============================================================================
			state_accept = 0;
			//==============================================================================
			affinePrev = affineResult;
			//==============================================================================
		}
		else // Use a new component to calculate the energy if current gives worse
		{
			state_accept = 1;
			//==============================================================================
			// Use previous components that were accepted and recalculate the previous distance pointer.
			affineResult = affinePrev;
			//==============================================================================
		}
		//==============================================================================
		// Adjust step size after a full cycle through components
		if (count_rejected >= 3) // Maximum rejections
		{
			// Adjust step size after a full cycle through components
			stepsize = (dataType)(stepsize / 2.0);
			//==============================================================================
			count_rejected = 0; // Reset
			//==============================================================================
		}
		//==============================================================================
		// Check if gone lower than acceptable minimum
		if (stepsize < 0.004)
		{
			stepsize = step_size; // Reset to orignal to start all over again
			count_steps_reset++;
		}
		//==============================================================================
		if ((state_accept == 1))
		{
			count_rejected++; // Increment
			//==============================================================================
		}
		else if (state_accept == 0)
		{
			count_rejected = 0; // Reset
			count_steps_reset = 0;
			//==============================================================================
		}
		//==============================================================================
#ifdef CONSOLE_OUTPUT
		printf("Energy Function calc. CPU time at iteration %4d: %e secs\n\n", iteration, secondCpuTime - firstCpuTime);
#endif
		// Print Pre-evaluate affine values
#ifdef CONSOLE_OUTPUT
		printf("Energy = %5.5lf, iteration %4d, Phi = %3.5lf, Theta = %3.5lf, Psi = %3.5lf, Sx = %2.5lf, Sy = %2.5lf, Sz = %2.5lf, Tx = %2.5lf, Ty = %2.5lf, Tz = %2.5lf\n",
			energyTmp, iteration, affineResult.rotation.x, affineResult.rotation.y, affineResult.rotation.z, affineResult.scaling.x, affineResult.scaling.y,
			affineResult.scaling.z, affineResult.translation.x, affineResult.translation.y, affineResult.translation.z);
#endif
		//==============================================================================
		// Check Stoping condition with tolerance and number of ierations
		if (energyTmp < tol || iteration == max_ter || count_steps_reset == max_resets)
		{
			//==============================================================================
			affineResult = affinePrev;
			stopCond = true;
			//==============================================================================
			printf("Total distance Function calc. CPU Time is: %e secs\n", distanceTotalCpuTime);
			printf("Total energy Function calc. CPU Time is: %e secs\n", energyTotalCpuTime);
			printf("Total Coordinate Descent Function calc. CPU Time is: %e secs\n", gradientTotalCpuTime);
			printf("Total transformation Function calc. CPU Time is: %e secs\n", transformationTotalCpuTime);
			//==============================================================================
			// Print the Calculated Transformation Parameters At the End of Registration
			printf("Energy = %8.8lf, iteration %4d, Phi = %3.8lf, Theta = %3.8lf, Psi = %3.8lf, Sx = %2.8lf, Sy = %2.8lf, Sz = %2.8lf, Tx = %2.8lf, Ty = %2.8lf, Tz = %2.8lf\n",
				prev_energy, iteration, affineResult.rotation.x, affineResult.rotation.y, affineResult.rotation.z, affineResult.scaling.x, affineResult.scaling.y,
				affineResult.scaling.z, affineResult.translation.x, affineResult.translation.y, affineResult.translation.z);
			//==============================================================================
			for (k = 0; k < imageHeight; k++)
			{
				free(destPtr[k]);
				free(transPtr[k]);
				free(distTransPtr[k]);
			}
			free(destPtr);
			destPtr = NULL;
			//
			free(transPtr);
			transPtr = NULL;
			//
			free(distTransPtr);
			distTransPtr = NULL;
			//==============================================================================
		}
		else
		{
			// Begin Record Time
#ifdef MEASURE_TIME
			firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
			//==============================================================================
			if (components > 3)
			{
				components = 1; // Reset to the first component to start all over again
			}
			//==============================================================================
			switchcomponent = components;
			//==============================================================================
			// USE_CLIP
			if (params.use_clipbox)
			{
				// affineTmp = gradCoorDescentCompClip(destPtr, distTransPtr, 1.0, &affineResult, imageHeight, imageLength, imageWidth, switchcomponent, bestFit);
				affineTmp = gradientCoorDinateDescentComp(destPtr, distTransPtr, 1.0, &affineResult, imageHeight, imageLength, imageWidth, switchcomponent);
			}
			else
			{
				affineTmp = gradientCoorDinateDescentComp(destPtr, distTransPtr, 1.0, &affineResult, imageHeight, imageLength, imageWidth, switchcomponent);
			}
			//==============================================================================
			switch (switchcomponent)
			{
			case 1:
				// Set new values for affine temp results
				// Rotation
				affineResult.rotation.x += rotation_weight * stepsize*affineTmp.rotation.x;
				affineResult.rotation.y += rotation_weight * stepsize*affineTmp.rotation.y;
				affineResult.rotation.z += rotation_weight * stepsize*affineTmp.rotation.z;
				//==============================================================================
				components++;
				break;
			case 2:
				// Scaling
				affineResult.scaling.x += scaling_weight * stepsize*affineTmp.scaling.x;
				affineResult.scaling.y += scaling_weight * stepsize*affineTmp.scaling.y;
				affineResult.scaling.z += scaling_weight * stepsize*affineTmp.scaling.z;
				//==============================================================================
				components++;
				break;
			case 3:
				//Translation
				affineResult.translation.x += translation_weight * stepsize*affineTmp.translation.x;
				affineResult.translation.y += translation_weight * stepsize*affineTmp.translation.y;
				affineResult.translation.z += translation_weight * stepsize*affineTmp.translation.z;
				//==============================================================================
				components++;
				break;
			default:
				break;
			}
			//==============================================================================
			affineTmp_prev = affineTmp;
			//==============================================================================
			// Increase Iteration
			iteration++;
		}
	}
	//==============================================================================
#ifdef MEASURE_TIME
	regStopCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
	// Store Time For Each Registration run
	regTotalCpuTimen = regStopCpuTime - regStartCpuTime;
#endif

#ifdef CONSOLE_OUTPUT
	printf("Total Registration Function calc. CPU Time is: %e secs\n\n", regTotalCpuTimen);
#endif
	//==============================================================================
	return affineResult;
	//==============================================================================
}
//==============================================================================
Affine_Parameter registrationStochastic3D(dataType ** fixedData, dataType ** movingData, Affine_Parameter initTransform, dataType step_size, dataType tol, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType centroid[3], Registration_Params params)
{
	//==============================================================================
	// Randomly generated points will be stored here
	Random3DPoints * generated_points = malloc(sizeof(Random3DPoints) * params.rand_points);
	// The distances calculated from the random generated points
	distanceCalculated * distances_calculated = malloc(sizeof(distanceCalculated) * params.rand_points);
	// The neighbour points distance values will be stored in this structs for each of the 3 directions
	xNbDistances * xfwd_dist = malloc(sizeof(xNbDistances) * params.rand_points);
	yNbDistances * yfwd_dist = malloc(sizeof(yNbDistances) * params.rand_points);
	zNbDistances * zfwd_dist = malloc(sizeof(zNbDistances) * params.rand_points);
	// The finite differences calc. from neighbour pts differences will be stored in this struct.
	Finite_Differences * fwd_vals = malloc(sizeof(Finite_Differences) * params.rand_points);
	//==============================================================================
	int k, i, j, l, x; // change from size_t to int (because of OpenMP)
	size_t dim2D = imageLength * imageWidth, maxSurfacePts = (size_t)(0.05 * dim2D * imageHeight);
	int iteration = 0;
	const dataType h = 1.0;
	dataType firstCpuTime, secondCpuTime, regStartCpuTime, regStopCpuTime, regTotalCpuTimen = 0.;
	dataType energyTotalCpuTime = 0., distanceTotalCpuTime = 0., gradientTotalCpuTime = 0., transformationTotalCpuTime = 0., conversionTotalCpuTime = 0., surfacePtsTotalCpuTime = 0., edgeDetectionTotalCpuTime = 0., generateRandomPtsTotalCpuTime = 0., distanceCalculateTotalCpuTime = 0., finiteDiffereceTotalCpuTime = 0.;
	//==============================================================================
	const size_t mem_alloc_2D_block = sizeof(dataType)*dim2D;
	//==============================================================================
	const dataType large_value = 50000;
	//==============================================================================
	Point3D * surface_points = malloc(sizeof(Point3D) * maxSurfacePts);
	//==============================================================================
	// Affine Parameters
	Affine_Parameter affineResult;
	//==============================================================================
	// Create new fixed dist. Pointers to be used
	dataType ** destPtr = (dataType **)malloc(sizeof(dataType *) * imageHeight); // distances for destination
	//==============================================================================
	dataType ** movInitPtr = (dataType **)malloc(sizeof(dataType *) * imageHeight); // distances for Moving
	// Calc. the narrow band areas for fixed and moving dist. fn's
	dataType ** fixedNBandPtr = (dataType **)malloc(sizeof(dataType *) * imageHeight); // narrow band area for fixed dista. fn
	dataType ** movingNBandPtr = (dataType **)malloc(sizeof(dataType *) * imageHeight); // narrow band area for moving dista. fn
	dataType centroidMovingBandArea[3];
	ClipBox coordFixed, coordMoving, bestFit, coordMovingTmp;
	if (params.use_clipbox)
	{
		for (i = 0; i < imageHeight; i++)
		{
			movInitPtr[i] = (dataType *)malloc(mem_alloc_2D_block);
		}
	}
	//==============================================================================
	// Initializations of Pointers
	for (i = 0; i < imageHeight; i++)
	{
		destPtr[i] = (dataType *)malloc(mem_alloc_2D_block);
	}
	//==============================================================================
	// Instantiate Affine Parameters
	affineResult.rotation = initTransform.rotation;
	affineResult.scaling = initTransform.scaling;
	affineResult.translation = initTransform.translation;
	//==============================================================================
	// Energy tmp optimal, stop boolean
	dataType energyTmp;
	bool stopCond = false;
	//==============================================================================
	//Edge detection fn for moving binary trnasformed image
	dataType ** edgeMovingPointer = (dataType **)malloc(sizeof(dataType*) * imageHeight); // distances for Transformed Ptr
	//==============================================================================
	if (!params.use_FSM)
	{
		for (i = 0; i < imageHeight; i++)
		{
			edgeMovingPointer[i] = malloc(mem_alloc_2D_block);
		}
	}
	//==============================================================================
	// Apply distance function between transPtr and distTrans
	// Begin Record Time
#ifdef MEASURE_TIME
	firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
	//==============================================================================
	fastSweepingFunction_3D(destPtr, fixedData, imageLength, imageWidth, imageHeight, params.h, large_value, params.imageForeground);
	//==============================================================================
	if (params.use_clipbox)
	{
		//==============================================================================
		fastSweepingFunction_3D(movInitPtr, movingData, imageLength, imageWidth, imageHeight, params.h, large_value, params.imageForeground);
		//==============================================================================
	}
	//==============================================================================
#ifdef MEASURE_TIME
	secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
	// Store the time
	distanceTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
	//==============================================================================
	if (params.use_clipbox)
	{
		//==============================================================================
		// Finding the clip box points for the fixed image dist. unsigned
		coordFixed = findClipBoxSingle(destPtr, imageHeight, imageLength, imageWidth);
		//==============================================================================
		// From unsigned dist. fn
		coordMoving = findClipBoxSingle(movInitPtr, imageHeight, imageLength, imageWidth);
		//==============================================================================
		if (params.binary_nband)
		{
			for (i = 0; i < imageHeight; i++)
			{
				fixedNBandPtr[i] = (dataType *)malloc(mem_alloc_2D_block);
				movingNBandPtr[i] = (dataType *)malloc(mem_alloc_2D_block);
			}
			// Fill the narrow band areas for fixed, moving respectively
			fillNarrowBandArea(destPtr, fixedNBandPtr, imageHeight, imageLength, imageWidth, params.imageForeground, params.imageBackground);
			fillNarrowBandArea(movInitPtr, movingNBandPtr, imageHeight, imageLength, imageWidth, params.imageForeground, params.imageBackground);
			// Centroid for moving narrow band area
			centroidImage(movingNBandPtr, centroidMovingBandArea, imageHeight, imageLength, imageWidth, params.imageBackground);
		}
		//==============================================================================
		// Free after
		for (k = 0; k < imageHeight; k++)
		{
			free(movInitPtr[k]);
		}
		free(movInitPtr);
		movInitPtr = NULL;
		//==============================================================================
	}
	//==============================================================================
#ifdef CONSOLE_OUTPUT
	printf("Distance calc before Registration CPU time: %e secs\n\n", secondCpuTime - firstCpuTime);
#endif
	//==============================================================================
	// Begin Registration of Distances between shapes
	// Start Timing the Registration Process
#ifdef MEASURE_TIME
	regStartCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
	//==============================================================================
	// Create a new shape Pointers to be used
	dataType ** transPtr = (dataType **)malloc(sizeof(dataType*) * imageHeight); // Transformed Ptr
	dataType ** transMovingPtr = (dataType **)malloc(sizeof(dataType*) * imageHeight); // Transformed Ptr
	dataType ** distTransPtr = (dataType **)malloc(sizeof(dataType*) * imageHeight); // distances for Transformed Ptr
	//==============================================================================
	for (i = 0; i < imageHeight; i++)
	{
		transPtr[i] = malloc(mem_alloc_2D_block);
		transMovingPtr[i] = malloc(mem_alloc_2D_block);
		distTransPtr[i] = malloc(mem_alloc_2D_block);
	}
	//==============================================================================
	while (!stopCond)
	{
		//==============================================================================
		// Timing The Transformation Inside the registration function
#ifdef MEASURE_TIME
		firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
		//==============================================================================
		transform3DImage(movingData, transPtr, affineResult.translation, affineResult.scaling, affineResult.rotation, imageHeight, imageLength, imageWidth, params.imageBackground, centroid, params.parallelize);
		if (params.binary_nband)
		{
			// Transform the moving narrow band area
			transform3DImage(movingNBandPtr, transMovingPtr, affineResult.translation, affineResult.scaling, affineResult.rotation, imageHeight, imageLength, imageWidth, params.imageBackground, centroidMovingBandArea, params.parallelize);
			// Copy back to movingNBandPtr
			dataType ** tmpPtr = NULL;
			tmpPtr = movingNBandPtr;
			movingNBandPtr = transMovingPtr;
			transMovingPtr = tmpPtr;
		}
#ifdef MEASURE_TIME
		secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
		// Store the time
		transformationTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
		//==============================================================================
		if (params.use_clipbox)
		{
			//==============================================================================
			// Transform the coordMoving clip box using calc. transform component results
			// Copy to coordMovingTmp
			coordMovingTmp = coordMoving;
			transformClip(&coordMovingTmp, affineResult.translation, affineResult.scaling, affineResult.rotation, centroid, imageHeight, imageLength, imageWidth);
			// Find the bestFit from transformed clip
			bestFit.k_min = min(coordFixed.k_min, coordMovingTmp.k_min);
			bestFit.i_min = min(coordFixed.i_min, coordMovingTmp.i_min);
			bestFit.j_min = min(coordFixed.j_min, coordMovingTmp.j_min);

			bestFit.k_max = max(coordFixed.k_max, coordMovingTmp.k_max);
			bestFit.i_max = max(coordFixed.i_max, coordMovingTmp.i_max);
			bestFit.j_max = max(coordFixed.j_max, coordMovingTmp.j_max);
			//==============================================================================
		}
		//==============================================================================
#ifdef CONSOLE_OUTPUT
		printf("Registration Transformation calc. CPU time at iteration %4d: %e secs\n", iteration, secondCpuTime - firstCpuTime);
#endif
		//==============================================================================
		// Apply distance function between transPtr and distTrans
		// Begin Record Time
#ifdef MEASURE_TIME
		firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
		//==============================================================================
		if (params.use_clipbox)
		{
			fSweeping3D(distTransPtr, transPtr, imageLength, imageWidth, imageHeight, params.h, large_value, params.imageForeground, bestFit);
		}
		else
		{
			fastSweepingFunction_3D(distTransPtr, transPtr, imageLength, imageWidth, imageHeight, params.h, large_value, params.imageForeground);
		}
		//==============================================================================
#ifdef MEASURE_TIME
		secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
		// Store the time
		distanceTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
		//==============================================================================
#ifdef CONSOLE_OUTPUT
		printf("Distance calc during Registration CPU time at iteration %4d: %e secs\n", iteration, secondCpuTime - firstCpuTime);
#endif
		//==============================================================================
		// Evaluate Energy Function - L2 Norm Between the two calc. distances
		// Begin Record Time
#ifdef MEASURE_TIME
		firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
		//==============================================================================
		if (params.use_clipbox)
		{
			if (params.binary_nband)
			{
				// For using narrowband areas
				energyTmp = energyFunctionClipBandArea(destPtr, distTransPtr, bestFit, imageLength, fixedNBandPtr, movingNBandPtr, params.imageForeground);
			}
			else
			{
				energyTmp = energyFunctionClip(destPtr, distTransPtr, bestFit, imageLength);
			}
		}
		else
		{
			energyTmp = energyFunction(destPtr, distTransPtr, imageHeight, imageLength, imageWidth, params.h);
		}
		//==============================================================================
#ifdef MEASURE_TIME
		secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
		// Store the time
		energyTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
		//==============================================================================
#ifdef CONSOLE_OUTPUT
		printf("Energy Function calc. CPU time at iteration %4d: %e secs\n\n", iteration, secondCpuTime - firstCpuTime);
#endif
		//==============================================================================
		// Print Pre-evaluate affine values
		if (params.displayRegistrationOutputs)
		{
			printf("Energy = %5.5lf, iteration %4d, Phi = %3.5lf, Theta = %3.5lf, Psi = %3.5lf, Sx = %2.5lf, Sy = %2.5lf, Sz = %2.5lf, Tx = %2.5lf, Ty = %2.5lf, Tz = %2.5lf\n",
				energyTmp, iteration, affineResult.rotation.x, affineResult.rotation.y, affineResult.rotation.z, affineResult.scaling.x, affineResult.scaling.y,
				affineResult.scaling.z, affineResult.translation.x, affineResult.translation.y, affineResult.translation.z);
		}
		//==============================================================================
		// Check Stoping condition with tolerance and number of ierations
		if (energyTmp < tol || iteration == params.max_iterations)
		{
			//==============================================================================
			stopCond = true;
			//==============================================================================
			free(surface_points);
			//==============================================================================
			// Free up the destination ptr before exiting
			for (k = 0; k < imageHeight; k++)
			{
				free(destPtr[k]);
				if (params.binary_nband)
				{
					// Free the band areas
					free(fixedNBandPtr[k]);
					free(movingNBandPtr[k]);
				}
				if (!params.use_FSM)
				{
					// free moving edge transfromed
					free(edgeMovingPointer[k]);
				}
				//
				free(transMovingPtr[k]);
				//
				free(distTransPtr[k]);
				//
				free(transPtr[k]); // free moving transfromed
			}
			free(destPtr);
			destPtr = NULL;
			if (params.binary_nband)
			{
				// Band areas
				free(fixedNBandPtr);
				free(movingNBandPtr);
				fixedNBandPtr = NULL;
				movingNBandPtr = NULL;
			}
			if (!params.use_FSM)
			{
				// free moving edge transfromed
				free(edgeMovingPointer);
				edgeMovingPointer = NULL;
			}
			//
			free(transMovingPtr);
			transMovingPtr = NULL;
			//
			free(distTransPtr);
			distTransPtr = NULL;
			//
			free(transPtr);
			transPtr = NULL;
			//==============================================================================
			printf("Total distance Function calc. CPU Time is: %e secs\n", distanceTotalCpuTime);
			printf("Total energy Function calc. CPU Time is: %e secs\n", energyTotalCpuTime);
			printf("Total SGD Function calc. CPU Time is: %e secs\n", gradientTotalCpuTime);
			printf("Total transformation Function calc. CPU Time is: %e secs\n", transformationTotalCpuTime);
			//==============================================================================
			// Print the Calculated Transformation Parameters At the End of Registration
			printf("Energy = %5.5lf, iteration %4d, Phi = %3.5lf, Theta = %3.5lf, Psi = %3.5lf, Sx = %2.5lf, Sy = %2.5lf, Sz = %2.5lf, Tx = %2.5lf, Ty = %2.5lf, Tz = %2.5lf\n",
				energyTmp, iteration, affineResult.rotation.x, affineResult.rotation.y, affineResult.rotation.z, affineResult.scaling.x, affineResult.scaling.y,
				affineResult.scaling.z, affineResult.translation.x, affineResult.translation.y, affineResult.translation.z);
			//==============================================================================
		}
		else
		{
			//==============================================================================
			// Call th SGD Method
			//==============================================================================
			// Set up the other parameters
			//==============================================================================
			// Forward difference parameters
			dataType xFwd, yFwd, zFwd;
			// Derivative component
			dataType componentX, componentY, componentZ;
			// Stores the difference between two distance pointers
			dataType distDifference;
			// Shorter Transformation names
			dataType phi = affineResult.rotation.x, theta = affineResult.rotation.y, psi = affineResult.rotation.z;
			dataType sx = affineResult.scaling.x, sy = affineResult.scaling.y, sz = affineResult.scaling.z;
			dataType tx = affineResult.translation.x, ty = affineResult.translation.y, tz = affineResult.translation.z;
			//==============================================================================
			// Angles to radians
			dataType _cos_phi = (dataType)cos(phi), _cos_psi = (dataType)cos(psi), _cos_theta = (dataType)cos(theta);
			dataType _cos_phi_neg = (dataType)(-1 * cos(phi)), _cos_psi_neg = (dataType)(-1 * cos(psi)), _cos_theta_neg = (dataType)(-1 * cos(theta));
			dataType _sin_phi = (dataType)sin(phi), _sin_psi = (dataType)sin(psi), _sin_theta = (dataType)sin(theta);
			dataType _sin_phi_neg = (dataType)(-1 * sin(phi)), _sin_psi_neg = (dataType)(-1 * sin(psi)), _sin_theta_neg = (dataType)(-1 * sin(theta));
			//==============================================================================
			dataType _cos_phi_psi = _cos_phi * _cos_psi; // cos(phi)*cos(psi)
			dataType _cos_phi_theta = _cos_phi * _cos_theta; // cos(phi)*cos(theta)
			dataType _cos_psi_theta = _cos_psi * _cos_theta; // cos(psi)*cos(theta)
			dataType _cos_phi_theta_psi = _cos_phi_psi * _cos_theta;
			dataType _sin_psi_theta = _sin_psi * _sin_theta;// sin(psi)*sin(theta)
			dataType _sin_phi_theta = _sin_phi * _sin_theta;// sin(phi)*sin(theta)
			dataType _sin_phi_psi = _sin_phi * _sin_psi;// sin(phi)*sin(psi)
			dataType _sin_phi_psi_theta = _sin_phi_psi * _sin_theta; // sin(phi)*sin(psi)*sin(theta)
			//==============================================================================
			dataType _sin_phi_neg_sin_psi = _sin_phi_neg * _sin_psi; // (-sin(phi)*sin(psi)
			dataType _cos_phi_psi_sin_theta = _cos_phi_psi * _sin_theta; // cos(phi)*cos(psi)*sin(theta)
			dataType _cos_psi_neg_sin_phi = _cos_psi_neg * _sin_phi; // (-cos(psi)*sin(phi)
			dataType _cos_phi_neg_sin_psi_theta = _cos_phi_neg * _sin_psi_theta;
			dataType _cos_phi_sin_psi_theta = _cos_phi * _sin_psi_theta; // cos(phi)*sin(psi)*sin(theta)
			dataType _cos_phi_sin_psi = _cos_phi * _sin_psi;// cos(phi)*sin(psi)
			dataType _cos_psi_sin_phi_theta = _cos_psi * _sin_phi_theta;// cos(psi)*sin(phi)*sin(theta)
			dataType _cos_theta_sin_phi = _cos_theta * _sin_phi; // cos(theta) * sin(phi)
			//==============================================================================
			// 2nd component Ry
			dataType _cos_psi_sin_theta = _cos_psi * _sin_theta;// cos(psi)*sin(theta)
			dataType _cos_psi_theta_sin_phi = _cos_psi_theta * _sin_phi; // cos(psi)*cos(theta)*sin(phi)
			dataType _cos_theta_sin_phi_psi = _cos_theta * _sin_phi_psi; // cos(theta)*sin(phi)*sin(psi)
			dataType _cos_phi_theta_sin_psi = _cos_phi_theta * _sin_psi;// cos(phi)*cos(theta)*sin(psi)
			dataType _cos_phi_sin_theta = _cos_phi * _sin_theta; // cos(phi)*sin(theta)
			//==============================================================================
			// 3nd component Rz
			dataType _cos_theta_sin_psi = _cos_theta * _sin_psi; // cos(theta) * sin(psi)
			dataType _cos_phi_neg_sin_psi = _cos_phi_neg * _sin_psi;// -cos(phi)*sin(psi)
			//==============================================================================
			// Set the scales - Directional
			dataType _cos_psi__sin_phi_theta = _cos_psi * _sin_phi_theta;
			//==============================================================================
			dataType _cos_psi_sin_phi = _cos_psi * _sin_phi;// cos(psi)*sin(phi)
			//==============================================================================
			dataType pVal, qVal; // Used in Finite differences
			dataType hh;
			//==============================================================================
			// Setting same if uniform scale and rotation
#ifdef UNIFORM
			sx = sy = sz;
			phi = theta = psi;
#endif // UNIFORM
			//==============================================================================
			// Scale denominator
			dataType inv_sx2 = (dataType)(1.0 / (sx*sx));
			dataType inv_sy2 = (dataType)(1.0 / (sy*sy));
			dataType inv_sz2 = (dataType)(1.0 / (sz*sz));
			//==============================================================================
			// Randomly generated x,y,z points
			//==============================================================================
			// Randomizing random number generation
			//srand(time(NULL));
			//==============================================================================
			bool loop = true;
			//==============================================================================
#ifdef MEASURE_TIME
			firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
			if (!params.use_FSM)
			{
				edgeDetection3dFunctionD(transPtr, edgeMovingPointer, imageLength, imageWidth, imageHeight, params.imageBackground, params.imageForeground, params.insideShapevalue);
			}
#ifdef MEASURE_TIME
			secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
			// Store the time recorded
			edgeDetectionTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
			//==============================================================================
			dataType getDist;
			//==============================================================================
			size_t k_2, i_2;//loop counter for z dimension
			//==============================================================================
			size_t ptsNum = 0;
			//==============================================================================
			size_t k_min, k_max, i_min, i_max, j_min, j_max;
			size_t i_2_max, i_2_min;
			//==============================================================================
			// Begin record time
			// Calc. the points and store them
#ifdef MEASURE_TIME
			firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
			if (!params.use_FSM)
			{
				//==============================================================================
				k_min = bestFit.k_min, k_max = bestFit.k_max + 1, i_min = bestFit.i_min, i_max = bestFit.i_max, j_min = bestFit.j_min, j_max = bestFit.j_max;
				i_2_max = x_new(i_max, j_max, imageLength), i_2_min = x_new(i_min, j_min, imageLength);
				//==============================================================================
				Point3D * tmpPt;
				//==============================================================================
				// Find the surcae points only within the clipbox
				for (k_2 = k_min; k_2 < k_max; k_2++) // z axis of the input surface or image
				{
					for (i_2 = i_2_min; i_2 < i_2_max; i_2++)// x-y axis of the input surface or image
					{
						//==============================================================================
						if (edgeMovingPointer[k_2][i_2] == params.imageForeground)
						{
							//==============================================================================
							tmpPt = &surface_points[ptsNum];
							//==============================================================================
							tmpPt->x = (dataType)(i_2 / imageLength);
							tmpPt->y = (dataType)(i_2 % imageLength);
							tmpPt->z = (dataType)k_2;
							ptsNum++;
						}
					}
				}
			}
			//==============================================================================
#ifdef MEASURE_TIME
			secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
			// Store the time recorded
			surfacePtsTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
			// Pass this surface_points to get distance fn
			//==============================================================================
			// Generate and store the random points
#ifdef MEASURE_TIME
			firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
			if (!params.use_FSM)
			{
				l = 0;
				Random3DPoints * tmpRdPts;
				do
				{
					// Generate random points
					if (params.use_clipbox)
					{
						// From clipbox
						k = bestFit.k_min + (rand() % bestFit.k_max), i = bestFit.i_min + (rand() % bestFit.i_max), j = bestFit.j_min + (rand() % bestFit.j_max);
					}
					else
					{
						// From whole domain
						k = rand() % imageHeight, i = rand() % imageLength, j = rand() % imageWidth;
					}
					// Check if inside the narrow band
					x = x_new(i, j, imageLength);
					if (params.binary_nband)
					{
						if (NFunctionBinary(fixedNBandPtr[k][x], movingNBandPtr[k][x], params.imageForeground) == 1) // Checks if inside the narrow band areas for random generated points
						{
							// Start loop

							//==============================================================================
							tmpRdPts = &generated_points[l];
							tmpRdPts->k = k;
							tmpRdPts->i = i;
							tmpRdPts->j = j;
							//==============================================================================
							l = l + 1;
							//==============================================================================
							// Check condition
							if (l == params.rand_points)
							{
								loop = false;
							}
						}
					}
					else
					{
						if (NFunction(destPtr[k][x], distTransPtr[k][x], params.imageForeground) == 1) // Checks if inside the narrow band areas for random generated points
						{
							// Start loop

							//==============================================================================
							tmpRdPts = &generated_points[l];
							tmpRdPts->k = k;
							tmpRdPts->i = i;
							tmpRdPts->j = j;
							//==============================================================================
							l = l + 1;
							//==============================================================================
							// Check condition
							if (l == params.rand_points)
							{
								loop = false;
							}
						}
					}
				} while ((loop) && (l <= params.rand_points));
			}
			//==============================================================================
#ifdef MEASURE_TIME
			secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
			// Store the time recorded
			generateRandomPtsTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
			//==============================================================================
			// Calculate the distances from generated points
			distanceCalculated * tmpDistances;
			xNbDistances * tmpXFwd;
			yNbDistances * tmpYFwd;
			zNbDistances * tmpZFwd;
			dataType destFixed;
#ifdef MEASURE_TIME
			firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
			//==============================================================================
			if (!params.use_FSM)
			{
				if (params.parallelize)
				{
					//==============================================================================
					// Parallelize
					omp_set_dynamic(0); // Disable dynamic adjustment of threads
					//omp_set_num_threads(omp_num_procs()); // Request as many threads as you have processors
					omp_set_num_threads(NUM_THREADS); // Request as many threads as you have processors
#pragma omp parallel
					{
#pragma omp for private(k , i , j, x, tmpDistances, tmpXFwd, tmpYFwd, tmpZFwd, getDist, destFixed) firstprivate(hh, pVal, qVal) schedule(static) nowait
						for (l = 0; l < params.rand_points; l++)
						{
							//==============================================================================
							k = generated_points[l].k;
							i = generated_points[l].i;
							j = generated_points[l].j;
							//==============================================================================
							// 2D to 1D representation for i, j
							x = x_new(i, j, imageLength);
							//==============================================================================
							tmpDistances = &distances_calculated[l];
							//==============================================================================
							tmpXFwd = &xfwd_dist[l];
							tmpYFwd = &yfwd_dist[l];
							tmpZFwd = &zfwd_dist[l];
							//==============================================================================
							getDist = getDistance(edgeMovingPointer, imageHeight, imageLength, dim2D, k, x, bestFit, surface_points, ptsNum, params.insideShapevalue, params.parallelize);
							destFixed = destPtr[k][x];
							tmpDistances->distDifference = (dataType)((destFixed - getDist) * 2.0);
							//==============================================================================
							nbPointsX(edgeMovingPointer, h, x, k, i, imageHeight, imageLength, imageWidth, &pVal, &qVal, &hh, bestFit, surface_points, ptsNum, params.insideShapevalue, params.parallelize);
							tmpXFwd->hvl = hh;
							tmpXFwd->pvl = pVal;
							tmpXFwd->qvl = qVal;
							//==============================================================================
							nbPointsY(edgeMovingPointer, h, k, i, j, imageHeight, imageLength, imageWidth, &pVal, &qVal, &hh, bestFit, surface_points, ptsNum, params.insideShapevalue, params.parallelize);
							tmpYFwd->hvl = hh;
							tmpYFwd->pvl = pVal;
							tmpYFwd->qvl = qVal;
							//==============================================================================
							nbPointsZ(edgeMovingPointer, h, x, k, imageHeight, imageLength, imageWidth, &pVal, &qVal, &hh, bestFit, surface_points, ptsNum, params.insideShapevalue, params.parallelize);
							tmpZFwd->hvl = hh;
							tmpZFwd->pvl = pVal;
							tmpZFwd->qvl = qVal;
							//==============================================================================
						}
						//==============================================================================
					}
				}
				else
				{
					//==============================================================================
					// Sequential
					for (l = 0; l < params.rand_points; l++)
					{
						k = generated_points[l].k;
						i = generated_points[l].i;
						j = generated_points[l].j;
						//==============================================================================
						// 2D to 1D representation for i, j
						x = x_new(i, j, imageLength);
						//==============================================================================
						tmpDistances = &distances_calculated[l];
						//==============================================================================
						tmpXFwd = &xfwd_dist[l];
						tmpYFwd = &yfwd_dist[l];
						tmpZFwd = &zfwd_dist[l];
						//==============================================================================
						getDist = getDistance(edgeMovingPointer, imageHeight, imageLength, dim2D, k, x, bestFit, surface_points, ptsNum, params.insideShapevalue, params.parallelize);
						destFixed = destPtr[k][x];
						tmpDistances->distDifference = (dataType)((destFixed - getDist) * 2.0);
						//==============================================================================
						nbPointsX(edgeMovingPointer, h, x, k, i, imageHeight, imageLength, imageWidth, &pVal, &qVal, &hh, bestFit, surface_points, ptsNum, params.insideShapevalue, params.parallelize);
						tmpXFwd->hvl = hh;
						tmpXFwd->pvl = pVal;
						tmpXFwd->qvl = qVal;
						//==============================================================================
						nbPointsY(edgeMovingPointer, h, k, i, j, imageHeight, imageLength, imageWidth, &pVal, &qVal, &hh, bestFit, surface_points, ptsNum, params.insideShapevalue, params.parallelize);
						tmpYFwd->hvl = hh;
						tmpYFwd->pvl = pVal;
						tmpYFwd->qvl = qVal;
						//==============================================================================
						nbPointsZ(edgeMovingPointer, h, x, k, imageHeight, imageLength, imageWidth, &pVal, &qVal, &hh, bestFit, surface_points, ptsNum, params.insideShapevalue, params.parallelize);
						tmpZFwd->hvl = hh;
						tmpZFwd->pvl = pVal;
						tmpZFwd->qvl = qVal;
						//==============================================================================
					}
					//==============================================================================
				}
			}
#ifdef MEASURE_TIME
			secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
			// Store the time recorded
			distanceCalculateTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
			//==============================================================================
			Finite_Differences * tmpFwd;
			dataType pv_tmp, qv_tmp, h_tmp;
			// Cacluate the finite differences
#ifdef MEASURE_TIME
			firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
			if (!params.use_FSM)
			{
				for (l = 0; l < params.rand_points; l++)
				{
					tmpFwd = &fwd_vals[l];
					//==============================================================================
					pv_tmp = xfwd_dist[l].pvl;
					qv_tmp = xfwd_dist[l].qvl;
					h_tmp = xfwd_dist[l].hvl;
					//==============================================================================
					tmpFwd->xFwd = finiteDifAll(pv_tmp, qv_tmp, h_tmp);
					//==============================================================================
					pv_tmp = yfwd_dist[l].pvl;
					qv_tmp = yfwd_dist[l].qvl;
					h_tmp = yfwd_dist[l].hvl;
					//==============================================================================
					tmpFwd->yFwd = finiteDifAll(pv_tmp, qv_tmp, h_tmp);
					//==============================================================================
					pv_tmp = zfwd_dist[l].pvl;
					qv_tmp = zfwd_dist[l].qvl;
					h_tmp = zfwd_dist[l].hvl;
					//==============================================================================
					tmpFwd->zFwd = finiteDifAll(pv_tmp, qv_tmp, h_tmp);
					//==============================================================================
				}
			}
			//==============================================================================
#ifdef MEASURE_TIME
			secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
			// Store the time recorded
			finiteDiffereceTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
			//==============================================================================
			// Calculate the gradient components
#ifdef MEASURE_TIME
			firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
			if (params.use_FSM)
			{
				l = 0;
				do
				{
					// Generate random points
					if (params.use_clipbox)
					{
						// From clipbox
						k = bestFit.k_min + (rand() % bestFit.k_max), i = bestFit.i_min + (rand() % bestFit.i_max), j = bestFit.j_min + (rand() % bestFit.j_max);
					}
					else
					{
						// From whole domain
						k = rand() % imageHeight, i = rand() % imageLength, j = rand() % imageWidth;
					}
					// 2D to 1D representation for i, j
					x = x_new(i, j, imageLength);
					// Check if inside the narrow band
					if (params.binary_nband)
					{
						if (NFunctionBinary(fixedNBandPtr[k][x], movingNBandPtr[k][x], params.imageForeground) == 1) // Checks if inside the narrow band areas for random generated points
						{
							// Start loop
							l = l + 1;
							//==============================================================================
							// Store the distance function difference
							distDifference = (dataType)((destPtr[k][x] - distTransPtr[k][x]) * 2.0);
							//==============================================================================
							// Directional component vector derivatives - i, j, k
							dataType tmpI = i / (dataType)imageLength, tmpJ = j / (dataType)imageWidth, tmpK = k / (dataType)imageHeight;
							//==============================================================================
							// Apply Forward Differences to the distTrans pointer
							xFwd = finiteDifX(distTransPtr, h, x, k, i, imageLength);
							yFwd = finiteDifY(distTransPtr, h, k, i, j, imageLength, imageWidth);
							zFwd = finiteDifZ(distTransPtr, h, x, k, i, imageLength, imageHeight);
							//==============================================================================
							// Evaluate Individual Gradient Components
							//==============================================================================
#ifdef DIRECTIONAL
							// Set the Rotation - Directional
							//==============================================================================
							componentX = yFwd * (((tmpI)*((_sin_phi_neg_sin_psi + _cos_phi_psi_sin_theta) / sy)) + ((tmpJ)*((_cos_psi_neg_sin_phi - _cos_phi_sin_psi_theta) / sy)) + ((-tmpK)*((_cos_phi_theta) / sy))) +
								zFwd * (((tmpI)*((_cos_phi_sin_psi + _cos_psi_sin_phi_theta) / sz)) + ((tmpJ)*((_cos_phi_psi - _sin_phi_psi_theta) / sz)) + ((-tmpK)*((_cos_theta_sin_phi) / sz)));
							//==============================================================================
							componentY = xFwd * (((-tmpI)*((_cos_psi_sin_theta) / sx)) + ((tmpJ)*((_sin_psi_theta) / sx)) + ((tmpK)*((_cos_theta) / sx))) +
								yFwd * (((tmpI)*((_cos_psi_theta_sin_phi) / sy)) + ((-tmpJ)*((_cos_theta_sin_phi_psi) / sy)) + ((tmpK)*((_sin_phi_theta) / sy))) +
								zFwd * (((-tmpI)*((_cos_phi_theta_psi) / sz)) + ((tmpJ)*((_cos_phi_theta_sin_psi) / sz)) + ((-tmpK)*((_cos_phi_sin_theta) / sz)));
							//==============================================================================
							componentZ = xFwd * (((-tmpI)*((_cos_theta_sin_psi) / sx)) + ((-tmpJ)*((_cos_psi_theta) / sx))) +
								yFwd * (((tmpI)*((_cos_phi_psi - _sin_phi_psi_theta) / sy)) + ((tmpJ)*((_cos_phi_neg_sin_psi - _cos_psi_sin_phi_theta) / sy))) +
								zFwd * (((tmpI)*((_cos_psi_sin_phi + _cos_phi_sin_psi_theta) / sz)) + ((tmpJ)*((_sin_phi_neg_sin_psi + _cos_phi_psi_sin_theta) / sz)));
							//==============================================================================
							// Set the Rotations - Directional
							affineResult.rotation.x += (params.rotation_weight * (step_size)*(componentX)*(distDifference)) / params.rand_points;
							affineResult.rotation.y += (params.rotation_weight * (step_size)*(componentY)*(distDifference)) / params.rand_points;
							affineResult.rotation.z += (params.rotation_weight * (step_size)*(componentZ)*(distDifference)) / params.rand_points;
							//==============================================================================
							// Set the scales - Directional
							componentX = xFwd * ((-tmpI)*((_cos_psi_theta) / (sx*sx)) + ((tmpJ)*((_cos_theta_sin_psi) / (sx*sx))) + ((-tmpK)*((_sin_theta) / (sx*sx))));
							//==============================================================================
							componentY = yFwd * (((-tmpI)*((_cos_phi_sin_psi + _cos_psi__sin_phi_theta) / (sy*sy))) + ((-tmpJ)*((_cos_phi_psi - _sin_phi_psi_theta) / (sy*sy))) + ((tmpK)*((_cos_theta_sin_phi) / (sy*sy))));
							//==============================================================================
							componentZ = zFwd * (((-tmpI)*((_sin_phi_psi - _cos_phi_psi_sin_theta) / (sz*sz))) + ((-tmpJ)*((_cos_psi_sin_phi + _cos_phi_sin_psi_theta) / (sz*sz))) + ((-tmpK)*((_cos_phi_theta) / (sz*sz))));
							//==============================================================================
							// Set the Scales - Directional Scales
							affineResult.scaling.x += (params.scaling_weight * (step_size)*(componentX)*(distDifference)) / params.rand_points;
							affineResult.scaling.y += (params.scaling_weight * (step_size)*(componentY)*(distDifference)) / params.rand_points;
							affineResult.scaling.z += (params.scaling_weight * (step_size)*(componentZ)*(distDifference)) / params.rand_points;
#endif // DIRECTIONAL
							//==============================================================================
							// Translation Parameters - Always DIRECTIONAL
							// Tx
							affineResult.translation.x += (params.translation_weight * step_size*(-xFwd)*(distDifference)) / params.rand_points;
							// Ty
							affineResult.translation.y += (params.translation_weight * step_size*(-yFwd)*(distDifference)) / params.rand_points;
							// Tz
							affineResult.translation.z += (params.translation_weight * step_size*(-zFwd)*(distDifference)) / params.rand_points;
							//==============================================================================
							// Check condition
							if (l == params.rand_points)
							{
								loop = false;
							}
						}
					}
					else
					{
						if (NFunction(destPtr[k][x], distTransPtr[k][x], NDelta) == 1) // Checks if inside the band for random generated points
						//if (NFunctionBinary(fixedNBandPtr[k][x], movingNBandPtr[k][x], foregound) == 1) // Checks if inside the narrow band areas for random generated points
						{
							// Start loop
							l = l + 1;
							//==============================================================================
							// Store the distance function difference
							distDifference = (dataType)((destPtr[k][x] - distTransPtr[k][x]) * 2.0);
							//==============================================================================
							// Directional component vector derivatives - i, j, k
							dataType tmpI = i / (dataType)imageLength, tmpJ = j / (dataType)imageWidth, tmpK = k / (dataType)imageHeight;
							//==============================================================================
							// Apply Forward Differences to the distTrans pointer
							xFwd = finiteDifX(distTransPtr, h, x, k, i, imageLength);
							yFwd = finiteDifY(distTransPtr, h, k, i, j, imageLength, imageWidth);
							zFwd = finiteDifZ(distTransPtr, h, x, k, i, imageLength, imageHeight);
							//==============================================================================
							// Evaluate Individual Gradient Components
							//==============================================================================
#ifdef DIRECTIONAL
							// Set the Rotation - Directional
							//==============================================================================
							componentX = yFwd * (((tmpI)*((_sin_phi_neg_sin_psi + _cos_phi_psi_sin_theta) / sy)) + ((tmpJ)*((_cos_psi_neg_sin_phi - _cos_phi_sin_psi_theta) / sy)) + ((-tmpK)*((_cos_phi_theta) / sy))) +
								zFwd * (((tmpI)*((_cos_phi_sin_psi + _cos_psi_sin_phi_theta) / sz)) + ((tmpJ)*((_cos_phi_psi - _sin_phi_psi_theta) / sz)) + ((-tmpK)*((_cos_theta_sin_phi) / sz)));
							//==============================================================================
							componentY = xFwd * (((-tmpI)*((_cos_psi_sin_theta) / sx)) + ((tmpJ)*((_sin_psi_theta) / sx)) + ((tmpK)*((_cos_theta) / sx))) +
								yFwd * (((tmpI)*((_cos_psi_theta_sin_phi) / sy)) + ((-tmpJ)*((_cos_theta_sin_phi_psi) / sy)) + ((tmpK)*((_sin_phi_theta) / sy))) +
								zFwd * (((-tmpI)*((_cos_phi_theta_psi) / sz)) + ((tmpJ)*((_cos_phi_theta_sin_psi) / sz)) + ((-tmpK)*((_cos_phi_sin_theta) / sz)));
							//==============================================================================
							componentZ = xFwd * (((-tmpI)*((_cos_theta_sin_psi) / sx)) + ((-tmpJ)*((_cos_psi_theta) / sx))) +
								yFwd * (((tmpI)*((_cos_phi_psi - _sin_phi_psi_theta) / sy)) + ((tmpJ)*((_cos_phi_neg_sin_psi - _cos_psi_sin_phi_theta) / sy))) +
								zFwd * (((tmpI)*((_cos_psi_sin_phi + _cos_phi_sin_psi_theta) / sz)) + ((tmpJ)*((_sin_phi_neg_sin_psi + _cos_phi_psi_sin_theta) / sz)));
							//==============================================================================
							// Set the Rotations - Directional
							affineResult.rotation.x += (params.rotation_weight * (step_size)*(componentX)*(distDifference)) / params.rand_points;
							affineResult.rotation.y += (params.rotation_weight * (step_size)*(componentY)*(distDifference)) / params.rand_points;
							affineResult.rotation.z += (params.rotation_weight * (step_size)*(componentZ)*(distDifference)) / params.rand_points;
							//==============================================================================
							// Set the scales - Directional
							componentX = xFwd * ((-tmpI)*((_cos_psi_theta) / (sx*sx)) + ((tmpJ)*((_cos_theta_sin_psi) / (sx*sx))) + ((-tmpK)*((_sin_theta) / (sx*sx))));
							//==============================================================================
							componentY = yFwd * (((-tmpI)*((_cos_phi_sin_psi + _cos_psi__sin_phi_theta) / (sy*sy))) + ((-tmpJ)*((_cos_phi_psi - _sin_phi_psi_theta) / (sy*sy))) + ((tmpK)*((_cos_theta_sin_phi) / (sy*sy))));
							//==============================================================================
							componentZ = zFwd * (((-tmpI)*((_sin_phi_psi - _cos_phi_psi_sin_theta) / (sz*sz))) + ((-tmpJ)*((_cos_psi_sin_phi + _cos_phi_sin_psi_theta) / (sz*sz))) + ((-tmpK)*((_cos_phi_theta) / (sz*sz))));
							//==============================================================================
							// Set the Scales - Directional Scales
							affineResult.scaling.x += (params.scaling_weight * (step_size)*(componentX)*(distDifference)) / params.rand_points;
							affineResult.scaling.y += (params.scaling_weight * (step_size)*(componentY)*(distDifference)) / params.rand_points;
							affineResult.scaling.z += (params.scaling_weight * (step_size)*(componentZ)*(distDifference)) / params.rand_points;
#endif // DIRECTIONAL
							//==============================================================================
							// Translation Parameters - Always DIRECTIONAL
							// Tx
							affineResult.translation.x += (params.translation_weight * step_size*(-xFwd)*(distDifference)) / params.rand_points;
							// Ty
							affineResult.translation.y += (params.translation_weight * step_size*(-yFwd)*(distDifference)) / params.rand_points;
							// Tz
							affineResult.translation.z += (params.translation_weight * step_size*(-zFwd)*(distDifference)) / params.rand_points;
							//==============================================================================
							// Check condition
							if (l == params.rand_points)
							{
								loop = false;
							}
						}
					}
				} while ((loop) && (l <= params.rand_points));
			}
			else
			{
				for (l = 0; l < params.rand_points; l++)
				{
					//==============================================================================
					yFwd = fwd_vals[l].yFwd;
					xFwd = fwd_vals[l].xFwd;
					zFwd = fwd_vals[l].zFwd;
					//==============================================================================
					distDifference = distances_calculated[l].distDifference;
					//==============================================================================
					k = generated_points[l].k;
					i = generated_points[l].i;
					j = generated_points[l].j;
					//==============================================================================
					// Directional component vector derivatives - i, j, k
					dataType tmpI = i / (dataType)imageLength, tmpJ = j / (dataType)imageWidth, tmpK = k / (dataType)imageHeight;
					//==============================================================================
					// Shorten the radian calcultaion
					//==============================================================================
#ifdef DIRECTIONAL
				// Set the Rotation - Directional
				//==============================================================================
					componentX = yFwd * (((tmpI)*((_sin_phi_neg_sin_psi + _cos_phi_psi_sin_theta) / sy)) + ((tmpJ)*((_cos_psi_neg_sin_phi - _cos_phi_sin_psi_theta) / sy)) + ((-tmpK)*((_cos_phi_theta) / sy))) +
						zFwd * (((tmpI)*((_cos_phi_sin_psi + _cos_psi_sin_phi_theta) / sz)) + ((tmpJ)*((_cos_phi_psi - _sin_phi_psi_theta) / sz)) + ((-tmpK)*((_cos_theta_sin_phi) / sz)));
					//==============================================================================
					componentY = xFwd * (((-tmpI)*((_cos_psi_sin_theta) / sx)) + ((tmpJ)*((_sin_psi_theta) / sx)) + ((tmpK)*((_cos_theta) / sx))) +
						yFwd * (((tmpI)*((_cos_psi_theta_sin_phi) / sy)) + ((-tmpJ)*((_cos_theta_sin_phi_psi) / sy)) + ((tmpK)*((_sin_phi_theta) / sy))) +
						zFwd * (((-tmpI)*((_cos_phi_theta_psi) / sz)) + ((tmpJ)*((_cos_phi_theta_sin_psi) / sz)) + ((-tmpK)*((_cos_phi_sin_theta) / sz)));
					//==============================================================================
					componentZ = xFwd * (((-tmpI)*((_cos_theta_sin_psi) / sx)) + ((-tmpJ)*((_cos_psi_theta) / sx))) +
						yFwd * (((tmpI)*((_cos_phi_psi - _sin_phi_psi_theta) / sy)) + ((tmpJ)*((_cos_phi_neg_sin_psi - _cos_psi_sin_phi_theta) / sy))) +
						zFwd * (((tmpI)*((_cos_psi_sin_phi + _cos_phi_sin_psi_theta) / sz)) + ((tmpJ)*((_sin_phi_neg_sin_psi + _cos_phi_psi_sin_theta) / sz)));
					//==============================================================================
					// Set the Rotations - Directional
					affineResult.rotation.x += (params.rotation_weight * (step_size)*(componentX)*(distDifference)) / params.rand_points;
					affineResult.rotation.y += (params.rotation_weight * (step_size)*(componentY)*(distDifference)) / params.rand_points;
					affineResult.rotation.z += (params.rotation_weight * (step_size)*(componentZ)*(distDifference)) / params.rand_points;
					//==============================================================================
					componentX = xFwd * ((-tmpI)*((_cos_psi_theta) / (sx*sx)) + ((tmpJ)*((_cos_theta_sin_psi) / (sx*sx))) + ((-tmpK)*((_sin_theta) / (sx*sx))));
					//==============================================================================
					componentY = yFwd * (((-tmpI)*((_cos_phi_sin_psi + _cos_psi__sin_phi_theta) / (sy*sy))) + ((-tmpJ)*((_cos_phi_psi - _sin_phi_psi_theta) / (sy*sy))) + ((tmpK)*((_cos_theta_sin_phi) / (sy*sy))));
					//==============================================================================
					componentZ = zFwd * (((-tmpI)*((_sin_phi_psi - _cos_phi_psi_sin_theta) / (sz*sz))) + ((-tmpJ)*((_cos_psi_sin_phi + _cos_phi_sin_psi_theta) / (sz*sz))) + ((-tmpK)*((_cos_phi_theta) / (sz*sz))));
					//==============================================================================
					// Set the Scales - Directional Scales
					affineResult.scaling.x += (params.scaling_weight * (step_size)*(componentX)*(distDifference)) / params.rand_points;
					affineResult.scaling.y += (params.scaling_weight * (step_size)*(componentY)*(distDifference)) / params.rand_points;
					affineResult.scaling.z += (params.scaling_weight * (step_size)*(componentZ)*(distDifference)) / params.rand_points;
#endif // DIRECTIONAL
					//==============================================================================
					// Translation Parameters - Always DIRECTIONAL
					// Tx
					affineResult.translation.x += (params.translation_weight * step_size*(-xFwd)*(distDifference)) / params.rand_points;
					// Ty
					affineResult.translation.y += (params.translation_weight * step_size*(-yFwd)*(distDifference)) / params.rand_points;
					// Tz
					affineResult.translation.z += (params.translation_weight * step_size*(-zFwd)*(distDifference)) / params.rand_points;
				}
			}
			//==============================================================================
#ifdef MEASURE_TIME
			secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
			// Store the time
			gradientTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
			//==============================================================================
#ifdef CONSOLE_OUTPUT
			printf("SGD Components calc. CPU time at iteration %4d: %e secs\n\n", iteration, secondCpuTime - firstCpuTime);
#endif
			//==============================================================================
			// Increase Iteration
			iteration++;
			//==============================================================================
		}
	}
	//==============================================================================
	// Stop Timing the Registration Process
#ifdef MEASURE_TIME
	regStopCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
	// Store Time For Each Registration run
	regTotalCpuTimen = regStopCpuTime - regStartCpuTime;
#endif
	//==============================================================================
#ifdef CONSOLE_OUTPUT
	printf("Total Registration Function calc. CPU Time is: %e secs\n\n", regTotalCpuTimen);
#endif
	//==============================================================================
	return affineResult;
	//==============================================================================
}
//==============================================================================
ClipBox findClipBoxSingle(dataType ** Source, size_t imageHeight, size_t imageLength, size_t imageWidth)
{
	ClipBox coord;
	coord.k_min = imageHeight, coord.i_min = imageLength, coord.j_min = imageWidth, coord.k_max = 0, coord.i_max = 0, coord.j_max = 0;
	size_t k, i, j;
	for (k = 0; k < imageHeight; k++)
	{
		for (i = 0; i < imageLength; i++)
		{
			for (j = 0; j < imageWidth; j++)
			{
				// 2D to 1D representation for i, j
				size_t x = x_new(i, j, imageLength);
				if (NFunctionOne(Source[k][x], NDelta) == 1)
				{
					// Find K clip
					if (k > coord.k_max)
					{
						coord.k_max = k;
					}
					if (k < coord.k_min)
					{
						coord.k_min = k;
					}
					// FInd I clip
					if (i > coord.i_max)
					{
						coord.i_max = i;
					}
					if (i < coord.i_min)
					{
						coord.i_min = i;
					}
					// Find J clip
					if (j > coord.j_max)
					{
						coord.j_max = j;
					}
					if (j < coord.j_min)
					{
						coord.j_min = j;
					}
				}
			}
		}
	}
	return coord;
}
//==============================================================================
void fillNarrowBandArea(dataType ** sourceDist, dataType ** bandContainer, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType insideValue, dataType outsideValue)
{
	size_t k, i, j;
	for (k = 0; k < imageHeight; k++)
	{
		for (i = 0; i < imageLength; i++)
		{
			for (j = 0; j < imageWidth; j++)
			{
				// 2D to 1D representation for i, j
				size_t xd = x_new(i, j, imageLength);
				if (NFunctionOne(sourceDist[k][xd], NDelta) == 1)
				{
					bandContainer[k][xd] = insideValue;
				}
				else
				{
					bandContainer[k][xd] = outsideValue;
				}
			}
		}
	}
}
//==============================================================================
ClipBox findClipBoxTwo(dataType ** destination, dataType ** source, size_t imageHeight, size_t imageLength, size_t imageWidth)
{
	ClipBox coord;
	coord.k_min = imageHeight, coord.i_min = imageLength, coord.j_min = imageWidth, coord.k_max = 0, coord.i_max = 0, coord.j_max = 0;
	size_t k, i, j;
	for (k = 0; k < imageHeight; k++)
	{
		for (i = 0; i < imageLength; i++)
		{
			for (j = 0; j < imageWidth; j++)
			{
				// 2D to 1D representation for i, j
				size_t x = x_new(i, j, imageLength);
				if (NFunction(destination[k][x], source[k][x], NDelta) == 1)
				{
					// Find K clip
					if (k > coord.k_max)
					{
						coord.k_max = k;
					}
					else if (k < coord.k_min)
					{
						coord.k_min = k;
					}
					// FInd I clip
					if (i > coord.i_max)
					{
						coord.i_max = i;
					}
					else if (i < coord.i_min)
					{
						coord.i_min = i;
					}
					// Find J clip
					if (j > coord.j_max)
					{
						coord.j_max = j;
					}
					else if (j < coord.j_min)
					{
						coord.j_min = j;
					}
				}
			}
		}
	}
	return coord;
}
//==============================================================================
void transformClip(ClipBox *bestfit, Point3D translation, Point3D scaling, Point3D rotation, dataType centroid[3], size_t imageHeight, size_t imageLength, size_t imageWidth)
{
	//==============================================================================
	size_t k_min = (*bestfit).k_min, k_max = (*bestfit).k_max;
	size_t i_min = (*bestfit).i_min, i_max = (*bestfit).i_max;
	size_t j_min = (*bestfit).j_min, j_max = (*bestfit).j_max;
	//==============================================================================
	// Floor Points - k min
	// i_min, j_min, k_min; // Corner 1
	// i_max, j_min, k_min; // Corner 2
	// i_min, j_max, k_min; // Corner 3
	// i_max, j_max, k_min; // Corner 4
	//==============================================================================
	// Ceiling Points - k max
	// i_min, j_min, k_max; // Corner 5
	// i_max, j_min, k_max; // Corner 6
	// i_min, j_max, k_max; // Corner 7
	// i_max, j_max, k_max; // Corner 8
	//==============================================================================
	// Transform the bottom corner point
	// Corner 1
	CoordPoints c1 = { k_min, i_min, j_min };
	c1 = transformPoint(&c1, translation, scaling, rotation, centroid, imageHeight, imageLength, imageWidth, 1);
	CoordPoints c2 = { k_min, i_max, j_min };
	c2 = transformPoint(&c2, translation, scaling, rotation, centroid, imageHeight, imageLength, imageWidth, 1);
	CoordPoints c3 = { k_min, i_min, j_max };
	c3 = transformPoint(&c3, translation, scaling, rotation, centroid, imageHeight, imageLength, imageWidth, 1);
	CoordPoints c4 = { k_min, i_max, j_max };
	c4 = transformPoint(&c4, translation, scaling, rotation, centroid, imageHeight, imageLength, imageWidth, 1);
	//==============================================================================
	//transform the top corner points
	// Corner 1
	CoordPoints c5 = { k_max, i_min, j_min };
	c5 = transformPoint(&c5, translation, scaling, rotation, centroid, imageHeight, imageLength, imageWidth, 2);
	CoordPoints c6 = { k_max, i_max, j_min };
	c6 = transformPoint(&c6, translation, scaling, rotation, centroid, imageHeight, imageLength, imageWidth, 2);
	CoordPoints c7 = { k_max, i_min, j_max };
	c7 = transformPoint(&c7, translation, scaling, rotation, centroid, imageHeight, imageLength, imageWidth, 2);
	CoordPoints c8 = { k_max, i_max, j_max };
	c8 = transformPoint(&c8, translation, scaling, rotation, centroid, imageHeight, imageLength, imageWidth, 2);
	//==============================================================================
	size_t min_ab, min_cd, min_ef, min_gh, min_1, min_2;
	// K min
	min_ab = c1.k < c2.k ? c1.k : c2.k; // c1 vs c2
	min_cd = c3.k < c4.k ? c3.k : c4.k; //  c3 vs c4
	min_1 = min_ab < min_cd ? min_ab : min_cd; // min. of cube floor point k

	min_ef = c4.k < c5.k ? c4.k : c5.k; // c4 vs c5
	min_gh = c7.k < c8.k ? c7.k : c8.k; //  c7 vs c8
	min_2 = min_ef < min_gh ? min_ef : min_gh; // min. of cube ceiling point k

	(*bestfit).k_min = min_1 < min_2 ? min_1 : min_2; // min. k in all 8 corners
	// I min
	min_ab = c1.i < c2.i ? c1.i : c2.i; // c1 vs c2
	min_cd = c3.i < c4.i ? c3.i : c4.i; //  c3 vs c4
	min_1 = min_ab < min_cd ? min_ab : min_cd; // min. of cube floor point i

	min_ef = c4.i < c5.i ? c4.i : c5.i; // c4 vs c5
	min_gh = c7.i < c8.i ? c7.i : c8.i; //  c7 vs c8
	min_2 = min_ef < min_gh ? min_ef : min_gh; // min. of cube ceiling point i

	(*bestfit).i_min = min_1 < min_2 ? min_1 : min_2; // min. i in all 8 corners
	// J min
	min_ab = c1.j < c2.j ? c1.j : c2.j; // c1 vs c2
	min_cd = c3.j < c4.j ? c3.j : c4.j; //  c3 vs c4
	min_1 = min_ab < min_cd ? min_ab : min_cd; // min. of cube floor point j

	min_ef = c4.j < c5.j ? c4.j : c5.j; // c4 vs c5
	min_gh = c7.j < c8.j ? c7.j : c8.j; //  c7 vs c8
	min_2 = min_ef < min_gh ? min_ef : min_gh; // min. of cube ceiling point j

	(*bestfit).j_min = min_1 < min_2 ? min_1 : min_2; // min. j in all 8 corners
	//==============================================================================
	size_t max_ab, max_cd, max_ef, max_gh, max_1, max_2;
	// Max End Points
	// K max
	max_ab = c1.k > c2.k ? c1.k : c2.k; // c1 vs c2
	max_cd = c3.k > c4.k ? c3.k : c4.k; //  c3 vs c4
	max_1 = max_ab > max_cd ? max_ab : max_cd; // max. of cube floor point k

	max_ef = c5.k > c6.k ? c5.k : c6.k; // c4 vs c5
	max_gh = c7.k > c8.k ? c7.k : c8.k; //  c7 vs c8
	max_2 = max_ef > max_gh ? max_ef : max_gh; // max. of cube ceiling point k

	(*bestfit).k_max = max_1 > max_2 ? max_1 : max_2; // max k in all 8 corners
	// I max
	max_ab = c1.i > c2.i ? c1.i : c2.i; // c1 vs c2
	max_cd = c3.i > c4.i ? c3.i : c4.i; //  c3 vs c4
	max_1 = max_ab > max_cd ? max_ab : max_cd; // max. of cube floor point i

	max_ef = c5.i > c6.i ? c5.i : c6.i; // c4 vs c5
	max_gh = c7.i > c8.i ? c7.i : c8.i; //  c7 vs c8
	max_2 = max_ef > max_gh ? max_ef : max_gh; // max. of cube ceiling point i

	(*bestfit).i_max = max_1 > max_2 ? max_1 : max_2; // max i in all 8 corners
	// J max
	max_ab = c1.j > c2.j ? c1.j : c2.j; // c1 vs c2
	max_cd = c3.j > c4.j ? c3.j : c4.j; //  c3 vs c4
	max_1 = max_ab > max_cd ? max_ab : max_cd; // max. of cube floor point k

	max_ef = c5.j > c6.j ? c5.j : c6.j; // c4 vs c5
	max_gh = c7.j > c8.j ? c7.j : c8.j; //  c7 vs c8
	max_2 = max_ef > max_gh ? max_ef : max_gh; // max. of cube ceiling point j

	(*bestfit).j_max = max_1 > max_2 ? max_1 : max_2; // max j in all 8 corners
	//==============================================================================
}
//==============================================================================
CoordPoints transformPoint(CoordPoints * inputPoints, Point3D translation, Point3D scaling, Point3D rotation, dataType centroid[3], size_t imageHeight, size_t imageLength, size_t imageWidth, int loc)
{
	//==============================================================================
	size_t bottom, left, begin;
	//==============================================================================
	dataType k_a, i_a, j_a; // Affine indices
	// Transformed
	dataType k_t, i_t, j_t; // Transformed indices
	// Temporary parameters
	dataType tmpX, tmpY, tmpZ;
	dataType cz = centroid[2], cx = centroid[0], cy = centroid[1];
	dataType theta = (rotation.y), psi = (rotation.z), phi = (rotation.x);
	//==============================================================================
	size_t k = (*inputPoints).k, i = (*inputPoints).i, j = (*inputPoints).j;
	//==============================================================================
	// 1. Move to origin for Min
	k_a = k - cz; // Move to origin Z
	i_a = i - cx; // Move to origin x
	j_a = j - cy; // Move to origin Y
	//==============================================================================
	// Apply scaling
	tmpZ = k_a / scaling.z;
	tmpX = i_a / scaling.x;
	tmpY = j_a / scaling.y;
	//==============================================================================
	// Apply Rotation
	i_t = x_rotate(tmpZ, tmpX, tmpY, theta, psi);
	j_t = y_rotate(tmpZ, tmpX, tmpY, theta, psi, phi);
	k_t = z_rotate(tmpZ, tmpX, tmpY, theta, psi, phi);
	//==============================================================================
	// Move back to centroid
	tmpX = i_t + cx;
	tmpY = j_t + cy;
	tmpZ = k_t + cz;
	//==============================================================================
	// Set the values
	i_t = tmpX;
	j_t = tmpY;
	k_t = tmpZ;
	//==============================================================================
	// Add translation - Translation already included in the points!
	i_t = i_t - translation.x;
	j_t = j_t - translation.y;
	k_t = k_t - translation.z;

	if (loc == 1) // Lower
	{
		bottom = (size_t)floor(k_t);
		bottom = max(bottom, 0);
		bottom = min(bottom, imageHeight - 1);
		(*inputPoints).k = bottom;
		// X
		size_t left = (size_t)floor(i_t);
		left = max(left, 0);
		left = min(left, imageLength - 1);
		(*inputPoints).i = left;
		// Y
		size_t begin = (size_t)floor(j_t);
		begin = max(begin, 0);
		begin = min(begin, imageWidth - 1);
		(*inputPoints).j = begin;
	}
	else if (loc == 2) // Upper
	{
		bottom = (int)ceil(k_t);
		bottom = max(bottom, 0);
		bottom = min(bottom, imageHeight - 1);
		(*inputPoints).k = bottom;
		// X
		left = (int)ceil(i_t);
		left = max(left, 0);
		left = min(left, imageLength - 1);
		(*inputPoints).i = left;
		// Y
		begin = (int)ceil(j_t);
		begin = max(begin, 0);
		begin = min(begin, imageWidth - 1);
		(*inputPoints).j = begin;
	}
	return *inputPoints;
}
//==============================================================================
dataType getDistance(dataType ** binaryImage, size_t imageHeight, size_t imageLength, size_t dim2D, const size_t k1, const size_t x1, ClipBox bestfitBox, Point3D * surface_points, size_t ptsNum, dataType insideShapevalue, bool parallelize)
{
	dataType pv = 10000;
	//==============================================================================
	dataType dx, dy, dz, dist;
	//==============================================================================
	dataType tmpX, tmpY, tmpZ;
	//==============================================================================
	tmpZ = (dataType)k1;
	//==============================================================================
	tmpX = (dataType)(x1 / imageLength);
	//==============================================================================
	tmpY = (dataType)(x1 % imageLength);
	//==============================================================================
	Point3D tmpXYZ;
	tmpXYZ.z = (dataType)k1;
	//==============================================================================
	tmpXYZ.x = (dataType)(x1 / imageLength);
	//==============================================================================
	tmpXYZ.y = (dataType)(x1 % imageLength);

	//==============================================================================
	int i;
	//==============================================================================
	if (!parallelize) // Run sequential code
	{
		//==============================================================================
		// Sequetial
		for (i = 0; i < ptsNum; i++)
		{
			//==============================================================================
			dz = tmpZ - surface_points[i].z; // difference between z axes of both images
			dx = tmpX - surface_points[i].x;// difference between x axes of both images
			dy = tmpY - surface_points[i].y;// difference between y axes of both images
			dist = dx * dx + dy * dy + dz * dz;
			//==============================================================================
			if (dist <= pv)
			{
				pv = dist;
			}
			//==============================================================================
		}
		//==============================================================================
	}
	else // Run parallelized
	{
		//==============================================================================
		// OpenMp
		int nthreads;
		dataType distances[NUM_THREADS][PAD];
		omp_set_dynamic(0); // Disable dynamic adjustment of threads
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel
		{
			int i, id, nthrds;
			dataType dx, dy, dz, dist, distPv;

			id = omp_get_thread_num();
			nthrds = omp_get_num_threads();

			if (id == 0) { nthreads = nthrds; }
			for (i = id, distPv = 10000; i < ptsNum; i = i + nthrds)
			{
				//==============================================================================
				//printf("Running on thread %d\n", id);
				//==============================================================================
				dz = tmpZ - surface_points[i].z; // difference between z axes of both images
				dx = tmpX - surface_points[i].x;// difference between x axes of both images
				dy = tmpY - surface_points[i].y;// difference between y axes of both images
				dist = dx * dx + dy * dy + dz * dz;
				//==============================================================================
				//distPv = min(distPv, dist);
				if (dist <= distPv)
				{
					distPv = dist;
				}
				//==============================================================================
				//distances[id][0] = distPv;
				//==============================================================================
			}
			distances[id][0] = distPv;
			//==============================================================================
		}
		for (size_t i = 0; i < nthreads; i++)
		{
			if (pv >= distances[i][0])
			{
				pv = distances[i][0];
			}
		}
	}
	//==============================================================================
	pv = (dataType)sqrt(pv);
	//==============================================================================
	// Check if the point is inside
	if (binaryImage[k1][x1] == insideShapevalue)
	{
		// pv = -1 * pv;
		pv = 0;
	}
	//==============================================================================
	return pv;
	//==============================================================================
}
//==============================================================================
// Returns points for the surface/object
size_t surfacePoints(dataType ** binaryImage, size_t imageLength, const unsigned char fgroundValue, ClipBox bestfitBox)
{
	//==============================================================================
	size_t k_2, i_2;//loop counter for z dimension
	//==============================================================================
	size_t ptsNum = 0;
	//==============================================================================
	// In the clip box volume
	size_t k_min = bestfitBox.k_min, k_max = bestfitBox.k_max + 1, i_min = bestfitBox.i_min, i_max = bestfitBox.i_max, j_min = bestfitBox.j_min, j_max = bestfitBox.j_max;
	size_t i_2_max = x_new(i_max, j_max, imageLength), i_2_min = x_new(i_min, j_min, imageLength);
	// Find the surcae points only within the clipbox
	for (k_2 = k_min; k_2 < k_max; k_2++) // z axis of the input surface or image
	{
		for (i_2 = i_2_min; i_2 < i_2_max; i_2++)// x-y axis of the input surface or image
		{
			if (binaryImage[k_2][i_2] == fgroundValue)
			{
				ptsNum++;
			}
		}
	}
	return ptsNum;
}
//==============================================================================
// Converts unsigned to signed distance map
void converToSignedDist(dataType ** distFn, dataType ** originalData, size_t imageHeight, size_t imageLength, size_t imageWidth, const unsigned char fgroundValue)
{
	for (size_t k = 0; k < imageHeight; k++)
	{
		for (size_t i = 0; i < imageLength; i++)
		{
			for (size_t j = 0; j < imageWidth; j++)
			{
				size_t xd = x_new(i, j, imageLength);
				if (originalData[k][xd] == fgroundValue)
				{
					distFn[k][xd] = -1 * distFn[k][xd];
				}
			}
		}
	}
}
//==============================================================================
dataType errorCalc(dataType ** aPtr, dataType ** bPtr, size_t height, size_t length, size_t width, dataType h_val)
{
	dataType error = 0.0, tmp;
	size_t i, j, k, xd, count = 0;
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				// 1D Conversion of row and column
				xd = x_new(i, j, length);
				// Error calculation
				tmp = (aPtr[k][xd] - bPtr[k][xd]) * h_val;
				error += (tmp*tmp);
				count++;
			}
		}
	}
	// Mean Square Error
	if (count == 0)
	{
		count = 1;
	}
	error = (dataType) ((error) / (2. * count));
	//return sqrtf(error); // // Root Mean Square Error
	return error; // // Mean Square Error
}
//==============================================================================