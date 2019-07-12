//==============================================================================
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
//==============================================================================
#include "shape_registration.h"
//==============================================================================
#include <omp.h>
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
dataType getDistance(dataType ** binaryImage, size_t imageHeight, size_t imageLength, size_t dim2D, const size_t k1, const size_t x1, const unsigned char fgroundValue, ClipBox bestfitBox, Point3D * surface_points, size_t ptsNum, dataType insideShapevalue, bool parallelize);
size_t surfacePoints(dataType ** binaryImage, size_t imageLength, const unsigned char fgroundValue, ClipBox bestfitBox);
//==============================================================================
void nbPointsX(dataType ** transformedBinaryData, dataType pixelSize, size_t x, size_t k, size_t i, size_t imageHeight, size_t imageLength, size_t imageWidth, T * pValue, T * qValue, T *hh, ClipBox bestfitBox, Point3D * surface_points, size_t ptsNum);
void nbPointsY(dataType ** transformedBinaryData, dataType pixelSize, size_t k, size_t i, size_t j, size_t imageHeight, size_t imageLength, size_t imageWidth, T * pValue, T * qValue, T *hh, ClipBox bestfitBox, Point3D * surface_points, size_t ptsNum);
void nbPointsZ(dataType ** transformedBinaryData, dataType pixelSize, size_t x, size_t k, size_t imageHeight, size_t imageLength, size_t imageWidth, T * pValue, T * qValue, T *hh, ClipBox bestfitBox, Point3D * surface_points, size_t ptsNum);
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
	transform3DImage(movingData, resultPtr, finalResults.translation, finalResults.scaling, finalResults.rotation, zDim, xDim, yDim, params.imageBackground, movingCentroid);
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
	Obj_Structure ** objectNthD = (Obj_Structure **)malloc(sizeof(Obj_Structure*)*imageHeight);
	for (i = 0; i < imageHeight; i++)
	{
		objectNthD[i] = (Obj_Structure *)malloc(sizeof(Obj_Structure) * (imageLength*imageWidth));
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
	Point3D *shapePoints = (Point3D *)malloc(sizeof(Point3D)*(imageHeight*imageLength*imageWidth));
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
					// Save the dimension with those values
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
inline int NFunctionBinary(dataType v1, dataType v2, dataType delta)
{
	//==============================================================================
	if (abs(v1) == delta || abs(v2) == delta)
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
	if (abs(v1) > delta)
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
				int x = x_new(i, j, imageLength);
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
				int x = x_new(i, j, imageLength);
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
dataType finiteDifAll(dataType pValue, dataType qValue, dataType hh)
{
	return (pValue - qValue) / hh;
}
//==============================================================================
Affine_Parameter gradientComponents(dataType ** destPtr, dataType ** distTrans, dataType h, Affine_Parameter * params, size_t imageHeight, size_t imageLength, size_t imageWidth)
{
	Affine_Parameter results;
	// Initialize the parameters
	size_t k, i, j, x;

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

	// Setting same if uniform scale and rotation
#ifdef UNIFORM
	sx = sy = sz;
	phi = theta = psi;
#endif // UNIFORM
	// Scale denominator
	dataType inv_sx2 = 1.0 / (sx*sx);
	dataType inv_sy2 = 1.0 / (sy*sy);
	dataType inv_sz2 = 1.0 / (sz*sz);

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
					// Store the distance function difference
					distDifference = (destPtr[k][x] - distTrans[k][x]) * 2.0;

					// Directional component vector derivatives - i, j, k
					dataType tmpI = i / (dataType)imageLength, tmpJ = j / (dataType)imageWidth, tmpK = k / (dataType)imageHeight;

					// Trignometry functions inside the component evaluation equation
					dataType a = cos(phi), b = sin(phi), ab = cos(phi)*sin(phi), aa = cos(phi)*cos(phi), bb = sin(phi)*sin(phi), aaa = a * aa, bbb = b * bb, aab = aa * b, abb = a * bb;

					// Apply Forward Differences to the distTrans pointer
					xFwd = finiteDifX(distTrans, h, x, k, i, imageLength);
					yFwd = finiteDifY(distTrans, h, k, i, j, imageLength, imageWidth);
					zFwd = finiteDifZ(distTrans, h, x, k, i, imageLength, imageHeight);
					// Evaluate Individual Gradient Components
#ifdef UNIFORM
					// Uniform Rotations Angles
					component = xFwd * ((-2.0*tmpI)*(ab / sx) + ((tmpJ)*((bb - aa) / sx)) + ((tmpK)*(a / sx))) +
						yFwd * (((tmpI)*((aa + 2.0 * aab - bb - bbb) / sy)) + ((tmpJ)*((-2.0 * ab - 3.0 * abb) / sy)) + ((tmpK)*((bb - aa) / sy))) +
						zFwd * (((tmpI)*((-aaa + 2.0 * ab + 2.0 * abb) / sz)) + ((tmpJ)*((aa + 2.0 * aa - bb - bbb) / sz)) + ((tmpK)*((-2.0 * ab) / sz)));
					// Set the Rotation - Uniform Angles
					results.rotation.x += (component)*(distDifference);
					results.rotation.y += (component)*(distDifference);
					results.rotation.z += (component)*(distDifference);
					// Scaling Parameters
					Uniform Scale Component
						component = xFwd * (((tmpI)*((-aa) * inv_sx2)) + ((tmpJ)*(ab * inv_sx2)) + ((tmpK)*((-b) * inv_sx2))) +
						yFwd * (((-tmpI)*((ab + abb)  * inv_sy2)) + ((-tmpJ)*((aa - bbb) * inv_sy2)) + ((tmpK)*(ab * inv_sy2))) +
						zFwd * (((-tmpI)*((-aab + bb) * inv_sz2)) + ((-tmpJ)*((ab + abb) * inv_sz2)) + ((-tmpK)*(aa * inv_sz2)));
					// Set the Scales - Uniform Scale
					results.scaling.x += (component)*(distDifference);
					results.scaling.y += (component)*(distDifference);
					results.scaling.z += (component)*(distDifference);
#endif // UNIFORM
#ifdef DIRECTIONAL
					// Rotation Components - Directionnal
					componentX = yFwd * (((tmpI)*((-sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)) / sy)) + ((tmpJ)*((-cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta)) / sy)) + ((-tmpK)*((cos(phi)*cos(theta)) / sy))) +
						zFwd * (((tmpI)*((cos(phi)*sin(psi) + cos(psi)*sin(phi)*sin(theta)) / sz)) + ((tmpJ)*((cos(phi)*cos(psi) - sin(phi)*sin(psi)*sin(theta)) / sz)) + ((-tmpK)*((cos(theta) * sin(phi)) / sz)));
					componentY = xFwd * (((-tmpI)*((cos(psi)*sin(theta)) / sx)) + ((tmpJ)*((sin(psi)*sin(theta)) / sx)) + ((tmpK)*((cos(theta)) / sx))) +
						yFwd * (((tmpI)*((cos(psi)*cos(theta)*sin(phi)) / sy)) + ((-tmpJ)*((cos(theta)*sin(phi)*sin(psi)) / sy)) + ((tmpK)*((sin(phi)*sin(theta)) / sy))) +
						zFwd * (((-tmpI)*((cos(phi)*cos(psi)*cos(theta)) / sz)) + ((tmpJ)*((cos(phi)*cos(theta)*sin(psi)) / sz)) + ((-tmpK)*((cos(phi)*sin(theta)) / sz)));
					componentZ = xFwd * (((-tmpI)*((cos(theta)*sin(psi)) / sx)) + ((-tmpJ)*((cos(psi)*cos(theta)) / sx))) +
						yFwd * (((tmpI)*((cos(phi)*cos(psi) - sin(phi)*sin(psi)*sin(theta)) / sy)) + ((tmpJ)*((-cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)) / sy))) +
						zFwd * (((tmpI)*((cos(psi)*sin(phi) + cos(phi)*sin(psi)*sin(theta)) / sz)) + ((tmpJ)*((-sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)) / sz)));
					// Set the Rotations - Directional
					results.rotation.x += (componentX)*(distDifference);
					results.rotation.y += (componentY)*(distDifference);
					results.rotation.z += (componentZ)*(distDifference);
					// Directional Scale Components
					componentX = xFwd * ((-tmpI)*((cos(psi)*cos(theta)) / (sx*sx)) + ((tmpJ)*((cos(theta)*sin(psi)) / (sx*sx))) + ((-tmpK)*((sin(theta)) / (sx*sx))));
					componentY = yFwd * (((-tmpI)*((cos(phi)*sin(psi) + cos(psi)*sin(phi)*sin(theta)) / (sy*sy))) + ((-tmpJ)*((cos(phi)*cos(psi) - sin(phi)*sin(psi)*sin(theta)) / (sy*sy))) + ((tmpK)*((cos(theta)*sin(phi)) / (sy*sy))));
					componentZ = zFwd * (((-tmpI)*((sin(phi)*sin(psi) - cos(phi)*cos(psi)*sin(theta)) / (sz*sz))) + ((-tmpJ)*((cos(psi)*sin(phi) + cos(phi)*sin(psi)*sin(theta)) / (sz*sz))) + ((-tmpK)*((cos(phi)*cos(theta)) / (sz*sz))));
					// Set the Scales - Directional Scales
					results.scaling.x += (componentX)*(distDifference);
					results.scaling.y += (componentY)*(distDifference);
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
	return results;
}
//==============================================================================
Affine_Parameter gradientComponentsClip(dataType ** destPtr, dataType ** distTrans, dataType h, Affine_Parameter * params, size_t imageHeight, size_t imageLength, size_t imageWidth, ClipBox bestFit)
{
	Affine_Parameter results;
	// Initialize the parameters
	int k, i, j, x, counter = 0;

	// Initialize the results
	results.rotation.x = 0.0, results.rotation.y = 0.0, results.rotation.z = 0.0;
	results.scaling.x = 0.0, results.scaling.y = 0.0, results.scaling.z = 0.0;
	results.translation.x = 0.0, results.translation.y = 0.0, results.translation.z = 0.0;

	// Forward difference parameters
	dataType xFwd, yFwd, zFwd;

	// Derivative component
	dataType component, componentX, componentY, componentZ;

	// Stores the difference between two distance pointers
	dataType distDifference;

	// Shorter Transformation names
	dataType phi = params->rotation.x, theta = params->rotation.y, psi = params->rotation.z;
	dataType sx = params->scaling.x, sy = params->scaling.y, sz = params->scaling.z;
	dataType tx = params->translation.x, ty = params->translation.y, tz = params->translation.z;

	// Setting same if uniform scale and rotation
#ifdef UNIFORM
	sx = sy = sz;
	phi = theta = psi;
#endif // UNIFORM
	// Scale denominator
	dataType inv_sx2 = 1.0 / (sx*sx);
	dataType inv_sy2 = 1.0 / (sy*sy);
	dataType inv_sz2 = 1.0 / (sz*sz);

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
					distDifference = (destPtr[k][x] - distTrans[k][x]) * 2.0;

					// Directional component vector derivatives - i, j, k
					dataType tmpI = i / (dataType)imageLength, tmpJ = j / (dataType)imageWidth, tmpK = k / (dataType)imageHeight;

					// Trignometry functions inside the component evaluation equation
					dataType a = cos(phi), b = sin(phi), ab = cos(phi)*sin(phi), aa = cos(phi)*cos(phi), bb = sin(phi)*sin(phi), aaa = a * aa, bbb = b * bb, aab = aa * b, abb = a * bb;

					// Apply Forward Differences to the distTrans pointer
					xFwd = finiteDifX(distTrans, h, x, k, i, imageLength);
					yFwd = finiteDifY(distTrans, h, k, i, j, imageLength, imageWidth);
					zFwd = finiteDifZ(distTrans, h, x, k, i, imageLength, imageHeight);
					// Evaluate Individual Gradient Components
#ifdef UNIFORM
					// Uniform Rotations Angles
					component = xFwd * ((-2.0*tmpI)*(ab / sx) + ((tmpJ)*((bb - aa) / sx)) + ((tmpK)*(a / sx))) +
						yFwd * (((tmpI)*((aa + 2.0 * aab - bb - bbb) / sy)) + ((tmpJ)*((-2.0 * ab - 3.0 * abb) / sy)) + ((tmpK)*((bb - aa) / sy))) +
						zFwd * (((tmpI)*((-aaa + 2.0 * ab + 2.0 * abb) / sz)) + ((tmpJ)*((aa + 2.0 * aa - bb - bbb) / sz)) + ((tmpK)*((-2.0 * ab) / sz)));
					// Set the Rotation - Uniform Angles
					results.rotation.x += (component)*(distDifference);
					results.rotation.y += (component)*(distDifference);
					results.rotation.z += (component)*(distDifference);
					// Scaling Parameters
					Uniform Scale Component
						component = xFwd * (((tmpI)*((-aa) * inv_sx2)) + ((tmpJ)*(ab * inv_sx2)) + ((tmpK)*((-b) * inv_sx2))) +
						yFwd * (((-tmpI)*((ab + abb)  * inv_sy2)) + ((-tmpJ)*((aa - bbb) * inv_sy2)) + ((tmpK)*(ab * inv_sy2))) +
						zFwd * (((-tmpI)*((-aab + bb) * inv_sz2)) + ((-tmpJ)*((ab + abb) * inv_sz2)) + ((-tmpK)*(aa * inv_sz2)));
					// Set the Scales - Uniform Scale
					results.scaling.x += (component)*(distDifference);
					results.scaling.y += (component)*(distDifference);
					results.scaling.z += (component)*(distDifference);
#endif // UNIFORM
#ifdef DIRECTIONAL
					// Rotation Components - Directionnal
					componentX = yFwd * (((tmpI)*((-sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)) / sy)) + ((tmpJ)*((-cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta)) / sy)) + ((-tmpK)*((cos(phi)*cos(theta)) / sy))) +
						zFwd * (((tmpI)*((cos(phi)*sin(psi) + cos(psi)*sin(phi)*sin(theta)) / sz)) + ((tmpJ)*((cos(phi)*cos(psi) - sin(phi)*sin(psi)*sin(theta)) / sz)) + ((-tmpK)*((cos(theta) * sin(phi)) / sz)));
					componentY = xFwd * (((-tmpI)*((cos(psi)*sin(theta)) / sx)) + ((tmpJ)*((sin(psi)*sin(theta)) / sx)) + ((tmpK)*((cos(theta)) / sx))) +
						yFwd * (((tmpI)*((cos(psi)*cos(theta)*sin(phi)) / sy)) + ((-tmpJ)*((cos(theta)*sin(phi)*sin(psi)) / sy)) + ((tmpK)*((sin(phi)*sin(theta)) / sy))) +
						zFwd * (((-tmpI)*((cos(phi)*cos(psi)*cos(theta)) / sz)) + ((tmpJ)*((cos(phi)*cos(theta)*sin(psi)) / sz)) + ((-tmpK)*((cos(phi)*sin(theta)) / sz)));
					componentZ = xFwd * (((-tmpI)*((cos(theta)*sin(psi)) / sx)) + ((-tmpJ)*((cos(psi)*cos(theta)) / sx))) +
						yFwd * (((tmpI)*((cos(phi)*cos(psi) - sin(phi)*sin(psi)*sin(theta)) / sy)) + ((tmpJ)*((-cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)) / sy))) +
						zFwd * (((tmpI)*((cos(psi)*sin(phi) + cos(phi)*sin(psi)*sin(theta)) / sz)) + ((tmpJ)*((-sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)) / sz)));
					// Set the Rotations - Directional
					results.rotation.x += (componentX)*(distDifference);
					results.rotation.y += (componentY)*(distDifference);
					results.rotation.z += (componentZ)*(distDifference);
					// Directional Scale Components
					componentX = xFwd * ((-tmpI)*((cos(psi)*cos(theta)) / (sx*sx)) + ((tmpJ)*((cos(theta)*sin(psi)) / (sx*sx))) + ((-tmpK)*((sin(theta)) / (sx*sx))));
					componentY = yFwd * (((-tmpI)*((cos(phi)*sin(psi) + cos(psi)*sin(phi)*sin(theta)) / (sy*sy))) + ((-tmpJ)*((cos(phi)*cos(psi) - sin(phi)*sin(psi)*sin(theta)) / (sy*sy))) + ((tmpK)*((cos(theta)*sin(phi)) / (sy*sy))));
					componentZ = zFwd * (((-tmpI)*((sin(phi)*sin(psi) - cos(phi)*cos(psi)*sin(theta)) / (sz*sz))) + ((-tmpJ)*((cos(psi)*sin(phi) + cos(phi)*sin(psi)*sin(theta)) / (sz*sz))) + ((-tmpK)*((cos(phi)*cos(theta)) / (sz*sz))));
					// Set the Scales - Directional Scales
					results.scaling.x += (componentX)*(distDifference);
					results.scaling.y += (componentY)*(distDifference);
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
#ifdef USE_CLIP
	dataType ** movInitPtr = (dataType **)malloc(sizeof(dataType *) * imageHeight); // distances for Moving
	for (i = 0; i < imageHeight; i++)
	{
		movInitPtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
	}
#endif // USE_CLIP
	//==============================================================================
	// Initializations of Pointers
	for (i = 0; i < imageHeight; i++)
	{
		destPtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
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
	fastSweepingFunction_3D(destPtr, fixedData, imageLength, imageWidth, imageHeight, params.h, 50000, params.imageForeground);
	//==============================================================================
#ifdef USE_CLIP
	// Initial dist. fn for moving image before adding any transformation
	fastSweepingFunction_3D(movInitPtr, movingData, imageLength, imageWidth, imageHeight, params.h, 50000, params.imageForeground);
#endif // USE_CLIP
	//==============================================================================
#ifdef MEASURE_TIME
	secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
	// Store the time
	distanceTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
	//==============================================================================
#ifdef USE_CLIP
	//==============================================================================
	// Finding the clip box points for the fixed image
	ClipBox coordFixed = findClipBoxSingle(destPtr, imageHeight, imageLength, imageWidth);
	//==============================================================================
	ClipBox coordMoving = findClipBoxSingle(movInitPtr, imageHeight, imageLength, imageWidth);
	//==============================================================================
	// Free after
	for (k = 0; k < imageHeight; k++)
	{
		free(movInitPtr[k]);
	}
	free(movInitPtr);
	movInitPtr = NULL;
	//==============================================================================
	ClipBox bestFit, coordMovingTmp; // Clipbox for bestFit of both fixed and moving images, Moving image clipbox
	//==============================================================================
#endif // USE_CLIP
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
		transPtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
		distTransPtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
	}
	//==============================================================================
	while (!stopCond)
	{
		//==============================================================================
		// Timing The Transformation Inside the registration function
#ifdef MEASURE_TIME
		firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
		transform3DImage(movingData, transPtr, affineResult.translation, affineResult.scaling, affineResult.rotation, imageHeight, imageLength, imageWidth, params.imageBackground, centroid);
		//==============================================================================
#ifdef MEASURE_TIME
		secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
		// Store the time
		transformationTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
		//==============================================================================
#ifdef USE_CLIP
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
#endif // USE_CLIP
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
#ifdef USE_CLIP
		fSweeping3D(distTransPtr, transPtr, imageLength, imageWidth, imageHeight, 1, 50000, params.imageForeground, bestFit);
#else
		fastSweepingFunction_3D(distTransPtr, transPtr, imageLength, imageWidth, imageHeight, params.h, 50000, params.imageForeground);
#endif // USE_CLIP		
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
#ifdef USE_CLIP
		// Evaluate Energy Function - L2 Norm Between the two calc. distances within the band and clipbox
		energyTmp = energyFunctionClip(destPtr, distTransPtr, bestFit, imageLength);
#else
		// Evaluate Energy Function - L2 Norm Between the two calc. distances
		energyTmp = energyFunction(destPtr, distTransPtr, imageHeight, imageLength, imageWidth, params.h);
#endif // USE_CLIP		
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
#ifdef USE_CLIP
			affineTmp = gradientComponentsClip(destPtr, distTransPtr, params.h, &affineResult, imageHeight, imageLength, imageWidth, bestFit);

#else
			affineTmp = gradientComponents(destPtr, distTransPtr, params.h, &affineResult, imageHeight, imageLength, imageWidth);
#endif // USE_CLIP			
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
			// Update the weights
			// updateWeights(&rotation_weight, &scaling_weight, &translation_weight, affineResult, energyTmp, lambda);
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
#ifdef UPDATEWEIGHT
			// Update Rotation weight
			updateWeight(&rotation_weight, initTransform.rotation, affineResult.rotation, lambda);
			// Update Scale weight
			updateWeight(&scaling_weight, initTransform.scaling, affineResult.scaling, lambda);
			// Update Translation weight
			updateWeight(&translation_weight, initTransform.translation, affineResult.translation, lambda);
#endif // UPDATEWEIGHT
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
	size_t k, i, j, l, x, dim2D = imageLength * imageWidth, iteration = 0, maxSurfacePts = 0.05 * dim2D * imageHeight;
	const dataType h = 1.0;
	dataType firstCpuTime, secondCpuTime, regStartCpuTime, regStopCpuTime, regTotalCpuTimen = 0.;
	dataType energyTotalCpuTime = 0., distanceTotalCpuTime = 0., gradientTotalCpuTime = 0., transformationTotalCpuTime = 0., conversionTotalCpuTime = 0., surfacePtsTotalCpuTime = 0., edgeDetectionTotalCpuTime = 0., generateRandomPtsTotalCpuTime = 0., distanceCalculateTotalCpuTime = 0., finiteDiffereceTotalCpuTime = 0.;
	//==============================================================================
	Point3D * surface_points = malloc(sizeof(Point3D) * maxSurfacePts);
	//==============================================================================
	// Affine Parameters
	Affine_Parameter affineResult;
	//==============================================================================
	// Create new fixed dist. Pointers to be used
	dataType ** destPtr = (dataType **)malloc(sizeof(dataType *) * imageHeight); // distances for destination
	//==============================================================================
#ifdef USE_CLIP
	dataType ** movInitPtr = (dataType **)malloc(sizeof(dataType *) * imageHeight); // distances for Moving
	for (i = 0; i < imageHeight; i++)
	{
		movInitPtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
	}
#endif // USE_CLIP
	//==============================================================================
	// Initializations of Pointers
	for (i = 0; i < imageHeight; i++)
	{
		destPtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
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
	for (i = 0; i < imageHeight; i++)
	{
		edgeMovingPointer[i] = malloc(sizeof(dataType) * dim2D);
	}
	//==============================================================================
	// Apply distance function between transPtr and distTrans
	// Begin Record Time
#ifdef MEASURE_TIME
	firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
	//==============================================================================
	fastSweepingFunction_3D(destPtr, fixedData, imageLength, imageWidth, imageHeight, params.h, 50000, params.imageForeground);
	//==============================================================================
#ifdef USE_CLIP
	//==============================================================================
	fastSweepingFunction_3D(movInitPtr, movingData, imageLength, imageWidth, imageHeight, params.h, 50000, params.imageForeground);
	//==============================================================================
#endif // USE_CLIP
	//==============================================================================
#ifdef MEASURE_TIME
	secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
	// Store the time
	distanceTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
	//==============================================================================
#ifdef USE_CLIP
	//==============================================================================
	// Finding the clip box points for the fixed image dist. unsigned
	ClipBox coordFixed = findClipBoxSingle(destPtr, imageHeight, imageLength, imageWidth);
	//==============================================================================
	// From unsigned dist. fn
	ClipBox coordMoving = findClipBoxSingle(movInitPtr, imageHeight, imageLength, imageWidth);
	//==============================================================================
	// Calc. the narrow band areas for fixed and moving dist. fn's
	dataType ** fixedNBandPtr = (dataType **)malloc(sizeof(dataType *) * imageHeight); // narrow band area for fixed dista. fn
	dataType ** movingNBandPtr = (dataType **)malloc(sizeof(dataType *) * imageHeight); // narrow band area for moving dista. fn
	for (i = 0; i < imageHeight; i++)
	{
		fixedNBandPtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
		movingNBandPtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
	}
	// Fill the narrow band areas for fixed, moving respectively
	fillNarrowBandArea(destPtr, fixedNBandPtr, imageHeight, imageLength, imageWidth, params.imageForeground, params.imageBackground);
	fillNarrowBandArea(movInitPtr, movingNBandPtr, imageHeight, imageLength, imageWidth, params.imageForeground, params.imageBackground);
	// Centroid for moving narrow band area
	dataType centroidMovingBandArea[3];
	centroidImage(movingNBandPtr, centroidMovingBandArea, imageHeight, imageLength, imageWidth, params.imageBackground);
	//==============================================================================
	// Free after
	for (k = 0; k < imageHeight; k++)
	{
		free(movInitPtr[k]);
	}
	free(movInitPtr);
	movInitPtr = NULL;
	//==============================================================================
	ClipBox bestFit, coordMovingTmp; // Clipbox for bestFit of both fixed and moving images, Moving image clipbox
	//==============================================================================
#endif // USE_CLIP
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
		transPtr[i] = malloc(sizeof(dataType) * dim2D);
		transMovingPtr[i] = malloc(sizeof(dataType) * dim2D);
		distTransPtr[i] = malloc(sizeof(dataType) * dim2D);
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
		transform3DImage(movingData, transPtr, affineResult.translation, affineResult.scaling, affineResult.rotation, imageHeight, imageLength, imageWidth, params.imageBackground, centroid);
		// Transform the moving narrow band area
		transform3DImage(movingNBandPtr, transMovingPtr, affineResult.translation, affineResult.scaling, affineResult.rotation, imageHeight, imageLength, imageWidth, params.imageBackground, centroidMovingBandArea);
		// Copy back to movingNBandPtr
		dataType ** tmpPtr = NULL;
		tmpPtr = movingNBandPtr;
		movingNBandPtr = transMovingPtr;
		transMovingPtr = tmpPtr;
#ifdef MEASURE_TIME
		secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
		// Store the time
		transformationTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
		//==============================================================================
#ifdef USE_CLIP
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
#endif // USE_CLIP
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
#ifdef USE_CLIP
		fSweeping3D(distTransPtr, transPtr, imageLength, imageWidth, imageHeight, params.h, 50000, params.imageForeground, bestFit);
#else
		fastSweepingFunction_3D(distTransPtr, transPtr, imageLength, imageWidth, imageHeight, params.h, 50000, params.imageForeground);
#endif // USE_CLIP
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
#ifdef USE_CLIP
		//energyTmp = energyFunctionClip(destPtr, distTransPtr, bestFit, imageLength);
		// For using narrowband areas
		energyTmp = energyFunctionClipBandArea(destPtr, distTransPtr, bestFit, imageLength, fixedNBandPtr, movingNBandPtr, params.imageForeground);
#else
		energyTmp = energyFunction(destPtr, distTransPtr, imageHeight, imageLength, imageWidth, params.h);
#endif // USE_CLIP
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
			printf("Energy = %5.5lf, iteration %4zd, Phi = %3.2lf, Theta = %3.2lf, Psi = %3.2lf, Sx = %2.2lf, Sy = %2.2lf, Sz = %2.2lf, Tx = %2.2lf, Ty = %2.2lf, Tz = %2.2lf\n",
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
				// Free the band areas
				free(fixedNBandPtr[k]);
				free(movingNBandPtr[k]);
				// free moving edge transfromed
				free(edgeMovingPointer[k]);
				//
				free(transMovingPtr[k]);
				//
				free(distTransPtr[k]);
				//
				free(transPtr[k]); // free moving transfromed
			}
			free(destPtr);
			destPtr = NULL;
			// Band areas
			free(fixedNBandPtr);
			free(movingNBandPtr);
			fixedNBandPtr = NULL;
			movingNBandPtr = NULL;
			// free moving edge transfromed
			free(edgeMovingPointer);
			edgeMovingPointer = NULL;
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
			printf("\n\n Final Results \n\n");
			printf("Total distance Function calc. CPU Time is: %e secs\n", distanceTotalCpuTime);
			printf("Total energy Function calc. CPU Time is: %e secs\n", energyTotalCpuTime);
			printf("Total SGD Function calc. CPU Time is: %e secs\n", gradientTotalCpuTime);
			printf("Total transformation Function calc. CPU Time is: %e secs\n", transformationTotalCpuTime);
			//==============================================================================
			// Print the Calculated Transformation Parameters At the End of Registration
			printf("Energy = %5.2lf, iteration %4zd, Phi = %3.2lf, Theta = %3.2lf, Psi = %3.2lf, Sx = %2.2lf, Sy = %2.2lf, Sz = %2.2lf, Tx = %2.2lf, Ty = %2.2lf, Tz = %2.2lf\n",
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
			dataType component, componentX, componentY, componentZ;
			// Stores the difference between two distance pointers
			dataType distDifference;
			// Shorter Transformation names
			dataType phi = affineResult.rotation.x, theta = affineResult.rotation.y, psi = affineResult.rotation.z;
			dataType sx = affineResult.scaling.x, sy = affineResult.scaling.y, sz = affineResult.scaling.z;
			dataType tx = affineResult.translation.x, ty = affineResult.translation.y, tz = affineResult.translation.z;
			//==============================================================================
			// Angles to radians
			dataType _cos_phi = cos(phi), _cos_psi = cos(psi), _cos_theta = cos(theta);
			dataType _cos_phi_neg = -1 * cos(phi), _cos_psi_neg = -1 * cos(psi), _cos_theta_neg = -1 * cos(theta);

			dataType _sin_phi = sin(phi), _sin_psi = sin(psi), _sin_theta = sin(theta);
			dataType _sin_phi_neg = -1 * sin(phi), _sin_psi_neg = -1 * sin(psi), _sin_theta_neg = -1 * sin(theta);
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
			dataType inv_sx2 = 1.0 / (sx*sx);
			dataType inv_sy2 = 1.0 / (sy*sy);
			dataType inv_sz2 = 1.0 / (sz*sz);
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
			edgeDetection3dFunctionD(transPtr, edgeMovingPointer, imageLength, imageWidth, imageHeight, params.imageBackground, params.imageForeground, params.insideShapevalue);
#ifdef MEASURE_TIME
			secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
			// Store the time recorded
			edgeDetectionTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
			//==============================================================================
			dataType getDist;
			//==============================================================================
			// Point3D ptr to store points
			//==============================================================================
			size_t k_1, i_1, k_2, i_2;//loop counter for z dimension
			//==============================================================================
			size_t ptsNum = 0;
			//==============================================================================
			size_t k_min = bestFit.k_min, k_max = bestFit.k_max + 1, i_min = bestFit.i_min, i_max = bestFit.i_max, j_min = bestFit.j_min, j_max = bestFit.j_max;
			size_t i_2_max = x_new(i_max, j_max, imageLength), i_2_min = x_new(i_min, j_min, imageLength);
			//==============================================================================
			//Point3D * surface_points = malloc(sizeof(Point3D) * ptsNum);
			// Begin record time
			// Calc. the points and store them
#ifdef MEASURE_TIME
			firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
			//==============================================================================
			Point3D * tmpPt;
			int rn;
			//==============================================================================
			// Find the surcae points only within the clipbox
			for (k_2 = k_min; k_2 < k_max; k_2++) // z axis of the input surface or image
			{
				for (i_2 = i_2_min; i_2 < i_2_max; i_2++)// x-y axis of the input surface or image
				{
					//if (transPtr[k_2][i_2] == foregound)
					//==============================================================================
					// Generate random numbers 1 - 100
					rn = 1 + (rand() % 100);
					//==============================================================================
					if (edgeMovingPointer[k_2][i_2] == params.imageForeground)
					{
						//==============================================================================
						// Pick only for even random numbers
						if (rn <= 50)
						{
							tmpPt = &surface_points[ptsNum];
							//==============================================================================
							tmpPt->x = (int)(i_2 / imageLength);
							tmpPt->y = (i_2 % imageLength);
							tmpPt->z = k_2;
							ptsNum++;
						}
						//==============================================================================
						// tmpPt = &surface_points[ptsNum];
						//==============================================================================
						/*tmpPt->x = (int)(i_2 / imageLength);
						tmpPt->y = (i_2 % imageLength);
						tmpPt->z = k_2;
						ptsNum++;*/
					}
				}
			}
			//==============================================================================
			// reduce number of selected points by 50%
			//ptsNum = 0.5 * ptsNum;
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
			l = 0;
			Random3DPoints * tmpRdPts;
			do
			{
				// Generate random points
#ifdef USE_CLIP
				// From clipbox
				k = bestFit.k_min + (rand() % bestFit.k_max), i = bestFit.i_min + (rand() % bestFit.i_max), j = bestFit.j_min + (rand() % bestFit.j_max);
#else
				// From whole domain
				k = rand() % imageHeight, i = rand() % imageLength, j = rand() % imageWidth;
#endif // USE_CLIP
				// Check if inside the narrow band
				x = x_new(i, j, imageLength);
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
			} while ((loop) && (l <= params.rand_points));
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
						getDist = getDistance(edgeMovingPointer, imageHeight, imageLength, dim2D, k, x, params.imageForeground, bestFit, surface_points, ptsNum, params.insideShapevalue, params.parallelize);
						destFixed = destPtr[k][x];
						tmpDistances->distDifference = (destFixed - getDist) * 2.0;
						//==============================================================================
						nbPointsX(edgeMovingPointer, h, x, k, i, imageHeight, imageLength, imageWidth, &pVal, &qVal, &hh, bestFit, surface_points, ptsNum);
						tmpXFwd->hvl = hh;
						tmpXFwd->pvl = pVal;
						tmpXFwd->qvl = qVal;
						//==============================================================================
						nbPointsY(edgeMovingPointer, h, k, i, j, imageHeight, imageLength, imageWidth, &pVal, &qVal, &hh, bestFit, surface_points, ptsNum);
						tmpYFwd->hvl = hh;
						tmpYFwd->pvl = pVal;
						tmpYFwd->qvl = qVal;
						//==============================================================================
						nbPointsZ(edgeMovingPointer, h, x, k, imageHeight, imageLength, imageWidth, &pVal, &qVal, &hh, bestFit, surface_points, ptsNum);
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
				getDist = getDistance(edgeMovingPointer, imageHeight, imageLength, dim2D, k, x, params.imageForeground, bestFit, surface_points, ptsNum, params.insideShapevalue, params.parallelize);
				destFixed = destPtr[k][x];
				tmpDistances->distDifference = (destFixed - getDist) * 2.0;
				//==============================================================================
				nbPointsX(edgeMovingPointer, h, x, k, i, imageHeight, imageLength, imageWidth, &pVal, &qVal, &hh, bestFit, surface_points, ptsNum);
				tmpXFwd->hvl = hh;
				tmpXFwd->pvl = pVal;
				tmpXFwd->qvl = qVal;
				//==============================================================================
				nbPointsY(edgeMovingPointer, h, k, i, j, imageHeight, imageLength, imageWidth, &pVal, &qVal, &hh, bestFit, surface_points, ptsNum);
				tmpYFwd->hvl = hh;
				tmpYFwd->pvl = pVal;
				tmpYFwd->qvl = qVal;
				//==============================================================================
				nbPointsZ(edgeMovingPointer, h, x, k, imageHeight, imageLength, imageWidth, &pVal, &qVal, &hh, bestFit, surface_points, ptsNum);
				tmpZFwd->hvl = hh;
				tmpZFwd->pvl = pVal;
				tmpZFwd->qvl = qVal;
				//==============================================================================
			}
			//==============================================================================
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
				int x = x_new(i, j, imageLength);
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
				int xd = x_new(i, j, imageLength);
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
				int x = x_new(i, j, imageLength);
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
	int bottom, left, begin;
	//==============================================================================
	dataType k_a, i_a, j_a; // Affine indices
	// Transformed
	dataType k_t, i_t, j_t; // Transformed indices
	// Temporary parameters
	dataType tmpX, tmpY, tmpZ, tmp;
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
		bottom = floorf(k_t);
		bottom = max(bottom, 0);
		bottom = min(bottom, imageHeight - 1);
		(*inputPoints).k = bottom;
		// X
		int left = floorf(i_t);
		left = max(left, 0);
		left = min(left, imageLength - 1);
		(*inputPoints).i = left;
		// Y
		int begin = floorf(j_t);
		begin = max(begin, 0);
		begin = min(begin, imageWidth - 1);
		(*inputPoints).j = begin;
	}
	else if (loc == 2) // Upper
	{
		bottom = ceilf(k_t);
		bottom = max(bottom, 0);
		bottom = min(bottom, imageHeight - 1);
		(*inputPoints).k = bottom;
		// X
		left = ceilf(i_t);
		left = max(left, 0);
		left = min(left, imageLength - 1);
		(*inputPoints).i = left;
		// Y
		begin = ceilf(j_t);
		begin = max(begin, 0);
		begin = min(begin, imageWidth - 1);
		(*inputPoints).j = begin;
	}
	return *inputPoints;
}
//==============================================================================
dataType getDistance(dataType ** binaryImage, size_t imageHeight, size_t imageLength, size_t dim2D, const size_t k1, const size_t x1, const unsigned char fgroundValue, ClipBox bestfitBox, Point3D * surface_points, size_t ptsNum, dataType insideShapevalue, bool parallelize)
{
	dataType pv = 10000;
	//==============================================================================
	dataType dx, dy, dz, dist;
	//==============================================================================
	dataType tmpX, tmpY, tmpZ;
	//==============================================================================
	tmpZ = k1;
	//==============================================================================
	tmpX = (int)(x1 / imageLength);
	//==============================================================================
	tmpY = (x1 % imageLength);
	//==============================================================================
	Point3D tmpXYZ;
	tmpXYZ.z = k1;
	//==============================================================================
	tmpXYZ.x = (int)(x1 / imageLength);
	//==============================================================================
	tmpXYZ.y = (x1 % imageLength);

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
	pv = sqrt(pv);
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
	size_t k_1, i_1, k_2, i_2;//loop counter for z dimension
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