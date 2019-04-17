//==============================================================================
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
//==============================================================================
#include "shape_registration.h"
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
	// Initializations of Pointers
	for (i = 0; i < imageHeight; i++)
	{
		destPtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
	}
	//==============================================================================
	// Initialize to same background - default value is 255, 0, 0
	// initialize3dArrayD(destPtr, imageLength, imageWidth, imageHeight, 0);
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
	fastSweepingFunction_3D(destPtr, fixedData, imageLength, imageWidth, imageHeight, 1, 50000, params.imageForeground);
	//==============================================================================
#ifdef MEASURE_TIME
	secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
	// Store the time
	distanceTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
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
	while (!stopCond)
	{
		// Create a new shape Pointers to be used
		dataType ** transPtr = (dataType **)malloc(sizeof(dataType*) * imageHeight); // Transformed Ptr
		dataType ** distTransPtr = (dataType **)malloc(sizeof(dataType*) * imageHeight); // distances for Transformed Ptr
		for (i = 0; i < imageHeight; i++)
		{
			transPtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
			distTransPtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
		}
		// initialize3dArrayD(transPtr, imageLength, imageWidth, imageHeight, params.imageBackground);
		// initialize3dArrayD(distTransPtr, imageLength, imageWidth, imageHeight, 0);
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
		fastSweepingFunction_3D(distTransPtr, transPtr, imageLength, imageWidth, imageHeight, 1, 50000, 0);
		//==============================================================================
		for (k = 0; k < imageHeight; k++)
		{
			free(transPtr[k]);
		}
		free(transPtr);
		transPtr = NULL;
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
		energyTmp = energyFunction(destPtr, distTransPtr, imageHeight, imageLength, imageWidth, params.h);
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
			}
			free(destPtr);
			destPtr = NULL;
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
			affineTmp = gradientComponents(destPtr, distTransPtr, 1.0, &affineResult, imageHeight, imageLength, imageWidth);
			//==============================================================================
			// Free
			for (k = 0; k < imageHeight; k++)
			{
				free(distTransPtr[k]);
			}
			free(distTransPtr);
			distTransPtr = NULL;
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
	size_t k, i, j, l, x, dim2D = imageLength * imageWidth, iteration = 0;
	const dataType h = 1.0;
	dataType firstCpuTime, secondCpuTime, regStartCpuTime, regStopCpuTime, regTotalCpuTimen = 0.;
	dataType energyTotalCpuTime = 0., distanceTotalCpuTime = 0., gradientTotalCpuTime = 0., transformationTotalCpuTime = 0.;
	//==============================================================================
	// Affine Parameters
	Affine_Parameter affineResult;
	//==============================================================================
	// Create new fixed dist. Pointers to be used
	dataType ** destPtr = (dataType **)malloc(sizeof(dataType *) * imageHeight); // distances for destination
	//==============================================================================
	// Initializations of Pointers
	for (i = 0; i < imageHeight; i++)
	{
		destPtr[i] = (dataType *)malloc(sizeof(dataType) * dim2D);
	}
	//==============================================================================
	// Initialize to same background - default value is 255, 0, 0
	// initialize3dArrayD(transPtr, imageLength, imageWidth, imageHeight, params.imageBackground);
	// initialize3dArrayD(distTransPtr, imageLength, imageWidth, imageHeight, 0);
	// initialize3dArrayD(destPtr, imageLength, imageWidth, imageHeight, 0);
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
	//rouyTourinFunction_3D(destPtr, source, 0.00001, imageLength, imageWidth, imageHeight, 0.4, 1);
	//bruteForceFunction_3D(destPtr, source, imageLength, imageWidth, imageHeight, 100000, 0);
	//fastMarching(destPtr, source, imageHeight, imageLength, imageWidth);
	fastSweepingFunction_3D(destPtr, fixedData, imageLength, imageWidth, imageHeight, params.h, 50000, params.imageForeground);
	//==============================================================================
#ifdef MEASURE_TIME
	secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
	// Store the time
	distanceTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
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
	while (!stopCond)
	{
		//==============================================================================
		// Create a new shape Pointers to be used
		dataType ** transPtr = (dataType **)malloc(sizeof(dataType*) * imageHeight); // Transformed Ptr
		dataType ** distTransPtr = (dataType **)malloc(sizeof(dataType*) * imageHeight); // distances for Transformed Ptr
		//==============================================================================
		for (i = 0; i < imageHeight; i++)
		{
			transPtr[i] = malloc(sizeof(dataType) * dim2D);
			distTransPtr[i] = malloc(sizeof(dataType) * dim2D);
		}
		//==============================================================================
		// Timing The Transformation Inside the registration function
#ifdef MEASURE_TIME
		firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
		//==============================================================================
		transform3DImage(movingData, transPtr, affineResult.translation, affineResult.scaling, affineResult.rotation, imageHeight, imageLength, imageWidth, params.imageBackground, centroid);
#ifdef MEASURE_TIME
		secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
		// Store the time
		transformationTotalCpuTime += secondCpuTime - firstCpuTime;
#endif
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
		//rouyTourinFunction_3D(distTransPtr, transPtr, 0.00001, imageLength, imageWidth, imageHeight, 0.4, 1);
		//bruteForceFunction_3D(distTransPtr, transPtr, imageLength, imageWidth, imageHeight, 100000, 0);
		//fastMarching(distTransPtr, transPtr, imageHeight, imageLength, imageWidth);
		fastSweepingFunction_3D(distTransPtr, transPtr, imageLength, imageWidth, imageHeight, params.h, 50000, params.imageForeground);
		//==============================================================================
		// Free transptr after creating the distTransptr
		for (k = 0; k < imageHeight; k++)
		{
			free(transPtr[k]);
		}
		free(transPtr);
		transPtr = NULL;
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
		energyTmp = energyFunction(destPtr, distTransPtr, imageHeight, imageLength, imageWidth, params.h);
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
			// Loop throguh some of the data
#ifdef MEASURE_TIME
			firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
#endif
			//==============================================================================
			// Randomly generated x,y,z points
			//==============================================================================
			// Randomizing random number generation
			//srand(time(NULL));
			//==============================================================================
			bool loop = true;
			l = 0;
			//==============================================================================
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
				// 2D to 1D representation for i, j
				x = x_new(i, j, imageLength);
				// Check if inside the narrow band
				if (NFunction(destPtr[k][x], distTransPtr[k][x], NDelta) == 1) // Checks if inside the band for random generated points
				{
					// Start loop
					l = l + 1;
					//==============================================================================
					// Store the distance function difference
					distDifference = (destPtr[k][x] - distTransPtr[k][x]) * 2.0;
					//==============================================================================
					// Directional component vector derivatives - i, j, k
					dataType tmpI = i / (dataType)imageLength, tmpJ = j / (dataType)imageWidth, tmpK = k / (dataType)imageHeight;
					//==============================================================================
					// Trignometry functions inside the component evaluation equation
					dataType a = cos(phi), b = sin(phi), ab = cos(phi)*sin(phi), aa = cos(phi)*cos(phi), bb = sin(phi)*sin(phi), aaa = a * aa, bbb = b * bb, aab = aa * b, abb = a * bb;
					//==============================================================================
					// Apply Forward Differences to the distTrans pointer
					xFwd = finiteDifX(distTransPtr, h, x, k, i, imageLength);
					yFwd = finiteDifY(distTransPtr, h, k, i, j, imageLength, imageWidth);
					zFwd = finiteDifZ(distTransPtr, h, x, k, i, imageLength, imageHeight);
					//==============================================================================
					// Evaluate Individual Gradient Components
					//==============================================================================

#ifdef UNIFORM
					// Uniform Rotations Angles
					component = xFwd * ((-2.0*tmpI)*(ab / sx) + ((tmpJ)*((bb - aa) / sx)) + ((tmpK)*(a / sx))) +
						yFwd * (((tmpI)*((aa + 2.0 * aab - bb - bbb) / sy)) + ((tmpJ)*((-2.0 * ab - 3.0 * abb) / sy)) + ((tmpK)*((bb - aa) / sy))) +
						zFwd * (((tmpI)*((-aaa + 2.0 * ab + 2.0 * abb) / sz)) + ((tmpJ)*((aa + 2.0 * aa - bb - bbb) / sz)) + ((tmpK)*((-2.0 * ab) / sz)));
					// Set the Rotation - Uniform Angles
					affineResult.rotation.x += rotation_weight * step_size*(component)*(distDifference);
					affineResult.rotation.y += rotation_weight * step_size*(component)*(distDifference);
					affineResult.rotation.z += rotation_weight * step_size*(component)*(distDifference);
					// Uniform Scale Component
					component = xFwd * (((tmpI)*((-aa) * inv_sx2)) + ((tmpJ)*(ab * inv_sx2)) + ((tmpK)*((-b) * inv_sx2))) +
						yFwd * (((-tmpI)*((ab + abb)  * inv_sy2)) + ((-tmpJ)*((aa - bbb) * inv_sy2)) + ((tmpK)*(ab * inv_sy2))) +
						zFwd * (((-tmpI)*((-aab + bb) * inv_sz2)) + ((-tmpJ)*((ab + abb) * inv_sz2)) + ((-tmpK)*(aa * inv_sz2)));
					// Set the Scales - Uniform Scale
					affineResult.scaling.x += scaling_weight * step_size*(component)*(distDifference);
					affineResult.scaling.y += scaling_weight * step_size*(component)*(distDifference);
					affineResult.scaling.z += scaling_weight * step_size*(component)*(distDifference);

#endif // UNIFORM
					//==============================================================================
#ifdef DIRECTIONAL
					// Set the Rotation - Directional
					componentX = yFwd * (((tmpI)*((-sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)) / sy)) + ((tmpJ)*((-cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta)) / sy)) + ((-tmpK)*((cos(phi)*cos(theta)) / sy))) +
						zFwd * (((tmpI)*((cos(phi)*sin(psi) + cos(psi)*sin(phi)*sin(theta)) / sz)) + ((tmpJ)*((cos(phi)*cos(psi) - sin(phi)*sin(psi)*sin(theta)) / sz)) + ((-tmpK)*((cos(theta) * sin(phi)) / sz)));
					componentY = xFwd * (((-tmpI)*((cos(psi)*sin(theta)) / sx)) + ((tmpJ)*((sin(psi)*sin(theta)) / sx)) + ((tmpK)*((cos(theta)) / sx))) +
						yFwd * (((tmpI)*((cos(psi)*cos(theta)*sin(phi)) / sy)) + ((-tmpJ)*((cos(theta)*sin(phi)*sin(psi)) / sy)) + ((tmpK)*((sin(phi)*sin(theta)) / sy))) +
						zFwd * (((-tmpI)*((cos(phi)*cos(psi)*cos(theta)) / sz)) + ((tmpJ)*((cos(phi)*cos(theta)*sin(psi)) / sz)) + ((-tmpK)*((cos(phi)*sin(theta)) / sz)));
					componentZ = xFwd * (((-tmpI)*((cos(theta)*sin(psi)) / sx)) + ((-tmpJ)*((cos(psi)*cos(theta)) / sx))) +
						yFwd * (((tmpI)*((cos(phi)*cos(psi) - sin(phi)*sin(psi)*sin(theta)) / sy)) + ((tmpJ)*((-cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)) / sy))) +
						zFwd * (((tmpI)*((cos(psi)*sin(phi) + cos(phi)*sin(psi)*sin(theta)) / sz)) + ((tmpJ)*((-sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)) / sz)));
					//==============================================================================
					// Set the Rotations - Directional
					affineResult.rotation.x += (params.rotation_weight * (step_size)*(componentX)*(distDifference)) / params.rand_points;
					affineResult.rotation.y += (params.rotation_weight * (step_size)*(componentY)*(distDifference)) / params.rand_points;
					affineResult.rotation.z += (params.rotation_weight * (step_size)*(componentZ)*(distDifference)) / params.rand_points;
					//==============================================================================
					// Set the scales - Directional
					componentX = xFwd * ((-tmpI)*((cos(psi)*cos(theta)) / (sx*sx)) + ((tmpJ)*((cos(theta)*sin(psi)) / (sx*sx))) + ((-tmpK)*((sin(theta)) / (sx*sx))));
					componentY = yFwd * (((-tmpI)*((cos(phi)*sin(psi) + cos(psi)*sin(phi)*sin(theta)) / (sy*sy))) + ((-tmpJ)*((cos(phi)*cos(psi) - sin(phi)*sin(psi)*sin(theta)) / (sy*sy))) + ((tmpK)*((cos(theta)*sin(phi)) / (sy*sy))));
					componentZ = zFwd * (((-tmpI)*((sin(phi)*sin(psi) - cos(phi)*cos(psi)*sin(theta)) / (sz*sz))) + ((-tmpJ)*((cos(psi)*sin(phi) + cos(phi)*sin(psi)*sin(theta)) / (sz*sz))) + ((-tmpK)*((cos(phi)*cos(theta)) / (sz*sz))));
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
			} while ((loop) && (l <= params.rand_points));
			//==============================================================================
			// Free up the copy pointers after calculating the gradients for the next iteration
			for (k = 0; k < imageHeight; k++)
			{
				free(distTransPtr[k]);
			}
			free(distTransPtr);
			distTransPtr = NULL;
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