/*
* Author : Konan Allaly
* Purpose : Updates for the INFLANET Project
*/

#include <stdlib.h>
#include <string.h>
#include <stdbool.h> 
#include <math.h>
#include "segmentation2d.h"

dataType l2norm(dataType* arrayPtr1, dataType* arrayPtr2, const size_t height, const size_t width, dataType h) {
	size_t i;

	dataType sumPower = 0.0, norm = 0.0;
	dataType hh = h * h;

	for (i = 0; i < height * width; i++) {
		sumPower += (dataType)((pow(arrayPtr1[i] - arrayPtr2[i], 2) * hh));
	}
	norm = sqrt(sumPower);

	return norm;
}

bool rescaleToZeroOne2d(dataType* imageDataPtr, const size_t height, const size_t width)
{
	//check if the memory was allocated successfully
	if (imageDataPtr == NULL)
		return false;

	size_t i;
	dataType max = 0, min = 100000, quotient, offset;

	//Determine minimum and maximum value
	for (i = 0; i < height * width; i++) {
		if (imageDataPtr[i] < min)
			min = imageDataPtr[i];
		if (imageDataPtr[i] > max)
			max = imageDataPtr[i];
	}

	quotient = (dataType)1. / (max - min);
	offset = min * quotient;
	//Rescale values to interval (0, 1)
	for (i = 0; i < height * width; i++) {
		imageDataPtr[i] = (quotient * imageDataPtr[i] - offset);
	}

	return true;
}

bool generateInitialSegmentationFunction(dataType* imageDataPtr, const size_t height, const size_t width, point2d* center, dataType v, dataType R)
{
	size_t i, j;
	int dx, dy;
	dataType norm_of_distance = 0.0, new_value = 0.0;

	if (imageDataPtr == NULL)
		return false;

	for (i = 0; i < height; i++) {
		dx = i - center->x;
		for (j = 0; j < width; j++) {
			dy = j - center->y;
			norm_of_distance = sqrt(dx * dx + dy * dy);
			new_value = (dataType)((1.0 / (norm_of_distance + v)) - (1.0 / (R + v)));
			if (norm_of_distance > R) {
				imageDataPtr[x_new(i, j, height)] = 0;
			}
			else {
				imageDataPtr[x_new(i, j, height)] = new_value;
			}
		}
	}

	rescaleToZeroOne2d(imageDataPtr, height, width);

	return true;
}

bool set2dDirichletBoundaryCondition(dataType* imageDataPtr, const size_t height, const size_t width) {
	size_t i, j;
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			if (i == 0 || i == height - 1 || j == 0 || j == width - 1) {
				imageDataPtr[x_new(i, j, height)] = 0.0;
			}
		}
	}
	return true;
}

bool computeNormOfGradientDiamondCells(dataType* arrayPtr, neighPtrs neigbours, const size_t height, const size_t width, dataType h) {

	size_t i, j, i_ext, j_ext, currentIndx;
	const size_t height_ext = height + 2, width_ext = width + 2;
	size_t dim2D_ext = height_ext * width_ext;

	dataType* extendedArray = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
	if (extendedArray == NULL)
		return false;

	copyDataTo2dExtendedArea(arrayPtr, extendedArray, height, width);
	reflection2D(extendedArray, height_ext, width_ext);

	dataType uP, uN, uNW, uNE, uS, uSW, uSE, uW, uE;
	dataType ux, uy;

	for (i = 0, i_ext = 1; i < height; i++, i_ext++) {
		for (j = 0, j_ext = 1; j < width; j++, j_ext++) {

			size_t iplus = i_ext + 1;
			size_t iminus = i_ext - 1;
			size_t jplus = j_ext + 1;
			size_t jminus = j_ext - 1;

			currentIndx = x_new(i, j, height);
			uP = extendedArray[x_new(i_ext, j_ext, height_ext)];
			uE = extendedArray[x_new(iplus, j_ext, height_ext)];
			uW = extendedArray[x_new(iminus, j_ext, height_ext)];
			uN = extendedArray[x_new(i_ext, jminus, height_ext)];
			uS = extendedArray[x_new(i_ext, jplus, height_ext)];
			uNE = extendedArray[x_new(iplus, jminus, height_ext)];
			uNW = extendedArray[x_new(iminus, jminus, height_ext)];
			uSE = extendedArray[x_new(iplus, jplus, height_ext)];
			uSW = extendedArray[x_new(iminus, jplus, height_ext)];

			//East
			ux = (uE - uP) / h;
			uy = (uNE + uN - uS - uSE) / (4.0 * h);
			neigbours.East[currentIndx] = sqrt(ux * ux + uy * uy);

			//West
			ux = (uP - uW) / h;
			uy = (uNW + uN - uSW - uS) / (4.0 * h);
			neigbours.West[currentIndx] = sqrt(ux * ux + uy * uy);

			//North
			ux = (uNE + uE - uNW - uW) / (4.0 * h);
			uy = (uN - uP) / h;
			neigbours.North[currentIndx] = sqrt(ux * ux + uy * uy);

			//South
			ux = (uSE + uE - uSW - uW) / (4.0 * h);
			uy = (uP - uS) / h;
			neigbours.South[currentIndx] = sqrt(ux * ux + uy * uy);
		}
	}

	free(extendedArray);

	return true;
}

bool epsilonRegularization(neighPtrs neighbours, const size_t height, const size_t width, dataType epsilon) {
	size_t i, j, currentIndx;
	dataType current = 0.0;
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			currentIndx = x_new(i, j, height);

			current = neighbours.East[currentIndx];
			neighbours.East[currentIndx] = (dataType)(sqrt(current * current + epsilon));

			current = neighbours.West[currentIndx];
			neighbours.West[currentIndx] = (dataType)(sqrt(current * current + epsilon));

			current = neighbours.North[currentIndx];
			neighbours.North[currentIndx] = (dataType)(sqrt(current * current + epsilon));

			current = neighbours.South[currentIndx];
			neighbours.South[currentIndx] = (dataType)(sqrt(current * current + epsilon));
		}
	}
	return true;
}

bool subsurf(Image_Data2D imageData, dataType* initialSegment, const char* segmentPath, const Filter_Parameters smooth_parms, Segmentation_Parameters seg_parms)
{
	size_t i, j, i_ext, j_ext;
	size_t height = imageData.height, width = imageData.width;
	const size_t height_ext = height + 2, width_ext = width + 2;
	size_t dim2D = height * width, dim2D_ext = height_ext * width_ext;

	dataType tau = seg_parms.tau, h = seg_parms.h;
	dataType tol = seg_parms.segTolerance, omega = seg_parms.omega_c;
	dataType coef_edge_detector = seg_parms.coef, eps = seg_parms.eps2;
	dataType coef_tau = tau / (h * h), gauss_seidel_coef = 0.0;

	size_t maxIter = seg_parms.maxNoGSIteration;

	dataType* segmentationPtr = (dataType*)malloc(sizeof(dataType) * dim2D);

	dataType* gaussSeidelPtr = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
	dataType* previousSolPtr = (dataType*)malloc(sizeof(dataType) * dim2D_ext);

	if (segmentationPtr == NULL || gaussSeidelPtr == NULL || previousSolPtr == NULL)
		return false;

	heatImplicit2dScheme(imageData, smooth_parms);

	dataType* uNorth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* uSouth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* uEast = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* uWest = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (uNorth == NULL || uSouth == NULL || uEast == NULL || uWest == NULL)
		return false;

	neighPtrs  U;
	U.West = uWest;
	U.East = uEast;
	U.North = uNorth;
	U.South = uSouth;

	dataType* gNorth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* gSouth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* gEast = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* gWest = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* gAverage = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (gNorth == NULL || gSouth == NULL || gEast == NULL || gWest == NULL || gAverage == NULL)
		return false;

	computeNormOfGradientDiamondCells(imageData.imageDataPtr, U, height, width, h);

	dataType current = 0.0;
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {

			size_t currentIndx = x_new(i, j, height);

			gEast[currentIndx] = gradientFunction(pow(U.East[currentIndx], 2), coef_edge_detector);
			gWest[currentIndx] = gradientFunction(pow(U.West[currentIndx], 2), coef_edge_detector);
			gNorth[currentIndx] = gradientFunction(pow(U.North[currentIndx], 2), coef_edge_detector);
			gSouth[currentIndx] = gradientFunction(pow(U.South[currentIndx], 2), coef_edge_detector);

			current = (dataType)(((gEast[currentIndx] + gWest[currentIndx] + gNorth[currentIndx] + gSouth[currentIndx]) / 4.0));
			gAverage[currentIndx] = current;
		}
	}

	//Array for name construction
	char name[350];
	char name_ending[100];
	Storage_Flags flags = { false,false };

	strcpy_s(name, sizeof name, segmentPath);
	sprintf_s(name_ending, sizeof(name_ending), "_edgeDetector.raw");
	strcat_s(name, sizeof(name), name_ending);
	store2dRawData(gAverage, height, width, name, flags);

	dataType* coefNorth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* coefSouth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* coefEast = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* coefWest = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (coefNorth == NULL || coefSouth == NULL || coefEast == NULL || coefWest == NULL)
		return false;

	dataType average_norm_gradient, u_average;

	copyDataToAnother2dArray(initialSegment, segmentationPtr, height, width);

	copyDataTo2dExtendedArea(initialSegment, previousSolPtr, height, width);
	set2dDirichletBoundaryCondition(previousSolPtr, height_ext, width_ext);
	copyDataTo2dExtendedArea(initialSegment, gaussSeidelPtr, height, width);
	set2dDirichletBoundaryCondition(gaussSeidelPtr, height_ext, width_ext);

	//segmentation loop
	size_t number_time_step = 0;
	dataType error_segmentation = 0.0;

	do {
		number_time_step++;

		//compute the coefficents
		computeNormOfGradientDiamondCells(segmentationPtr, U, height, width, h);
		epsilonRegularization(U, height, width, eps);
		for (i = 0; i < height; i++) {
			for (j = 0; j < width; j++) {
				size_t currentIndx = x_new(i, j, height);

				average_norm_gradient = (dataType)((U.East[currentIndx] + U.West[currentIndx] + U.North[currentIndx] + U.South[currentIndx]) / 4.0);
				u_average = sqrt(average_norm_gradient * average_norm_gradient + eps * eps);

				coefEast[currentIndx] = coef_tau * u_average * gEast[currentIndx] * (1.0 / U.East[currentIndx]);
				coefNorth[currentIndx] = coef_tau * u_average * gNorth[currentIndx] * (1.0 / U.North[currentIndx]);
				coefWest[currentIndx] = coef_tau * u_average * gWest[currentIndx] * (1.0 / U.West[currentIndx]);
				coefSouth[currentIndx] = coef_tau * u_average * gSouth[currentIndx] * (1.0 / U.South[currentIndx]);
			}
		}

		//gauss seidel for segmentation function
		size_t cpt = 0;
		dataType error_gauss_seidel = 0.0;
		do {
			cpt++;
			for (i = 0, i_ext = 1; i < height; i++, i_ext++) {
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) {

					size_t currentIndx = x_new(i, j, height);
					size_t currentIndx_ext = x_new(i_ext, j_ext, height_ext);

					gauss_seidel_coef = (dataType)((previousSolPtr[currentIndx_ext] + coefEast[currentIndx] * gaussSeidelPtr[x_new(i_ext + 1, j_ext, height_ext)] + coefNorth[currentIndx] * gaussSeidelPtr[x_new(i_ext, j_ext - 1, height_ext)]
						+ coefWest[currentIndx] * gaussSeidelPtr[x_new(i_ext - 1, j_ext, height_ext)] + coefSouth[currentIndx] * gaussSeidelPtr[x_new(i_ext, j_ext + 1, height_ext)])
						/ (1 + coefEast[currentIndx] + coefNorth[currentIndx] + coefWest[currentIndx] + coefSouth[currentIndx]));

					gaussSeidelPtr[currentIndx_ext] = gaussSeidelPtr[currentIndx_ext] + omega * (gauss_seidel_coef - gaussSeidelPtr[currentIndx_ext]);
				}
			}

			error_gauss_seidel = 0.0;
			for (i = 0, i_ext = 1; i < height; i++, i_ext++) {
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) {

					size_t currentIndx = x_new(i, j, height);
					size_t currentIndx_ext = x_new(i_ext, j_ext, height_ext);

					error_gauss_seidel += (dataType)(pow((1 + coefEast[currentIndx] + coefNorth[currentIndx] + coefWest[currentIndx] + coefSouth[currentIndx]) * gaussSeidelPtr[currentIndx_ext]
						- (coefEast[currentIndx] * gaussSeidelPtr[x_new(i_ext + 1, j_ext, height_ext)] + coefNorth[currentIndx] * gaussSeidelPtr[x_new(i_ext, j_ext - 1, height_ext)]
							+ coefWest[currentIndx] * gaussSeidelPtr[x_new(i_ext - 1, j_ext, height_ext)] + coefSouth[currentIndx] * gaussSeidelPtr[x_new(i_ext, j_ext + 1, height_ext)]) - previousSolPtr[currentIndx_ext], 2) * h * h);
				}
			}
		} while (cpt < maxIter && error_gauss_seidel > tol);

		//rescall to data range 0-1
		rescaleToZeroOne2d(gaussSeidelPtr, height_ext, width_ext);

		//compute L2-norm
		error_segmentation = l2norm(gaussSeidelPtr, previousSolPtr, height_ext, width_ext, h);

		printf("Step  %zd : residu = %e \n", number_time_step, error_segmentation);

		//Dirichlet Boundary condition
		set2dDirichletBoundaryCondition(gaussSeidelPtr, height_ext, width_ext);

		//copy
		copyDataToAnother2dArray(gaussSeidelPtr, previousSolPtr, height_ext, width_ext);

		//copy to reduce array
		copyDataTo2dReducedArea(segmentationPtr, gaussSeidelPtr, height, width);

		//save the solution
		if (number_time_step % seg_parms.mod == 0) {
			strcpy_s(name, sizeof name, segmentPath);
			sprintf_s(name_ending, sizeof(name_ending), "_seg_func_%03zd.raw", number_time_step);
			strcat_s(name, sizeof(name), name_ending);
			store2dRawData(segmentationPtr, height, width, name, flags);
		}

	} while (number_time_step <= seg_parms.maxNoOfTimeSteps && error_segmentation > tol);

	free(uNorth);
	free(uSouth);
	free(uEast);
	free(uWest);

	free(gNorth);
	free(gSouth);
	free(gEast);
	free(gWest);
	free(gAverage);

	free(coefNorth);
	free(coefSouth);
	free(coefEast);
	free(coefWest);

	free(segmentationPtr);
	free(gaussSeidelPtr);
	free(previousSolPtr);

	return true;
}

bool gsubsurf(Image_Data2D imageData, dataType* initialSegment, const char* segmentPath, const Filter_Parameters smooth_parms, Segmentation_Parameters seg_parms)
{
	size_t i, j, i_ext, j_ext;
	const size_t height = imageData.height, width = imageData.width;
	const size_t height_ext = height + 2, width_ext = width + 2;
	size_t dim2D = height * width, dim2D_ext = height_ext * width_ext;

	dataType tau = seg_parms.tau, h = seg_parms.h;
	dataType tol = seg_parms.segTolerance, omega = seg_parms.omega_c;
	dataType coef_edge_detector = seg_parms.coef, eps = seg_parms.eps2;
	dataType coef_tau = tau / (h * h), gauss_seidel_coef = 0.0;
	dataType diff = seg_parms.coef_dif, adv = seg_parms.coef_conv;
	size_t maxIter = seg_parms.maxNoGSIteration;

	dataType* segmentationPtr = (dataType*)malloc(sizeof(dataType) * dim2D);

	dataType* gaussSeidelPtr = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
	dataType* previousSolPtr = (dataType*)malloc(sizeof(dataType) * dim2D_ext);

	if (segmentationPtr == NULL || gaussSeidelPtr == NULL || previousSolPtr == NULL)
		return false;

	dataType* uNorth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* uSouth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* uEast = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* uWest = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (uNorth == NULL || uSouth == NULL || uEast == NULL || uWest == NULL)
		return false;

	dataType* edgeDetectorPtr = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (edgeDetectorPtr == NULL)
		return false;

	neighPtrs  uCoef;
	uCoef.West = uWest;
	uCoef.East = uEast;
	uCoef.North = uNorth;
	uCoef.South = uSouth;

	dataType* vNorth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* vSouth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* vEast = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* vWest = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (vNorth == NULL || vSouth == NULL || vEast == NULL || vWest == NULL)
		return false;

	dataType* coefNorth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* coefSouth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* coefEast = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* coefWest = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (coefNorth == NULL || coefSouth == NULL || coefEast == NULL || coefWest == NULL)
		return false;

	dataType current, average_norm_gradient, average_gFunction;
	dataType u_average;

	heatImplicit2dScheme(imageData, smooth_parms);

	//compute g function
	computeNormOfGradientDiamondCells(imageData.imageDataPtr, uCoef, height, width, h);
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			size_t currentIndx = x_new(i, j, height);
			average_gFunction = (dataType)((uCoef.East[currentIndx] + uCoef.West[currentIndx] + uCoef.North[currentIndx] + uCoef.South[currentIndx]) / 4.0);
			edgeDetectorPtr[currentIndx] = gradientFunction(pow(average_gFunction, 2), coef_edge_detector);
		}
	}

	//Array for name construction
	char name[350];
	char name_ending[100];
	Storage_Flags flags = { false,false };

	strcpy_s(name, sizeof name, segmentPath);
	sprintf_s(name_ending, sizeof(name_ending), "_edgeDetector.raw");
	strcat_s(name, sizeof(name), name_ending);
	store2dRawData(edgeDetectorPtr, height, width, name, flags);

	//compute gradient of edge detector function
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			size_t currentIndx = x_new(i, j, height);

			if (i == 0) {
				vEast[currentIndx] = -adv * ((edgeDetectorPtr[x_new(i + 1, j, height)] - edgeDetectorPtr[currentIndx]) / (2 * h));
				vWest[currentIndx] = -vEast[currentIndx];
			}
			else {
				if (i == height - 1) {
					vEast[currentIndx] = -adv * ((edgeDetectorPtr[currentIndx] - edgeDetectorPtr[x_new(i - 1, j, height)]) / (2 * h));
					vWest[currentIndx] = -vEast[currentIndx];
				}
				else {
					vEast[currentIndx] = -adv * ((edgeDetectorPtr[x_new(i + 1, j, height)] - edgeDetectorPtr[x_new(i - 1, j, height)]) / (2 * h));
					vWest[currentIndx] = -vEast[currentIndx];
				}
			}
			if (j == 0) {
				vSouth[currentIndx] = -adv * ((edgeDetectorPtr[x_new(i, j + 1, height)] - edgeDetectorPtr[currentIndx]) / (2 * h));
				vNorth[currentIndx] = -vSouth[currentIndx];
			}
			else {
				if (j == width - 1) {
					vSouth[currentIndx] = -adv * ((edgeDetectorPtr[currentIndx] - edgeDetectorPtr[x_new(i, j - 1, height)]) / (2 * h));
					vNorth[currentIndx] = -vSouth[currentIndx];
				}
				else {
					vSouth[currentIndx] = -adv * ((edgeDetectorPtr[x_new(i, j + 1, height)] - edgeDetectorPtr[x_new(i, j - 1, height)]) / (2 * h));
					vNorth[currentIndx] = -vSouth[currentIndx];
				}
			}
		}
	}

	copyDataToAnother2dArray(initialSegment, segmentationPtr, height, width);
	copyDataTo2dExtendedArea(initialSegment, gaussSeidelPtr, height, width);
	copyDataTo2dExtendedArea(initialSegment, previousSolPtr, height, width);

	set2dDirichletBoundaryCondition(gaussSeidelPtr, height_ext, width_ext);
	set2dDirichletBoundaryCondition(previousSolPtr, height_ext, width_ext);

	//segmentation loop
	size_t number_time_step = 0;
	dataType error_segmentation = 0.0;
	do {
		number_time_step++;

		computeNormOfGradientDiamondCells(segmentationPtr, uCoef, height, width, h);
		epsilonRegularization(uCoef, height, width, eps);

		for (i = 0; i < height; i++) {
			for (j = 0; j < width; j++) {
				size_t currentIndx = x_new(i, j, height);
				average_norm_gradient = (dataType)((uCoef.East[currentIndx] + uCoef.West[currentIndx] + uCoef.North[currentIndx] + uCoef.South[currentIndx]) / 4.0);
				u_average = sqrt(average_norm_gradient * average_norm_gradient + eps * eps);

				//Explicit for advection
				coefEast[currentIndx] = (dataType)(coef_tau * (-min(vEast[currentIndx], 0) + diff * edgeDetectorPtr[currentIndx] * u_average * (1.0 / uCoef.East[currentIndx])));
				coefNorth[currentIndx] = (dataType)(coef_tau * (-min(vNorth[currentIndx], 0) + diff * edgeDetectorPtr[currentIndx] * u_average * (1.0 / uCoef.North[currentIndx])));
				coefWest[currentIndx] = (dataType)(coef_tau * (-min(vWest[currentIndx], 0) + diff * edgeDetectorPtr[currentIndx] * u_average * (1.0 / uCoef.West[currentIndx])));
				coefSouth[currentIndx] = (dataType)(coef_tau * (-min(vSouth[currentIndx], 0) + diff * edgeDetectorPtr[currentIndx] * u_average * (1.0 / uCoef.South[currentIndx])));
			}
		}

		//gauss seidel for segmentation function
		size_t cpt = 0;

		dataType error_gauss_seidel = 0.0;
		do {
			cpt++;
			for (i = 0, i_ext = 1; i < height; i++, i_ext++) {
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) {

					size_t iplus = i_ext + 1;
					size_t iminus = i_ext - 1;
					size_t jplus = j_ext + 1;
					size_t jminus = j_ext - 1;

					size_t currentIndx = x_new(i, j, height);
					size_t currentIndx_ext = x_new(i_ext, j_ext, height_ext);

					//explicit for advection term
					gauss_seidel_coef = (dataType)((previousSolPtr[currentIndx_ext] + coefEast[currentIndx] * gaussSeidelPtr[x_new(iplus, j_ext, height_ext)]
						+ coefWest[currentIndx] * gaussSeidelPtr[x_new(iminus, j_ext, height_ext)] + coefNorth[currentIndx] * gaussSeidelPtr[x_new(i_ext, jminus, height_ext)]
						+ coefSouth[currentIndx] * gaussSeidelPtr[x_new(i_ext, jplus, height_ext)]) / (1 + coefEast[currentIndx] + coefWest[currentIndx] +
							coefNorth[currentIndx] + coefSouth[currentIndx]));
					gaussSeidelPtr[currentIndx_ext] = gaussSeidelPtr[currentIndx_ext] + omega * (gauss_seidel_coef - gaussSeidelPtr[currentIndx_ext]);
				}
			}

			error_gauss_seidel = 0.0;
			for (i = 0, i_ext = 1; i < height; i++, i_ext++) {
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) {

					size_t iplus = i_ext + 1;
					size_t iminus = i_ext - 1;
					size_t jplus = j_ext + 1;
					size_t jminus = j_ext - 1;

					size_t currentIndx = x_new(i, j, height);
					size_t currentIndx_ext = x_new(i_ext, j_ext, height_ext);

					//Explicit for advection term
					error_gauss_seidel += (dataType)(pow((1 + coefEast[currentIndx] + coefWest[currentIndx] + coefNorth[currentIndx] + coefSouth[currentIndx]) * gaussSeidelPtr[currentIndx_ext]
						- (coefEast[currentIndx] * gaussSeidelPtr[x_new(iplus, j_ext, height_ext)]
							+ coefWest[currentIndx] * gaussSeidelPtr[x_new(iminus, j_ext, height_ext)] + coefNorth[currentIndx] * gaussSeidelPtr[x_new(i_ext, jminus, height_ext)]
							+ coefSouth[currentIndx] * gaussSeidelPtr[x_new(i_ext, jplus, height_ext)]) - previousSolPtr[currentIndx_ext], 2) * h * h);
				}
			}

		} while (cpt < maxIter && error_gauss_seidel > 0.000001);

		//rescall to data range 0-1
		rescaleToZeroOne2d(gaussSeidelPtr, height_ext, width_ext);

		//compute L2-norm
		error_segmentation = l2norm(gaussSeidelPtr, previousSolPtr, height_ext, width_ext, h);

		printf("Step %zd , residual = %e \n", number_time_step, error_segmentation);

		set2dDirichletBoundaryCondition(gaussSeidelPtr, height_ext, width_ext);

		copyDataToAnother2dArray(gaussSeidelPtr, previousSolPtr, height_ext, width_ext);

		//copy to reduce array
		copyDataTo2dReducedArea(segmentationPtr, gaussSeidelPtr, height, width);

		//save the solution
		if (number_time_step % seg_parms.mod == 0) {
			strcpy_s(name, sizeof name, segmentPath);
			sprintf_s(name_ending, sizeof(name_ending), "_seg_func_%03zd.raw", number_time_step);
			strcat_s(name, sizeof(name), name_ending);
			store2dRawData(segmentationPtr, height, width, name, flags);
		}

	} while (number_time_step <= seg_parms.maxNoOfTimeSteps && error_segmentation > seg_parms.segTolerance);

	free(uNorth);
	free(uSouth);
	free(uEast);
	free(uWest);

	free(edgeDetectorPtr);

	free(vNorth);
	free(vSouth);
	free(vEast);
	free(vWest);

	free(coefNorth);
	free(coefSouth);
	free(coefEast);
	free(coefWest);

	free(segmentationPtr);
	free(gaussSeidelPtr);
	free(previousSolPtr);

	return true;
}
