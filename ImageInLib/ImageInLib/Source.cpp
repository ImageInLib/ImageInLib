#pragma warning(disable : 6011)
#pragma warning(disable : 4244)
#pragma warning(disable : 4477)
#pragma warning(disable : 4267)
#pragma warning(disable : 4700)

#include <stdlib.h>
#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <time.h>

#include "common_vtk.h"
#include "file.h"
#include "common_functions.h"
#include "../src/data_initialization.h"
#include "filtering.h"
#include "../src/data_load.h"
#include "../src/data_storage.h"
#include "../src/thresholding.h"
#include "Labelling.h"
#include "template_functions.h"
#include "../src/endianity_bl.h"
#include "../include/morphological_change.h"
#include "../src/shape_registration.h"
#include "../src/segmentation3D_subsurf.h"
#include "../src/segmentation3d_gsubsurf.h"
#include "../src/image_difference.h"
#include "../src/noise_generator.h"
#include "../src/image_norm.h"
#include "../src/imageInterpolation.h"
#include "distanceForPathFinding.h"

#define thresmin 995
#define thresmax 1213

#define epsilonFast 0.001


int main() {

	size_t i, j, k, x;

	//image Dimensions
	const size_t Width = 512;
	const size_t Length = 512;
	const size_t Height = 406;
	const size_t dim2D = Width * Length;
	
	//Preparation of image data pointers
	dataType** imageData = (dataType**)malloc(Height * sizeof(dataType*));
	short** image = (short**)malloc(Height * sizeof(short*));
	for (k = 0; k < Height; k++) {
		imageData[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		image[k] = (short*)malloc(dim2D * sizeof(short));
	}

	if (imageData == NULL || image == NULL) {
		return false;
	}

	std::string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/";
	std::string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";
	std::string inputImagePath = inputPath +  "patient2.raw";

	if (load3dArrayRAW<short>(image, Length, Width, Height, inputImagePath.c_str()) == false)
	{
		printf("inputImagePath does not exist\n");
	}

	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				x = x_new(i, j, Length);
				imageData[k][x] = (dataType)image[k][x];
			}
		}
	}

	rescaleNewRange(imageData, Length, Width, Height, 0, 1);
	//Filtering
	FilterMethod method =  MEAN_CURVATURE_FILTER; // NONLINEAR_HEATEQUATION_IMPLICIT; //
	Filter_Parameters filterParm; filterParm = { 1.2, 1.0, 0.1, 1000, 1.5, 0.0004, 0.0001, 0.01, 1, 1, 1000 };
	Image_Data toBeFiltered;
	toBeFiltered.imageDataPtr = imageData; toBeFiltered.height = Height; toBeFiltered.length = Length; toBeFiltered.width = Width;
	filterImage(toBeFiltered, filterParm, method);
	rescaleNewRange(imageData, Length, Width, Height, 0, 4000);
	 
	//std::string filteredImagePath = outputPath + "loaded.raw";
	//store3dRawData<dataType>(imageData, Length, Width, Height, filteredImagePath.c_str());

	//Extract sagital slice
	//----------------------------------
	dataType* OneSliceImage = (dataType*)malloc(Height * Width * sizeof(dataType));
	dataType* rotatedImage = (dataType*)malloc(Height * Width * sizeof(dataType));
	dataType* distanceMap = (dataType*)malloc(Height * Width * sizeof(dataType));
	dataType* potentialPtr = (dataType*)malloc(Height * Width * sizeof(dataType));
	dataType* gradientPtr = (dataType*)malloc(Height * Width * sizeof(dataType));
	dataType* maskThresh = (dataType*)malloc(Height * Width * sizeof(dataType));
	size_t cst = 238; i = 0;
	for (k = 0; k < Height; k++) {
		for (j = 0; j < Width; j++) {
			x = x_new(cst, j, Length);
			OneSliceImage[i] = imageData[k][x];
			distanceMap[i] = 0;
			i++;
		}
	}

	for (i = 0; i < Height; i++) {
		for (j = 0; j < Width; j++) {
			x = x_new(i, j, Height);
			rotatedImage[x] = OneSliceImage[x_new(Height - i - 1, Width - j - 1, Height)];
			potentialPtr[x] = 0;
			gradientPtr[x] = 0;
		}
	}

	////Manual thresholding
	//for (k = 0; k < Height; k++) {
	//	for (j = 0; j < Width; j++) {
	//		if (rotatedImage[x_new(k, j, Height)] >= thresmin && rotatedImage[x_new(k, j, Height)] <= thresmax) {
	//			maskThresh[x_new(k, j, Height)] = 1;
	//		}
	//		else {
	//			maskThresh[x_new(k, j, Height)] = 0;
	//		}
	//	}
	//}

	//std::string thresOutPut = outputPath + "initialImage.raw";
	//store2dRawData<dataType>(rotatedImage, Height, Width, thresOutPut.c_str());

	Point2D * startPoint = (Point2D*)malloc(sizeof(Point2D));
	//startPoint->x = 0; startPoint->y = 0;
	startPoint->x = 246; startPoint->y = 277; // x = 238
	//startPoint->x = 209; startPoint->y = 297; // x = 235
	Point2D* endPoint = (Point2D*)malloc(2 * sizeof(Point2D));
	endPoint->x = 250; endPoint->y = 135; // x = 238
	//endPoint->x = 263; endPoint->y = 78;  // x = 235

	Point2D* pointsPath = (Point2D*)malloc(2 * sizeof(Point2D));
	pointsPath[0].x = 246; pointsPath[0].y = 277;
	pointsPath[1].x = 250; pointsPath[1].y = 135;

	////Manual rescalling
	//dataType scale_factor = 1.0 / 4000.0;
	//for (i = 0; i < Height; i++) {
	//	for (j = 0; j < Width; j++) {
	//		x = x_new(j, i, Width);
	//		rotatedImage[x] = scale_factor * (rotatedImage[x] - 4000) + 1;
	//	}
	//}
	
	fastMarching2d(rotatedImage, distanceMap, potentialPtr, Height, Width, pointsPath);
	//fastMarching2d(maskThresh, distanceMap, potentialPtr, Height, Width, pointsPath);

	//std::string outPutDistance = outputPath + "distance.raw";
	//store2dRawData<dataType>(distanceMap, Height, Width, outPutDistance.c_str());

	//std::string potentialPath = outputPath + "potential.raw";
	//store2dRawData<dataType>(potentialPtr, Height, Width, potentialPath.c_str());

	shortestPath2d(distanceMap, gradientPtr, Height, Width, 1.0, pointsPath);
	 
	//std::string outPath = outputPath + "path.raw";
	//store2dRawData<dataType>(gradientPtr, Height, Width, outPath.c_str());


	////========================Add start and Point =================
	//distanceMap[x_new(startPoint->x, startPoint->y - 1, Width)] = 0;
	//distanceMap[x_new(startPoint->x, startPoint->y - 1, Width)] = 0;
	//distanceMap[x_new(startPoint->x, startPoint->y + 1, Width)] = 0;
	//distanceMap[x_new(startPoint->x - 1, startPoint->y, Width)] = 0;
	//distanceMap[x_new(startPoint->x - 1, startPoint->y - 1, Width)] = 0;
	//distanceMap[x_new(startPoint->x - 1, startPoint->y + 1, Width)] = 0;
	//distanceMap[x_new(startPoint->x + 1, startPoint->y, Width)] = 0;
	//distanceMap[x_new(startPoint->x + 1, startPoint->y - 1, Width)] = 0;
	//distanceMap[x_new(startPoint->x + 1, startPoint->y + 1, Width)] = 0;
	//distanceMap[x_new(endPoint->x, endPoint->y, Width)] = 0;
	//distanceMap[x_new(endPoint->x, endPoint->y - 1, Width)] = 0;
	//distanceMap[x_new(endPoint->x, endPoint->y + 1, Width)] = 0;
	//distanceMap[x_new(endPoint->x - 1, endPoint->y - 1, Width)] = 0;
	//distanceMap[x_new(endPoint->x - 1, endPoint->y + 1, Width)] = 0;
	//distanceMap[x_new(endPoint->x - 1, endPoint->y, Width)] = 0;
	//distanceMap[x_new(endPoint->x + 1, endPoint->y - 1, Width)] = 0;
	//distanceMap[x_new(endPoint->x + 1, endPoint->y + 1, Width)] = 0;
	//distanceMap[x_new(endPoint->x + 1, endPoint->y, Width)] = 0;
	////=============================================================

	for (i = 0; i < Height; i++) {
		for (j = 0; j < Width; j++) {
			x = x_new(j, i, Width);
			if (gradientPtr[x] == 1) {
				distanceMap[x] = 0;
				rotatedImage[x] = 0;
				maskThresh[x] = 0;
			}
		}
	}
	//outPutDistance = outputPath + "distance_plus_path_.raw";
	//store2dRawData<dataType>(distanceMap, Height, Width, outPutDistance.c_str());

	std::string originalImage = outputPath + "Initial_plus_path.raw";
	store2dRawData<dataType>(rotatedImage, Height, Width, originalImage.c_str());
	//store2dRawData<dataType>(maskThresh, Height, Width, originalImage.c_str());

	//computeImageGradient(distanceMap, gradientPtr, Height, Width, 1.0);
	//std::string outPutGradient = outputPath + "gradient.raw";
	//store2dRawData<dataType>(gradientPtr, Height, Width, outPutGradient.c_str());
	//outPutDistance = outputPath + "potential.raw";
	//store2dRawData<dataType>(potentialPtr, Height, Width, outPutDistance.c_str());

	//computeImageGradient(distanceMap, gradientPtr, Height, Width, 1.0);
	//std::string outPutDistance = outputPath + "gradient_DM.raw";
	//store2dRawData<dataType>(gradientPtr, Height, Width, outPutDistance.c_str());

	//free memory
	for (k = 0; k < Height; k++) {
		free(imageData[k]); free(image[k]);
	}
	free(imageData); free(image);

	free(OneSliceImage); free(rotatedImage); free(distanceMap); 
	free(potentialPtr);  free(gradientPtr);
	free(startPoint); free(endPoint);
	free(pointsPath); free(maskThresh);

	return EXIT_SUCCESS;
}
