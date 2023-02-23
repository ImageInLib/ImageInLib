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

	////image Dimensions
	//const size_t Width = 512;
	//const size_t Length = 512;
	//const size_t Height = 406;
	//const size_t dim2D = Width * Length;

	//-------------Real 3D image -------------------------
	
	////Preparation of image data pointers
	//dataType** imageData = (dataType**)malloc(Height * sizeof(dataType*));
	//dataType** potential3D = (dataType**)malloc(Height * sizeof(dataType*));
	//short** image = (short**)malloc(Height * sizeof(short*));
	//for (k = 0; k < Height; k++) {
	//	imageData[k] = (dataType*)malloc(dim2D * sizeof(dataType));
	//	potential3D[k] = (dataType*)malloc(dim2D * sizeof(dataType));
	//	image[k] = (short*)malloc(dim2D * sizeof(short));
	//}
	//if (imageData == NULL || image == NULL || potential3D == NULL) {
	//	return false;
	//}
	//std::string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/";
	//std::string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";
	//std::string inputImagePath = inputPath +  "patient2.raw";
	//if (load3dArrayRAW<short>(image, Length, Width, Height, inputImagePath.c_str()) == false)
	//{
	//	printf("inputImagePath does not exist\n");
	//}
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			x = x_new(i, j, Length);
	//			imageData[k][x] = (dataType)image[k][x];
	//			potential3D[k][x] = 0;
	//		}
	//	}
	//}

	//Point3d* seed = (Point3d*)malloc(2 * sizeof(Point3d));
	//seed[0].x = 246; seed[0].y = 277, seed[0].z = 238;
	//seed[1].x = 250; seed[1].y = 135, seed[1].z = 238;

	//compute3dPotential(imageData, potential3D, Length, Width, Height, seed);
	//std::string loadedImagePath = outputPath + "potential.raw";
	//store3dRawData<dataType>(potential3D, Length, Width, Height, loadedImagePath.c_str());

	//-------------Filtering-------------------------

	//rescaleNewRange(imageData, Length, Width, Height, 0, 1);
	//FilterMethod method = NONLINEAR_HEATEQUATION_IMPLICIT; // MEAN_CURVATURE_FILTER; //
	//Filter_Parameters filterParm; filterParm = { 1.2, 1.0, 0.1, 1000, 1.5, 0.0004, 0.0001, 0.01, 1, 1, 1000 };
	//Image_Data toBeFiltered;
	//toBeFiltered.imageDataPtr = imageData; toBeFiltered.height = Height; toBeFiltered.length = Length; toBeFiltered.width = Width;
	//filterImage(toBeFiltered, filterParm, method);
	//rescaleNewRange(imageData, Length, Width, Height, 0, 4000);
	//std::string filteredImagePath = outputPath + "filtered.raw";
	//store3dRawData<dataType>(imageData, Length, Width, Height, filteredImagePath.c_str());

	//---------------Real 2D image -----------------------

	//Extract Slice of interest

	//dataType* real2dImageData = (dataType*)malloc(dim2D * sizeof(dataType));
	//dataType* loadingPtr = (dataType*)malloc(dim2D * sizeof(dataType));
	//dataType* maskThresh = (dataType*)malloc(dim2D * sizeof(dataType));

	//size_t cst = 238; i = 0;
	//for (k = 0; k < Height; k++) {
	//	for (j = 0; j < Width; j++) {
	//		x = x_new(cst, j, Length);
	//		loadingPtr[i] = imageData[k][x];
	//		i++;
	//	}
	//}

	//for (i = 0; i < Height; i++) {
	//	for (j = 0; j < Width; j++) {
	//		real2dImageData[x_new(i, j, Height)] = loadingPtr[x_new(Height - i - 1, Width - j - 1, Height)];
	//	}
	//}

	//std::string real2dImagePath = outputPath + "loaded2dImage.raw";
	//store2dRawData<dataType>(real2dImageData, Height, Width, real2dImagePath.c_str());

	//Point2D* startAndEnd = (Point2D*)malloc(2 * sizeof(Point2D));
	//startAndEnd[1].x = 246; startAndEnd[1].y = 277;
	//startAndEnd[0].x = 250; startAndEnd[0].y = 135;

	////Manual thresholding
	//for (k = 0; k < Height; k++) {
	//	for (j = 0; j < Width; j++) {
	//		if (real2dImageData[x_new(k, j, Height)] >= thresmin && real2dImageData[x_new(k, j, Height)] <= thresmax) {
	//			maskThresh[x_new(k, j, Height)] = 1;
	//		}
	//		else {
	//			maskThresh[x_new(k, j, Height)] = 0;
	//		}
	//	}
	//}

	//-----------Create artificial image--------------------

	const size_t Height = 200;
	const size_t Width = 200;
	const size_t dim2D = Width * Height;
	 
	dataType* artificial2dImage = (dataType*)malloc(dim2D * sizeof(dataType));
	std::string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";

	//Draw snake
	for (i = 0; i < Height; i++) {
		for (j = 0; j < Width; j++) {
			artificial2dImage[x_new(j, i, Width)] = 0;
		}
	}
	size_t i1 = 50, j1 = 100, radius = 50;
	for (i = 0; i < 100; i++) {
		for (j = j1; j < ( j1 + radius) ; j++) {
			if (sqrt((i1 - i) * (i1 - i) + (j1 - j) * (j1 - j)) <= radius) {
				artificial2dImage[x_new(j, i, Width)] = 1;
			}
			if (sqrt((i1 - i) * (i1 - i) + (j1 - j) * (j1 - j)) <= 20) {
				artificial2dImage[x_new(j, i, Width)] = 0;
			}
		}
	}
	size_t i2 = 120;
	for (i = 70; i < 170; i++) {
		for (j = 50; j < 100; j++) {
			if (sqrt((i2 - i) * (i2 - i) + (j1 - j) * (j1 - j)) <= radius) {
				artificial2dImage[x_new(j, i, Width)] = 1;
			}
			if (sqrt((i2 - i) * (i2 - i) + (j1 - j) * (j1 - j)) <= 20) {
				artificial2dImage[x_new(j, i, Width)] = 0;
			}
		}
	}
	Point2D* startAndEnd = (Point2D*)malloc(2 * sizeof(Point2D));
	startAndEnd[1].x = 113; startAndEnd[1].y = 18;
	startAndEnd[0].x = 79; startAndEnd[0].y = 148;

	////Draw L
	//for (i = 0; i < Height; i++) {
	//	for (j = 0; j < Width; j++) {
	//		x = x_new(j, i, Width);
	//		artificial2dImage[x] = 0;
	//	}
	//}
	//for (i = 0; i < Height; i++) {
	//	for (j = 0; j < Width; j++) {
	//		x = x_new(j, i, Width);
	//		if ((j >= 10 && j <= 40) && (i >= 10 && i <= 190)) {
	//			artificial2dImage[x] = 1;
	//		}
	//		if ((j >= 10 && j <= 190) && (i >= 160 && i <= 190)) {
	//			artificial2dImage[x] = 1;
	//		}
	//	}
	//}
	//Point2D* startAndEnd = (Point2D*)malloc(2 * sizeof(Point2D));
	//startAndEnd[1].x = 25; startAndEnd[1].y = 15;
	//startAndEnd[0].x = 181; startAndEnd[0].y = 175;
	
	//std::string pathArtificialImage = outputPath + "artificial.raw";
	//store2dRawData<dataType>(artificial2dImage, Height, Width, pathArtificialImage.c_str());

	//------------------------------------------------------

	dataType* distanceMap = (dataType*)malloc(dim2D * sizeof(dataType));
	dataType* potentialPtr = (dataType*)malloc(dim2D * sizeof(dataType));
	dataType* pathPtr = (dataType*)malloc(dim2D * sizeof(dataType));


	for (i = 0; i < Height; i++) {
		for (j = 0; j < Width; j++) {
			x = x_new(j, i, Width);
			potentialPtr[x] = 0;
			distanceMap[x] = 0;
			pathPtr[x] = 0;
		}
	}

	fastMarching2d(artificial2dImage, distanceMap, potentialPtr, Height, Width, startAndEnd);
	//fastMarching2d(real2dImageData, distanceMap, potentialPtr, Height, Width, startAndEnd);
	//fastMarching2d(maskThresh, distanceMap, potentialPtr, Height, Width, startAndEnd);

	std::string outPath = outputPath + "distanceMap.raw";
	store2dRawData<dataType>(distanceMap, Height, Width, outPath.c_str());
	//outPath = outputPath + "Potential.raw";
	//store2dRawData<dataType>(potentialPtr, Height, Width, outPath.c_str());

	//std::string distanceOutPut = outputPath + "distanceMap.raw";
	//store2dRawData<dataType>(distanceMap, Height, Width, distanceOutPut.c_str());

	shortestPath2d(distanceMap, pathPtr, Height, Width, 1.0, startAndEnd);

	for (i = 0; i < Height; i++) {
		for (j = 0; j < Width; j++) {
			x = x_new(j, i, Width);
			if (pathPtr[x] == 1) {
				distanceMap[x] = 0;
				//real2dImageData[x] = 0;
				artificial2dImage[x] = 0;
			}
		}
	}

	outPath = outputPath + "ImagePlusPath.raw";
	//store2dRawData<dataType>(real2dImageData, Height, Width, outPath.c_str());
	store2dRawData<dataType>(artificial2dImage, Height, Width, outPath.c_str());

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

	////free memory
	//for (k = 0; k < Height; k++) {
	//	free(imageData[k]); free(image[k]); free(potential3D[k]);
	//}
	//free(imageData); free(image); free(potential3D);

	//free(real2dImageData); free(loadingPtr); free(maskThresh);

	//free(artificial2dImage);
	free(distanceMap); free(potentialPtr); free(pathPtr);
	free(startAndEnd);

	//free(seed);

	return EXIT_SUCCESS;
}
