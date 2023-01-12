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

#define originalMean 71.1245
#define offSet 1024
#define standarDeviation 22.001
#define thresmin 995
#define thresmax 1213
#define minimalSize 2000



int main() {

	size_t i, j, k;

	//image Dimensions
	const size_t Width = 512;
	const size_t Length = 512;
	const size_t Height = 406;
	const size_t dim2D = Width * Length;
	/*The Width and Lenght are the same for all the patients
	Height : patient1b = patient3 = 508, patient2 = 406, patient6 = 880, patient7 = 607*/
	
	//Preparation of image data pointers
	dataType** imageData = (dataType**)malloc(Height * sizeof(dataType*));
	short** image = (short**)malloc(Height * sizeof(short*));
	dataType** liverContainer = (dataType**)malloc(Height * sizeof(dataType*));
	short** liver = (short**)malloc(Height * sizeof(short*));
	for (k = 0; k < Height; k++) {
		imageData[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		image[k] = (short*)malloc(dim2D * sizeof(short));
		liverContainer[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		liver[k] = (short*)malloc(dim2D * sizeof(short));
	}

	if (imageData == NULL || image == NULL || liverContainer == NULL || liver == NULL) 
		return false;

	std::string inputPath = "input/";	// "C:/Users/Konan Allaly/Documents/Tests/input/";
	std::string outputPath = "output/";	// "C:/Users/Konan Allaly/Documents/Tests/output/";
	std::string inputImagePath = inputPath + "patient2.raw";
	std::string inputShapePath = inputPath + "liver_p2.raw";

	if (load3dArrayRAW<short>(image, Length, Width, Height, inputImagePath.c_str()) == false)
	{
		printf("inputImagePath does not exist\n");
	}

	if (load3dArrayRAW<short>(liver, Length, Width, Height /*zDim*/, inputShapePath.c_str()) == false)
	{
		printf("inputShapePath does not exist\n");
	}

	//Copy
	//Original data type is short, so I load it with short pointer and copy in dataType=float pointer

	size_t x;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				x = x_new(i, j, Length);
				imageData[k][x_new(i, j, Length)] = (dataType)image[k][x];
				liverContainer[k][x_new(i, j, Length)] = (dataType)liver[k][x];
			}
		}
	}

	//Vtk_File_Info * savingInfo = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	//savingInfo->spacing[0] = 1.171875; savingInfo->spacing[1] = 1.171875; savingInfo->spacing[2] = 1.171875;
	//savingInfo->origin[0] = 0.0; savingInfo->origin[1] = 0.0; savingInfo->origin[2] = 0.0;
	//savingInfo->dimensions[0] = Length; savingInfo->dimensions[1] = Width; savingInfo->dimensions[2] = Height;
	//savingInfo->vDataType = dta_Flt; savingInfo->operation = copyTo;
	//vtkDataForm dataForm = dta_binary;
	//savingInfo->dataPointer = imageData;

	////std::string outputPathVTK = outputPath + "loadedImage.vtk";
	//const char* pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/loadedImage.vtk";
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);
	//pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/loadedLiver.vtk";
	//savingInfo->dataPointer = liverContainer;
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);

	//--------------------------------------------------------------------------------------------------
	//Image Interpolation
	dataType k_spacingNew = 1.171875, k_spacingOld = 2.5;
	//const size_t zDim = (size_t)((k_spacingOld / k_spacingNew) * Height);
	/*spacing patient1b, patient2, patient3 , x = y = 1.171875, z = 2.5
	spacing patient 6, patient7, x = y = 0.9765625*/
	const size_t zDim = 866;
	//Interpolated Image container
	dataType** resampledImage = (dataType**)malloc(zDim * sizeof(dataType*));
	dataType** resampledLiver = (dataType**)malloc(zDim * sizeof(dataType*));
	for (k = 0; k < zDim; k++) {
		resampledImage[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		resampledLiver[k] = (dataType*)malloc(dim2D * sizeof(dataType));
	}
	if (resampledImage == NULL || resampledLiver == NULL) return false;

	linear2dInterpolation(imageData, resampledImage, Length, Width, Height, k_spacingOld, k_spacingNew);
	linear2dInterpolation(liverContainer, resampledLiver, Length, Width, Height, k_spacingOld, k_spacingNew);

	//savingInfo->dimensions[2] = zDim;
	//pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/resampledImage.vtk";
	//savingInfo->dataPointer = resampledImage;
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);
	//pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/resampledLiver.vtk";
	//savingInfo->dataPointer = resampledLiver;
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);
	 
	//-----------------------------------------------------------------------------------------------------

	////Save as .raw file
	//std::string outputInterpolatedImagePath = outputPath + "interpolatedImage_LI.raw";
	//store3dRawData<dataType>(resampledImage, Length, Width, zDim, outputInterpolatedImagePath.c_str());

	////Save as .vtk file
	//Vtk_File_Info * savingInfo = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	//savingInfo->spacing[0] = 1.0; savingInfo->spacing[1] = 1.0; savingInfo->spacing[2] = 1.0;
	//savingInfo->origin[0] = 0.0; savingInfo->origin[1] = 0.0; savingInfo->origin[2] = 0.0;
	//savingInfo->dimensions[0] = Length; savingInfo->dimensions[1] = Width; savingInfo->dimensions[2] = zDim;
	//savingInfo->vDataType = dta_Flt; savingInfo->operation = copyTo;
	//vtkDataForm dataForm = dta_binary;
	//savingInfo->dataPointer = resampledImageData;
	//const char* pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/resampledImage2.vtk";
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);
	//------------------------------------------------------------------------------------------------------

	//Manual croppping
	 
	// Stating points for cropping with original data
	//size_t kMin = 155, iMin = 180, jMin = 126, kn, in, jn;
	//const size_t heightNew = 80, lengthNew = 180, widthNew = 180;

	//kMin depends on the interpolated image
	size_t kMin = 332, iMin = 126, jMin = 180, kn, in, jn;
	const size_t heightNew = 170, lengthNew = 180, widthNew = 180, dim2dNew = lengthNew * widthNew;

	//size_t kMin = 332, iMin = 180, jMin = 126, kn, in, jn;
	//const size_t heightNew = 170, lengthNew = 180, widthNew = 180, dim2dNew = lengthNew * widthNew;

	dataType** croppedImage = (dataType**)malloc(sizeof(dataType*) * heightNew);
	dataType** croppedLiver = (dataType**)malloc(sizeof(dataType*) * heightNew);
	dataType** maskThreshold = (dataType**)malloc(sizeof(dataType*) * heightNew);
	dataType** distanceMap = (dataType**)malloc(sizeof(dataType*) * heightNew);
	dataType** initialSegment = (dataType**)malloc(sizeof(dataType*) * heightNew);
	for (k = 0; k < heightNew; k++) {
		croppedImage[k] = (dataType*)malloc(sizeof(dataType) * dim2dNew);
		croppedLiver[k] = (dataType*)malloc(sizeof(dataType) * dim2dNew);
		distanceMap[k] = (dataType*)malloc(sizeof(dataType) * dim2dNew);
		maskThreshold[k] = (dataType*)malloc(sizeof(dataType) * dim2dNew);
		initialSegment[k] = (dataType*)malloc(sizeof(dataType) * dim2dNew);
	}
	if (croppedImage == NULL || croppedLiver == NULL || distanceMap == NULL || maskThreshold == NULL || initialSegment == NULL) return false;

	const dataType initialSegmentationValue = 1.0;

	for (k = 0, kn = kMin; k < heightNew; k++, kn++) {
		for (i = 0, in = iMin; i < lengthNew; i++, in++) {
			for (j = 0, jn = jMin; j < widthNew; j++, jn++) {
				croppedImage[k][x_new(i, j, lengthNew)] = resampledImage[kn][x_new(in, jn, Length)];
				maskThreshold[k][x_new(i, j, lengthNew)] = croppedImage[k][x_new(i, j, lengthNew)];
				croppedLiver[k][x_new(i, j, lengthNew)] = resampledLiver[kn][x_new(in, jn, Length)];
				initialSegment[k][x_new(i, j, lengthNew)] = 0;
			}
		}
	}

	//Vtk_File_Info* savingInfo = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	//savingInfo->spacing[0] = 1.171875; savingInfo->spacing[1] = 1.171875; savingInfo->spacing[2] = 1.171875;
	//savingInfo->origin[0] = 0.0; savingInfo->origin[1] = 0.0; savingInfo->origin[2] = 0.0;
	//savingInfo->dimensions[0] = lengthNew; savingInfo->dimensions[1] = widthNew; savingInfo->dimensions[2] = heightNew;
	//savingInfo->vDataType = dta_Flt; savingInfo->operation = copyTo;
	//vtkDataForm dataForm = dta_binary;
	//savingInfo->dataPointer = croppedImage;
	//const char* pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/croppedImage.vtk";
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);
	//pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/croppedLiver.vtk";
	//savingInfo->dataPointer = croppedLiver;
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);

	//------------------------------------------------------------------------------------------------------------

	//Fast sweeping to find the point with the higest distance
	thresholding3dFunctionN(maskThreshold, lengthNew, widthNew, heightNew, thresmin, thresmax, 0, 1);
	fastSweepingFunction_3D(distanceMap, maskThreshold, lengthNew, widthNew, heightNew, 1, 100000000, 0);
	dataType distanceMax = -1;
	for (k = 0; k < heightNew; k++) {
		for (i = 0; i < lengthNew; i++) {
			for (j = 0; j < widthNew; j++) {
				if (distanceMap[k][x_new(i, j, lengthNew)] >= distanceMax) {
					distanceMax = distanceMap[k][x_new(i, j, lengthNew)];
				}
			}
		}
	}
	int i_max = 0, j_max = 0, k_max = 0;
	for (k = 0; k < heightNew; k++) {
		for (i = 0; i < lengthNew; i++) {
			for (j = 0; j < widthNew; j++) {
				if (distanceMap[k][x_new(i, j, lengthNew)] == distanceMax) {
					i_max = (int)i; j_max = (int)j; k_max = (int)k;
				}
			}
		}
	}
	printf("Maximal distance : %f \n", distanceMax);
	printf("Coordinates of the highest distance x = %d, y = %d and z = %d \n", i_max, j_max, k_max);

	//dataType diff = distanceMax - 10, distanceMin = distanceMax;
	//int i_min = 0, j_min = 0, k_min = 0;
	//for (k = 0; k < heightNew; k++) {
	//	for (i = 0; i < lengthNew; i++) {
	//		for (j = 0; j < widthNew; j++) {
	//			if (distanceMap[k][x_new(i, j, lengthNew)] <= diff) {
	//				distanceMap[k][x_new(i, j, lengthNew)] = 0;
	//			}
	//			if (distanceMap[k][x_new(i, j, lengthNew)] != 0 && distanceMap[k][x_new(i, j, lengthNew)] <= distanceMin) {
	//				distanceMin = distanceMap[k][x_new(i, j, lengthNew)];
	//				i_min = (int)i; j_min = (int)j; k_min = (int)k;
	//			}
	//		}
	//	}
	//}
	//printf("Second distance : %f \n", distanceMin);
	//printf("Coordinates of the second Point x = %d, y = %d and z = %d \n", i_min, j_min, k_min);

	//pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/thresholded.vtk";
	//savingInfo->dataPointer = maskThreshold;
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);
	//pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/distanceMap.vtk";
	//savingInfo->dataPointer = distanceMap;
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);

	////------------------------------------------------------------------------------------------------------------
	////Labelling
	//int** labelArray = (int**)malloc(sizeof(int*) * heightNew);
	//bool** statusArray = (bool**)malloc(sizeof(bool*) * heightNew);
	//for (k = 0; k < heightNew; k++) {
	//	labelArray[k] = (int*)malloc(sizeof(int) * lengthNew * widthNew);
	//	statusArray[k] = (bool*)malloc(sizeof(bool) * lengthNew * widthNew);
	//}
	//if (labelArray == NULL || statusArray == NULL) return false;
	////Initialization
	//for (k = 0; k < heightNew; k++) {
	//	for (i = 0; i < lengthNew; i++) {
	//		for (j = 0; j < widthNew; j++) {
	//			labelArray[k][x_new(i, j, lengthNew)] = 0;
	//			statusArray[k][x_new(i, j, lengthNew)] = false;
	//		}
	//	}
	//}
	//labelling3D(maskThreshold, labelArray, statusArray, lengthNew, widthNew, heightNew, 1);
	////Number of object voxels
	//int objectVoxels = 0;
	//for (k = 0; k < heightNew; k++) {
	//	for (i = 0; i < lengthNew; i++) {
	//		for (j = 0; j < widthNew; j++) {
	//			if (maskThreshold[k][x_new(i, j, lengthNew)] == 1) {
	//				objectVoxels++;
	//			}
	//		}
	//	}
	//}
	////Counting
	//int* countingArray = (int*)malloc(sizeof(int) * objectVoxels);
	//if (countingArray == NULL) return false;
	//for (k = 0; k < heightNew; k++) {
	//	for (i = 0; i < lengthNew; i++) {
	//		for (j = 0; j < widthNew; j++) {
	//			if (labelArray[k][x_new(i, j, lengthNew)] > 0) {
	//				countingArray[labelArray[k][x_new(i,j,lengthNew)]]++;
	//			}
	//		}
	//	}
	//}
	////biggest region
	//int maxElement = 0;
	//for (i = 0; i < objectVoxels; i++) {
	//	if (countingArray[i] > maxElement) {
	//		maxElement = countingArray[i];
	//	}
	//}
	////Saving the biggest element
	//for (k = 0; k < heightNew; k++) {
	//	for (i = 0; i < lengthNew; i++) {
	//		for (j = 0; j < widthNew; j++) {
	//			if (countingArray[labelArray[k][x_new(i, j, lengthNew)]] != maxElement) {
	//				maskThreshold[k][x_new(i, j, lengthNew)] = 0;
	//			}
	//		}
	//	}
	//}

	//----------------------------------------------------------------------------------------------------
	//Saving
	 
	//rescaleNewRange(initialSegment, lengthNew, widthNew, heightNew, 0, 1);
	//store3dRawData<dataType>(initialSegment, lengthNew, widthNew, heightNew, "C:/Users/Konan Allaly/Documents/Tests/output/initialSegRescall.raw");

	//Vtk_File_Info * savingInfo = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	//savingInfo->spacing[0] = k_spacingNew; savingInfo->spacing[1] = k_spacingNew; savingInfo->spacing[2] = k_spacingNew;
	//savingInfo->origin[0] = 0.0; savingInfo->origin[1] = 0.0; savingInfo->origin[2] = 0.0;
	//savingInfo->dimensions[0] = lengthNew; savingInfo->dimensions[1] = widthNew; savingInfo->dimensions[2] = heightNew;
	//savingInfo->vDataType = dta_Flt; savingInfo->operation = copyTo;
	//vtkDataForm dataForm = dta_binary;
	//savingInfo->dataPointer = initialSegment;
	//const char* pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/_seg_func_000.vtk";
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);

	//savingInfo->dataPointer = croppedLiver;
	//pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/initialSeg.vtk";
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);

	//savingInfo->dataPointer = croppedImage;
	//pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/croppedVolume.vtk";
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);

	//----------------------------------------------------------------------------------------------------
	//Artificial Image

	////Create artificial image : if we are working with artificial image , everything befor this part is not needed
	//const size_t heightNew = 64, lengthNew = 150, widthNew = 150;
	//size_t i_max = lengthNew / 2, j_max = widthNew / 2, k_max = heightNew / 2;
	//size_t i_n = i_max - 27, j_n = j_max, k_n = k_max;
	//size_t i_m = i_max + 27, j_m = j_max, k_m = k_max;
	//dataType** Image = (dataType**)malloc(sizeof(dataType*) * heightNew);
	//dataType** initialSegment = (dataType**)malloc(sizeof(dataType*) * heightNew);
	//for (k = 0; k < heightNew; k++) {
	//	Image[k] = (dataType*)malloc(sizeof(dataType) * lengthNew * widthNew);
	//  initialSegment[k] = (dataType*)malloc(sizeof(dataType) * lengthNew * widthNew);
	//}
	//if (Image == NULL || initialSegment == NULL) return false;

	//// Initialization
	//for (k = 0; k < heightNew; k++) {
	//	for (i = 0; i < lengthNew; i++) {
	//		for (j = 0; j < widthNew; j++) {
	//			if (sqrt((k - k_n) * (k - k_n) + (i - i_n) * (i - i_n) + (j - j_n) * (j - j_n)) <= 30) {
	//				Image[k][x_new(i, j, lengthNew)] = 1;
	//			}
	//			if (sqrt((k - k_m) * (k - k_m) + (i - i_m) * (i - i_m) + (j - j_m) * (j - j_m)) <= 30) {
	//				Image[k][x_new(i, j, lengthNew)] = 1;
	//			}
	//		}
	//	}
	//}

	////add salt and pepper noise
	//saltAndPepper3dNoise_D(Image, lengthNew, widthNew, heightNew, 0.5, 1);

    //-----------------------------------------------------------------------------------------------------------
	//Segmentation parameters for artificial image
	//size_t numb_centers = 1; Point3D* centerSeg = (Point3D*)malloc(sizeof(Point3D) * numb_centers);
	//centerSeg->x = i_max; centerSeg->y = j_max; centerSeg->z = k_max; //---> used for one center
	//// used for multiple centers
	////centerSeg[0].x = i_n; centerSeg[0].y = j_max; centerSeg[0].z = k_max;
	////centerSeg[1].x = i_m; centerSeg[1].y = j_max; centerSeg[1].z = k_max;
	//generateInitialSegmentationFunctionForMultipleCentres(initialSegment, lengthNew, widthNew, heightNew, centerSeg, 0.5, 60, numb_centers);
	//Image_Data segment; segment.height = heightNew; segment.length = lengthNew; segment.width = widthNew; segment.imageDataPtr = Image;
	//Segmentation_Parameters segmentParameters; segmentParameters.coef = 10; segmentParameters.eps2 = 1e-6; segmentParameters.gauss_seidelTolerance = 1e-3;
	//segmentParameters.h = 1.171875; segmentParameters.maxNoGSIteration = 100; segmentParameters.maxNoOfTimeSteps = 500; segmentParameters.mod = 2;
	//segmentParameters.numberOfTimeStep = 500; segmentParameters.omega_c = 1.5; segmentParameters.segTolerance = 1e-4; segmentParameters.tau = 4;
	//Filter_Parameters filterParameters; filterParameters.coef = 1e-6; filterParameters.edge_detector_coefficient = 100; filterParameters.eps2 = 1e-6;
	//filterParameters.h = 1.171875; filterParameters.maxNumberOfSolverIteration = 100; filterParameters.omega_c = 1.1; filterParameters.p = 1;
	//filterParameters.sigma = 1e-3; filterParameters.timeStepSize = 1.2; filterParameters.timeStepsNum = 1; filterParameters.tolerance = 1e-3;
	//unsigned char outputPathPtr[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/";
	//subsurfSegmentation(segment, initialSegment, segmentParameters, filterParameters, centerSeg, numb_centers, outputPathPtr);
	//generalizedSubsurfSegmentation(segment, initialSegment, segmentParameters, filterParameters, centerSeg, numb_centers, outputPathPtr, 1.0, 1.0);
	//--------------------------------------------------------------------------------------------------------------

	//------------------------------------------------------------------------------------------------------
	//Segmentation parameters for real image

	size_t numb_centers = 2; Point3D* centerSeg = (Point3D*)malloc(sizeof(Point3D) * numb_centers);
	//centerSeg->x = i_max; centerSeg->y = j_max; centerSeg->z = k_max; //---> used for one center
	////used for multiple centers
	centerSeg[0].x = i_max; centerSeg[0].y = j_max; centerSeg[0].z = k_max;
	centerSeg[1].x = i_max + 50; centerSeg[1].y = j_max - 70; centerSeg[1].z = k_max;
	//centerSeg[1].x = i_max - 70; centerSeg[1].y = j_max + 50; centerSeg[1].z = k_max;

	//If we want to start with the segmentatation function originally implemented in the library
	generateInitialSegmentationFunctionForMultipleCentres(initialSegment, lengthNew, widthNew, heightNew, centerSeg, 0.5, 15, numb_centers);
	//for (k = 0; k < heightNew; k++) {
	//	for (i = 0; i < lengthNew; i++) {
	//		for (j = 0; j < widthNew; j++) {
	//			if (croppedLiver[k][x_new(i, j, lengthNew)] != 0) {
	//				initialSegment[k][x_new(i, j, lengthNew)] = 1.0;
	//			}
	//		}
	//	}
	//}

	//Save the initial segmentation function
	Vtk_File_Info* savingInfo = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	savingInfo->spacing[0] = k_spacingNew; savingInfo->spacing[1] = k_spacingNew; savingInfo->spacing[2] = k_spacingNew;
	savingInfo->origin[0] = 0.0; savingInfo->origin[1] = 0.0; savingInfo->origin[2] = 0.0;
	savingInfo->dimensions[0] = lengthNew; savingInfo->dimensions[1] = widthNew; savingInfo->dimensions[2] = heightNew;
	savingInfo->vDataType = dta_Flt; savingInfo->operation = copyTo;
	vtkDataForm dataForm = dta_binary;
	savingInfo->dataPointer = initialSegment;
	const char* pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/_seg_func_000.vtk";
	storeVtkFile(pathsaveVTK, savingInfo, dataForm);

	//If we want start by the liver model just comment the previous line
	Image_Data segment; segment.height = heightNew; segment.length = lengthNew; segment.width = widthNew; segment.imageDataPtr = croppedImage;
	rescaleNewRange(segment.imageDataPtr, lengthNew, widthNew, heightNew, 0, 1);
	Segmentation_Parameters segmentParameters; segmentParameters.coef = 10000; segmentParameters.eps2 = 1e-6; segmentParameters.gauss_seidelTolerance = 1e-3;
	segmentParameters.h = k_spacingNew; segmentParameters.maxNoGSIteration = 100; segmentParameters.maxNoOfTimeSteps = 5000; segmentParameters.mod = 10;
	segmentParameters.numberOfTimeStep = 5000; segmentParameters.omega_c = 1.5; segmentParameters.segTolerance = 1e-4; segmentParameters.tau = 4;
	Filter_Parameters filterParameters; filterParameters.coef = 1e-6; filterParameters.edge_detector_coefficient = 100; filterParameters.eps2 = 1e-6;
	filterParameters.h = k_spacingNew; filterParameters.maxNumberOfSolverIteration = 100; filterParameters.omega_c = 1.5; filterParameters.p = 1;
	filterParameters.sigma = 1e-3; filterParameters.timeStepSize = 1.2; filterParameters.timeStepsNum = 1; filterParameters.tolerance = 1e-3;

	unsigned char outputPathPtr[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/";
	//subsurfSegmentation(segment, initialSegment, segmentParameters, filterParameters, centerSeg, numb_centers, outputPathPtr);
	generalizedSubsurfSegmentation(segment, initialSegment, segmentParameters, filterParameters, centerSeg, numb_centers, outputPathPtr, 1.0, 1.0);

	//------------------------------------------------------------------------------------------------

	//free memory
	for (k = 0; k < Height; k++) {
		free(imageData[k]); free(image[k]);
		free(liverContainer[k]); free(liver[k]);
	}
	free(imageData); free(image); free(liverContainer); free(liver);

	for (k = 0; k < zDim; k++) {
		free(resampledImage[k]); free(resampledLiver[k]);
	}
	free(resampledImage); free(resampledLiver);

	for (k = 0; k < heightNew; k++) {
		free(croppedImage[k]); free(croppedLiver[k]); free(distanceMap[k]); free(maskThreshold[k]); free(initialSegment[k]);
		//free(labelArray[k]); free(statusArray[k]);
	}
	free(croppedImage); free(croppedLiver); free(distanceMap); free(maskThreshold); free(initialSegment);
	//free(labelArray); free(statusArray);
	 
	free(centerSeg);
	free(savingInfo);
	//free(countingArray);

	return EXIT_SUCCESS;
}
