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
#include "../src/image_difference.h"
#include "../src/noise_generator.h"
#include "../src/image_norm.h"
#include "resamplingVolume.h"

#define originalMean 71.1245
#define offSet 1024
#define standarDeviation 22.001
#define thresmin 995 // 990
#define thresmax 1213 // 1250
#define minimalSize 2000


int main() {

	size_t i, j, k;
	const size_t Width = 512;
	const size_t Length = 512;
	const size_t Height = 406;
	const size_t dim2D = Width * Length;

	dataType** imageData = (dataType**)malloc(Height * sizeof(dataType*));
	short** image = (short**)malloc(Height * sizeof(short*));
	for (k = 0; k < Height; k++) {
		imageData[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		image[k] = (short*)malloc(dim2D * sizeof(short));
	}
	if (imageData == NULL || image == NULL) return false;
	const char* pathLoad = "C:/Users/Konan Allaly/Documents/Tests/input/patient2.raw";
	load3dArrayRAW<short>(image, Length, Width, Height, pathLoad);
	//Copy
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				imageData[k][x_new(i, j, Length)] = (dataType)image[k][x_new(i, j, Length)];
			}
		}
	}

	//store3dRawData<dataType>(imageData, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/loadedP2.raw");
	//store3dRawData<dataType>(edgePtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/loadPatient2Liver.raw");

	dataType k_spacingNew = 1.171875, k_spacingOld = 2.5;
	//const size_t zDim = (size_t)((k_spacingOld / k_spacingNew) * Height);
	const size_t zDim = 866;

	dataType** resampledImageData = (dataType**)malloc(zDim * sizeof(dataType*));
	dataType** liverContainer = (dataType**)malloc(zDim * sizeof(dataType*));
	short** liver = (short**)malloc(zDim * sizeof(short*));
	for (k = 0; k < zDim; k++) {
		resampledImageData[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		liverContainer[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		liver[k] = (short*)malloc(dim2D * sizeof(short));
	}
	if (resampledImageData == NULL || liverContainer == NULL || liver == NULL) return false;

	pathLoad = "C:/Users/Konan Allaly/Documents/Tests/input/liverP222.raw";
	load3dArrayRAW<short>(liver, Length, Width, zDim, pathLoad);
	for (k = 0; k < zDim; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				liverContainer[k][x_new(i, j, Length)] = (dataType)liver[k][x_new(i, j, Length)];
			}
		}
	}

	//linear2dInterpolation(imageData, resampledImageData, Length, Width, Height, k_spacingOld, k_spacingNew);
	//nearestNeighborInterpolation(imageData, resampledImageData, Length, Width, Height, k_spacingOld, k_spacingNew);
	//store3dRawData<dataType>(resampledImageData, Length, Width, zDim, "C:/Users/Konan Allaly/Documents/Tests/output/interpolatedImage_LI.raw");

	//store3dRawData<dataType>(resampledImageData, Length, Width, zDim, "C:/Users/Konan Allaly/Documents/Tests/output/resampled2.raw");
	//Vtk_File_Info * savingInfo = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	//savingInfo->spacing[0] = 1.17; savingInfo->spacing[1] = 1.17; savingInfo->spacing[2] = 1.17;
	//savingInfo->origin[0] = 0.0; savingInfo->origin[1] = 0.0; savingInfo->origin[2] = 0.0;
	//savingInfo->dimensions[0] = Length; savingInfo->dimensions[1] = Width; savingInfo->dimensions[2] = zDim;
	//savingInfo->vDataType = dta_Flt; savingInfo->operation = copyTo;
	//vtkDataForm dataForm = dta_binary;
	//savingInfo->dataPointer = resampledImageData;
	//const char* pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/resampledImage2.vtk";
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);
	////thresholding3dFunctionN(imageData, Length, Width, Height, thresmin, thresmax, 0, 1);

	//fastSweepingFunction_3D(edgePtr, imageData, Length, Width, Height, 1, 100000000, 0);
	//dataType distanceMax = -1;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (edgePtr[k][x_new(i, j, Length)] >= distanceMax) {
	//				distanceMax = edgePtr[k][x_new(i, j, Length)];
	//			}
	//		}
	//	}
	//}
	//dataType cptMax = 0;
	//int i_max = 0, j_max = 0, k_max = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (edgePtr[k][x_new(i, j, Length)] == distanceMax) {
	//				cptMax++;
	//				i_max = (int)i; j_max = (int)j; k_max = (int)k;
	//			}
	//		}
	//	}
	//}
	//printf("Maximal distance : %f \n", distanceMax);
	//printf("Coordinates of the highest distance x = %d, y = %d and z = %d \n", i_max, j_max, k_max);

	//Vtk_File_Info * savingInfo = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	//savingInfo->spacing[0] = 1.0; savingInfo->spacing[1] = 1.0; savingInfo->spacing[2] = 1.0;
	//savingInfo->origin[0] = 0.0; savingInfo->origin[1] = 0.0; savingInfo->origin[2] = 0.0;
	//savingInfo->dimensions[0] = Length; savingInfo->dimensions[1] = Width; savingInfo->dimensions[2] = Height;
	//savingInfo->vDataType = dta_Flt; savingInfo->operation = copyTo;
	//vtkDataForm dataForm = dta_binary;
	//savingInfo->dataPointer = imageData;
	//const char* pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/loaded.vtk";
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);
	////savingInfo->dataPointer = edgePtr;
	////pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/distanceMap7.vtk";
	////storeVtkFile(pathsaveVTK, savingInfo, dataForm);

	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if ( sqrt((i - i_max)* (i - i_max) + (j - j_max) * (j - j_max) + (k - k_max) * (k - k_max)) < 10 ) {
	//				imageData[k][x_new(i, j, Length)] = 1;
	//			}
	//			else {
	//				imageData[k][x_new(i, j, Length)] = 0;
	//			}
	//		}
	//	}
	//}
	//savingInfo->dataPointer = imageData;
	//pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/ballCenter7.vtk";
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);

	////K-mean equalization
	//size_t K = 100;
	//dataType minData = 1000000, maxData = -1000000;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (imageData[k][x_new(i, j, Length)] < minData) {
	//				minData = imageData[k][x_new(i, j, Length)];
	//			}
	//			if (imageData[k][x_new(i, j, Length)] > maxData) {
	//				maxData = imageData[k][x_new(i, j, Length)];
	//			}
	//		}
	//	}
	//}
	//printf("min data = %f and max data = %f \n", minData, maxData);
	//dataType sizeOfGroup = (maxData - minData) / K, dataRangeMin = minData, mean[100], dataRangeMax;
	//size_t l, cpt, sum;
	//for (l = 0; l < K; l++) {
	//	cpt = 0; sum = 0; dataRangeMax = dataRangeMin + sizeOfGroup;
	//	for (k = 0; k < Height; k++) {
	//		for (i = 0; i < Length; i++) {
	//			for (j = 0; j < Width; j++) {
	//				if (imageData[k][x_new(i, j, Length)] >= dataRangeMin && imageData[k][x_new(i, j, Length)] <= dataRangeMax) {
	//					sum = sum + imageData[k][x_new(i, j, Length)];
	//					cpt++;
	//				}
	//			}
	//		}
	//	}
	//	mean[l] = sum / cpt;
	//	for (k = 0; k < Height; k++) {
	//		for (i = 0; i < Length; i++) {
	//			for (j = 0; j < Width; j++) {
	//				if (imageData[k][x_new(i, j, Length)] >= dataRangeMin && imageData[k][x_new(i, j, Length)] <= dataRangeMax) {
	//					imageData[k][x_new(i, j, Length)] = mean[l];
	//				}
	//			}
	//		}
	//	}
	//	dataRangeMin = dataRangeMax;
	//}
	////store3dRawData<dataType>(imageData, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/thresh100Mean.raw");

	//rescaleNewRange(imageData, Length, Width, Height, 0, 1);
	//Filter_Parameters smoothPa; smoothPa.h = 1.0; smoothPa.edge_detector_coefficient = 1000; smoothPa.maxNumberOfSolverIteration = 300;
	//smoothPa.omega_c = 1.1; smoothPa.timeStepSize = 1.2; smoothPa.p = 1.0; smoothPa.tolerance = 1e-4;
	//imageGradient(imageData, "C:/Users/Konan Allaly/Documents/Tests/output/", Length, Width, Height, 1.0);
	//edgesDetector(imageData, "C:/Users/Konan Allaly/Documents/Tests/output/", Length, Width, Height, smoothPa);

	//thresholding3dFunctionN(imageData, Length, Width, Height, thresmin, thresmax, 0, 1);
	//edgeDetection3dFunctionD(imageData, edgePtr, Length, Width, Height, 0, 4000, 2000);
	//store3dRawData<dataType>(edgePtr, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/TestEdges2.raw");

	//Manual croppping
	//size_t kMin = 155, iMin = 180, jMin = 126, kn, in, jn;
	//const size_t heightNew = 80, lengthNew = 180, widthNew = 180;

	size_t kMin = 332, iMin = 180, jMin = 126, kn, in, jn;
	const size_t heightNew = 170, lengthNew = 180, widthNew = 180;
	//size_t kMin = 388, iMin = 210, jMin = 148, kn, in, jn;
	//const size_t heightNew = 200, lengthNew = 210, widthNew = 210;
	dataType** croppedImage = (dataType**)malloc(sizeof(dataType*) * heightNew);
	dataType** croppedLiver = (dataType**)malloc(sizeof(dataType*) * heightNew);
	dataType** maskThreshold = (dataType**)malloc(sizeof(dataType*) * heightNew);
	dataType** distanceMap = (dataType**)malloc(sizeof(dataType*) * heightNew);
	dataType** initialSegment = (dataType**)malloc(sizeof(dataType*) * heightNew);
	for (k = 0; k < heightNew; k++) {
		croppedImage[k] = (dataType*)malloc(sizeof(dataType) * lengthNew * widthNew);
		croppedLiver[k] = (dataType*)malloc(sizeof(dataType) * lengthNew * widthNew);
		distanceMap[k] = (dataType*)malloc(sizeof(dataType) * lengthNew * widthNew);
		maskThreshold[k] = (dataType*)malloc(sizeof(dataType) * lengthNew * widthNew);
		initialSegment[k] = (dataType*)malloc(sizeof(dataType) * lengthNew * widthNew);
	}
	if (croppedImage == NULL || croppedLiver == NULL || distanceMap == NULL || maskThreshold == NULL || initialSegment == NULL) return false;
	initialize3dArrayD(croppedImage, lengthNew, widthNew, heightNew, 0);
	initialize3dArrayD(croppedLiver, lengthNew, widthNew, heightNew, 0);
	initialize3dArrayD(distanceMap, lengthNew, widthNew, heightNew, 0);
	initialize3dArrayD(maskThreshold, lengthNew, widthNew, heightNew, 0);
	initialize3dArrayD(initialSegment, lengthNew, widthNew, heightNew, 0);
	for (k = 0, kn = kMin; k < heightNew; k++, kn++) {
		for (i = 0, in = iMin; i < lengthNew; i++, in++) {
			for (j = 0, jn = jMin; j < widthNew; j++, jn++) {
				croppedImage[k][x_new(i, j, lengthNew)] = resampledImageData[kn][x_new(in, jn, Length)];
				maskThreshold[k][x_new(i, j, lengthNew)] = croppedImage[k][x_new(i, j, lengthNew)];
				croppedLiver[k][x_new(i, j, lengthNew)] = liverContainer[kn][x_new(in, jn, Length)];
				if (croppedLiver[k][x_new(i, j, lengthNew)] != 0) {
					initialSegment[k][x_new(i, j, lengthNew)] = croppedImage[k][x_new(i, j, lengthNew)];
					//maskThreshold[k][x_new(i, j, lengthNew)] = croppedImage[k][x_new(i, j, lengthNew)];
				}
			}
		}
	}

	////autre initialisation
	//for (k = 0; k < heightNew; k++) {
	//	for (i = 0; i < lengthNew; i++) {
	//		for (j = 0; j < widthNew; j++) {
	//			croppedImage[k][x_new(i, j, lengthNew)] = maskThreshold[k][x_new(i, j, lengthNew)];
	//		}
	//	}
	//}

	//Vtk_File_Info * savingInfo = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	//savingInfo->spacing[0] = 1.171875; savingInfo->spacing[1] = 1.171875; savingInfo->spacing[2] = 1.171875;
	//savingInfo->origin[0] = 0.0; savingInfo->origin[1] = 0.0; savingInfo->origin[2] = 0.0;
	//savingInfo->dimensions[0] = lengthNew; savingInfo->dimensions[1] = widthNew; savingInfo->dimensions[2] = heightNew;
	//savingInfo->vDataType = dta_Flt; savingInfo->operation = copyTo;
	//vtkDataForm dataForm = dta_binary;
	//savingInfo->dataPointer = croppedImage;
	//const char* pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/croppedImage_Liver.vtk";
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);
	//store3dRawData<dataType>(croppedImage, lengthNew, widthNew, heightNew, "C:/Users/Konan Allaly/Documents/Tests/output/croppedImage_LI.raw");

	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);
	//savingInfo->dataPointer = liverContainer;
	//pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/liverModel.vtk";
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);

	//thresholding3dFunctionN(initialSeg, lengthNew, widthNew, heightNew, thresmin, thresmax, 0, 1);
	//savingInfo->dataPointer = initialSeg;
	//pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/thresholded_LI.vtk";

	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);
	//////store3dRawData<dataType>(croppedImage, lengthNew, widthNew, heightNew, "C:/Users/Konan Allaly/Documents/Tests/output/Volume.raw");
	//////store3dRawData<dataType>(initialSeg, lengthNew, widthNew, heightNew, "C:/Users/Konan Allaly/Documents/Tests/output/Model.raw");

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
	dataType cptMax = 0;
	int i_max = 0, j_max = 0, k_max = 0;
	for (k = 0; k < heightNew; k++) {
		for (i = 0; i < lengthNew; i++) {
			for (j = 0; j < widthNew; j++) {
				if (distanceMap[k][x_new(i, j, lengthNew)] == distanceMax) {
					cptMax++;
					i_max = (int)i; j_max = (int)j; k_max = (int)k;
				}
			}
		}
	}
	printf("Maximal distance : %f \n", distanceMax);
	printf("Coordinates of the highest distance x = %d, y = %d and z = %d \n", i_max, j_max, k_max);

	//savingInfo->dataPointer = distanceMap;
	//pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/distanceMap_LI.vtk";
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);

	////Create artificial image
	//const size_t heightNew = 64, lengthNew = 150, widthNew = 150;
	//size_t i_max = lengthNew / 2, j_max = widthNew / 2, k_max = heightNew / 2;
	//size_t i_n = i_max - 27, j_n = j_max, k_n = k_max;
	//size_t i_m = i_max + 27, j_m = j_max, k_m = k_max;
	//dataType** Image = (dataType**)malloc(sizeof(dataType*) * heightNew);
	//for (k = 0; k < heightNew; k++) {
	//	Image[k] = (dataType*)malloc(sizeof(dataType) * lengthNew * widthNew);
	//}
	//if (Image == NULL) return false;
	//initialize3dArrayD(Image, lengthNew, widthNew, heightNew, 0);

	////centroid Liver
	//dataType * cenTroidLiver = (dataType*)malloc(3 * sizeof(dataType));
	//centroidImage(liverContainer, cenTroidLiver, heightNew, lengthNew, widthNew, 0);
	//dataType k_max = cenTroidLiver[2], i_max = cenTroidLiver[0], j_max = cenTroidLiver[1];

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

	//add salt and pepper noise
	//saltAndPepper3dNoise_D(croppedImage, lengthNew, widthNew, heightNew, 0.5, 1);

	//Vtk_File_Info * savingInfo = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	//savingInfo->spacing[0] = 1.171875; savingInfo->spacing[1] = 1.171875; savingInfo->spacing[2] = 1.171875;
	//savingInfo->origin[0] = 0.0; savingInfo->origin[1] = 0.0; savingInfo->origin[2] = 0.0;
	//savingInfo->dimensions[0] = lengthNew; savingInfo->dimensions[1] = widthNew; savingInfo->dimensions[2] = heightNew;
	//savingInfo->vDataType = dta_Flt; savingInfo->operation = copyTo;
	//vtkDataForm dataForm = dta_binary;
	//savingInfo->dataPointer = croppedImage;
	//const char* pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/croppedVolume.vtk";
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);
	//savingInfo->dataPointer = maskThreshold;
	//pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/thresh.vtk";
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);
	//savingInfo->dataPointer = croppedLiver;
	//pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/liver.vtk";
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);
	//savingInfo->dataPointer = distanceMap;
	//pathsaveVTK = "C:/Users/Konan Allaly/Documents/Tests/output/distanceMap.vtk";
	//storeVtkFile(pathsaveVTK, savingInfo, dataForm);

	size_t numb_centers = 1; Point3D* centerSeg = (Point3D*)malloc(sizeof(Point3D) * numb_centers);
	//centerSeg[0].x = i_n; centerSeg[0].y = j_max; centerSeg[0].z = k_max;
	//centerSeg[1].x = i_m; centerSeg[1].y = j_max; centerSeg[1].z = k_max;
	centerSeg->x = i_max; centerSeg->y = j_max; centerSeg->z = k_max;
	generateInitialSegmentationFunctionForMultipleCentres(initialSegment, lengthNew, widthNew, heightNew, centerSeg, 0.5, 60, numb_centers);

	Image_Data segment; segment.height = heightNew; segment.length = lengthNew; segment.width = widthNew; segment.imageDataPtr = croppedImage;
	rescaleNewRange(segment.imageDataPtr, lengthNew, widthNew, heightNew, 0, 1);
	rescaleNewRange(initialSegment, lengthNew, widthNew, heightNew, 0, 1);
	Segmentation_Parameters segmentParameters; segmentParameters.coef = 200; segmentParameters.eps2 = 1e-6; segmentParameters.gauss_seidelTolerance = 1e-3;
	segmentParameters.h = 1.171875; segmentParameters.maxNoGSIteration = 100; segmentParameters.maxNoOfTimeSteps = 500; segmentParameters.mod = 2;
	segmentParameters.numberOfTimeStep = 500; segmentParameters.omega_c = 1.5; segmentParameters.segTolerance = 1e-4;
	dataType Tmax = segmentParameters.maxNoOfTimeSteps; segmentParameters.tau = 4; //1.171875 * 1.171875 * 1.171875;
	Filter_Parameters filterParameters; filterParameters.coef = 1e-6; filterParameters.edge_detector_coefficient = 100; filterParameters.eps2 = 1e-6;
	filterParameters.h = 1.171875; filterParameters.maxNumberOfSolverIteration = 100; filterParameters.omega_c = 1.1; filterParameters.p = 1;
	filterParameters.sigma = 1e-3; filterParameters.timeStepSize = 1.2; filterParameters.timeStepsNum = 1; filterParameters.tolerance = 1e-3;

	//const FilterMethod Fmethod = MEAN_CURVATURE_FILTER; //NONLINEAR_HEATEQUATION_IMPLICIT;
	//filterImage(segment, filterParameters, Fmethod);
	//store3dRawData<dataType>(segment.imageDataPtr, lengthNew, widthNew, heightNew, "C:/Users/Konan Allaly/Documents/Tests/output/filterdMC.raw");
	//filterImage(model, filterParameters, Fmethod);
	//store3dRawData<dataType>(model.imageDataPtr, lengthNew, widthNew, heightNew, "C:/Users/Konan Allaly/Documents/Tests/output/filterdModel.raw");

	//for (k = 0; k < heightNew; k++) {
	//	for (i = 0; i < lengthNew; i++) {
	//		for (j = 0; j < widthNew; j++) {
	//			if (model.imageDataPtr[k][x_new(i, j, lengthNew)] != 0) {
	//				model.imageDataPtr[k][x_new(i, j, lengthNew)] = segment.imageDataPtr[k][x_new(i, j, lengthNew)];
	//			}
	//		}
	//	}
	//}

	unsigned char outputPathPtr[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/";
	//subsurfSegmentation(segment, initialSegment, segmentParameters, filterParameters, centerSeg, numb_centers, outputPathPtr);
	//generalizedSubsurfSegmentation(segment, initialSegment, segmentParameters, filterParameters, centerSeg, numb_centers, outputPathPtr, 1.0, 1.0);

	//thresholding3dFunctionN(imageData, Length, Width, Height, thresmin, thresmax, minData, maxData);
	//store3dRawData<dataType>(imageData, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/threshWinth100Mean.raw");
	//erosion3dHeighteenNeigbours(imageData, Length, Width, Height, maxData, minData);
	//store3dRawData<dataType>(imageData, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/Eroded.raw");

	////Labelling
	//int** labelArray = (int**)malloc(Height * sizeof(int*));
	//bool** status = (bool**)malloc(Height * sizeof(bool*));
	//for (k = 0; k < Height; k++) {
	//	labelArray[k] = (int*)malloc(dim2D * sizeof(int));
	//	status[k] = (bool*)malloc(dim2D * sizeof(bool));
	//}
	//if (labelArray == NULL || status == NULL) return false;
	////initialization
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			labelArray[k][x_new(i, j, Length)] = 0;
	//			status[k][x_new(i, j, Length)] = false;
	//		}
	//	}
	//}
	//Filter_Parameters smoothPa; smoothPa.h = 1.0; smoothPa.edge_detector_coefficient = 1000; smoothPa.maxNumberOfSolverIteration = 300;
	//smoothPa.omega_c = 1.1; smoothPa.timeStepSize = 1.2; smoothPa.p = 1.0; smoothPa.tolerance = 1e-4;
	//double start = clock();
	////labelling3D(imageData, labelArray, status, Length, Width, Height, maxData);
	//regionGrowing(imageData, labelArray, status, Length, Width, Height, thresmin, thresmax, smoothPa);
	//double finish = clock();
	//printf("Execution time for the new code : %.3lf \n", (finish - start) / CLOCKS_PER_SEC);
	////store3dRawData<int>(labelArray, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/segmented.raw");
	////number of region cells
	//int numberOfRegionsCells = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (imageData[k][x_new(i, j, Length)] >= thresmin && imageData[k][x_new(i, j, Length)] <= thresmax) {
	//				numberOfRegionsCells++;
	//			}
	//		}
	//	}
	//}
	//printf("Number of Regions Cells : %d \n", numberOfRegionsCells);
	////Counting
	//int* countingArray = (int*)malloc(numberOfRegionsCells * sizeof(int));
	//if (countingArray == NULL) return false;
	//for (i = 0; i < numberOfRegionsCells; i++) countingArray[i] = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (labelArray[k][x_new(i, j, Length)] > 0) {
	//				countingArray[labelArray[k][x_new(i, j, Length)]]++;
	//			}
	//		}
	//	}
	//}
	////Remove small regions 
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (countingArray[labelArray[k][x_new(i, j, Length)]] < minimalSize) {
	//				countingArray[labelArray[k][x_new(i, j, Length)]] = 0;
	//			}
	//		}
	//	}
	//}
	////Number of regions
	//int numberOfRegions = 0;
	//for (i = 0; i < numberOfRegionsCells; i++) {
	//	if (countingArray[i] > 0) {
	//		numberOfRegions++;
	//	}
	//}
	//printf("Number of regions = %d \n", numberOfRegions);
	////Save biggest region 
	//int maxVolume = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (countingArray[labelArray[k][x_new(i, j, Length)]] > maxVolume) {
	//				maxVolume = countingArray[labelArray[k][x_new(i, j, Length)]];
	//			}
	//		}
	//	}
	//}
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (countingArray[labelArray[k][x_new(i, j, Length)]] != maxVolume) {
	//				imageData[k][x_new(i, j, Length)] = minData;
	//			}
	//		}
	//	}
	//}
	//store3dRawData<dataType>(imageData, Length, Width, Height, "C:/Users/Konan Allaly/Documents/Tests/output/bigestRegionLabelling.raw");

	//free memory
	for (k = 0; k < Height; k++) {
		free(imageData[k]); free(image[k]);
		//free(labelArray[k]); free(status[k]);
		//free(savingInfo->dataPointer[k]);
	}
	free(imageData); free(image);
	//free(labelArray); free(status);
	//free(savingInfo->dataPointer); 
	//free(countingArray);

	for (k = 0; k < zDim; k++) {
		free(resampledImageData[k]); free(liver[k]); free(liverContainer[k]);
	}
	free(resampledImageData); free(liver); free(liverContainer);

	for (k = 0; k < heightNew; k++) {
		free(croppedImage[k]); free(croppedLiver[k]); free(distanceMap[k]); free(maskThreshold[k]); free(initialSegment[k]);
	}
	free(croppedImage); free(croppedLiver); free(distanceMap); free(maskThreshold); free(initialSegment);

	free(centerSeg);
	//free(savingInfo);

	return EXIT_SUCCESS;
}
