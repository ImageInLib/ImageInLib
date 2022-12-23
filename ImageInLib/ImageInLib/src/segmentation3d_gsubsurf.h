#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

#include "common_functions.h"
#include "heat_equation.h"
#include "filter_params.h"

	bool generalizedSubsurfSegmentation(Image_Data inputImageData, dataType** segFunct, Segmentation_Parameters segParameters, Filter_Parameters explicit_lhe_Parameters,
		Point3D* centers, size_t no_of_centers, unsigned char* outputPathPtr, dataType a, dataType b);
	bool generalizedGFunctionForImageToBeSegmented(Image_Data inputImageData, dataType** extendedCoefPtr, dataType** edgeGradientPtr, Gradient_Pointers GPtrs, Gradient_Pointers VPtrs,
		Segmentation_Parameters segParameters, Filter_Parameters explicit_lhe_Parameters);
	bool generalizedGaussSeidelCoefficients(dataType** extendedCoefPtr, Segment_Image_Data inputImageData, dataType** edgeGradientPtr,
		Coefficient_Pointers CoefPtrs, Segmentation_Parameters segParameters);
	bool generalizedSubsurfSegmentationTimeStep(dataType** prevSol_extPtr, dataType** gauss_seidelPtr, Segment_Image_Data inputImageData, Gradient_Pointers GPtrs, Gradient_Pointers VPtrs,
		Segmentation_Parameters segParameters, Coefficient_Pointers CoefPtrs, Point3D* centers, size_t no_of_centers, dataType a, dataType b);

#ifdef __cplusplus
}
#endif