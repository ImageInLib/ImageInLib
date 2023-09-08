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

	/// <summary>The function perform segmentation by generalized subjective surface method (GSUBSURF). K.Mikula, M.Remesikova, Finite volume schemes for the generalized subjective surface equation in image segmentation, Kybernetika, Vol. 45, No. 4 (2009) pp. 646-656</summary>
	/// <param name="inputImageData">Structure for the image to be segmented</param>
	/// <param name="initialSegment">Initial segment</param>
	/// <param name="segParameters">Segmentation parameters</param>
	/// <param name="implicit_lhe_Parameters">Smoothing parameters</param>
	/// <param name="centers">centers to be considered for the segmentation</param>
	/// <param name="no_of_centers">number of centers in the segment</param>
	/// <param name="outputPathPtr">segmentation path to save the results</param>
	/// <returns></returns>
	bool generalizedSubsurfSegmentation(Image_Data inputImageData, dataType** initialSegment, Segmentation_Parameters segParameters, FilterParameters implicit_lhe_Parameters,
		Point3D* centers, size_t no_of_centers, unsigned char* outputPathPtr);

	bool generalizedGFunctionForImageToBeSegmented(Image_Data inputImageData, dataType** edgeGradientPtr, Gradient_Pointers VPtrs,
		Segmentation_Parameters segParameters, FilterParameters implicit_lhe_Parameters, dataType coef_conv);

	/// <summary>
	/// 
	/// </summary>
	/// <param name="inputImageData"></param>
	/// <param name="edgeGradientPtr"></param>
	/// <param name="CoefPtrs"></param>
	/// <param name="VPtrs"></param>
	/// <param name="segParameters"></param>
	/// <param name="coef_dif"></param>
	/// <returns></returns>
	bool generalizedGaussSeidelCoefficients(Segment_Image_Data inputImageData, dataType** edgeGradientPtr, Coefficient_Pointers CoefPtrs, 
		Gradient_Pointers VPtrs, Segmentation_Parameters segParameters, dataType coef_dif);

	/// <summary>
	/// 
	/// </summary>
	/// <param name="prevSol_extPtr"></param>
	/// <param name="gauss_seidelPtr"></param>
	/// <param name="inputImageData"></param>
	/// <param name="segParameters"></param>
	/// <param name="CoefPtrs"></param>
	/// <param name="centers"></param>
	/// <param name="no_of_centers"></param>
	/// <returns></returns>
	bool generalizedSubsurfSegmentationTimeStep(dataType** prevSol_extPtr, dataType** gauss_seidelPtr, Segment_Image_Data inputImageData,
		Segmentation_Parameters segParameters, Coefficient_Pointers CoefPtrs, Point3D* centers, size_t no_of_centers);

#ifdef __cplusplus
}
#endif