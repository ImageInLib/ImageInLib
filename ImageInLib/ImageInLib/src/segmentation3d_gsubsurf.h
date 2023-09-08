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

	/// <summary>The function performs segmentation by generalized subjective surface method (GSUBSURF). K.Mikula, M.Remesikova, Finite volume schemes for the generalized subjective surface equation in image segmentation, Kybernetika, Vol. 45, No. 4 (2009) pp. 646-656</summary>
	/// <param name="inputImageData"> : Structure for the image to be segmented</param>
	/// <param name="initialSegment"> : Initial segment</param>
	/// <param name="segParameters"> : Segmentation parameters</param>
	/// <param name="implicit_lhe_Parameters"> : Smoothing parameters</param>
	/// <param name="centers"> : centers to be considered for the segmentation</param>
	/// <param name="no_of_centers"> : number of centers considered during the segmentation</param>
	/// <param name="outputPathPtr"> : path to save the regmentation results</param>
	/// <returns></returns>
	bool generalizedSubsurfSegmentation(Image_Data inputImageData, dataType** initialSegment, Segmentation_Parameters segParameters, FilterParameters implicit_lhe_Parameters,
		Point3D* centers, size_t no_of_centers, unsigned char* outputPathPtr);

	/// <summary>
	/// The function performs the smoothing of the image to be segmented and compute the edge detector 
	/// </summary>
	/// <param name="inputImageData">: Structure for the image to be segmented</param>
	/// <param name="edgeGradientPtr"> : pointer for the edge detector</param>
	/// <param name="VPtrs"> : pointers for the gradients of the edge detector</param>
	/// <param name="segParameters">: Segmentation parameters</param>
	/// <param name="implicit_lhe_Parameters">: Smoothing parameters</param>
	/// <param name="coef_conv"> : convection coefficient in GSUBSURF</param>
	/// <returns></returns>
	bool generalizedGFunctionForImageToBeSegmented(Image_Data inputImageData, dataType** edgeGradientPtr, Gradient_Pointers VPtrs,
		Segmentation_Parameters segParameters, FilterParameters implicit_lhe_Parameters, dataType coef_conv);

	/// <summary>
	/// The function perform the gauss seidel itearation to solve the GSUBSURF linear system
	/// </summary>
	/// <param name="inputImageData">: Structure for the image to be segmented</param>
	/// <param name="edgeGradientPtr">: pointer for the edge detector</param>
	/// <param name="CoefPtrs"> : coefficient of the linear system</param>
	/// <param name="VPtrs">: pointers for the gradients of the edge detector</param>
	/// <param name="segParameters">: Segmentation parameters</param>
	/// <param name="coef_dif">: diffusion coefficient in GSUBSURF</param>
	/// <returns></returns>
	bool generalizedGaussSeidelCoefficients(Segment_Image_Data inputImageData, dataType** edgeGradientPtr, Coefficient_Pointers CoefPtrs, 
		Gradient_Pointers VPtrs, Segmentation_Parameters segParameters, dataType coef_dif);

	/// <summary>
	/// The function performs segmentation by generalized subjective surface method (GSUBSURF) for one time step
	/// </summary>
	/// <param name="prevSol_extPtr"> : previous solution</param>
	/// <param name="gauss_seidelPtr"> : current solution</param>
	/// <param name="inputImageData">Structure for the image to be segmented</param>
	/// <param name="segParameters">: Segmentation parameters</param>
	/// <param name="CoefPtrs"> : coefficient of the linear system</param>
	/// <param name="centers"> : centers of the segmentation</param>
	/// <param name="no_of_centers"> : number of the segmentation centers</param>
	/// <returns></returns>
	bool generalizedSubsurfSegmentationTimeStep(dataType** prevSol_extPtr, dataType** gauss_seidelPtr, Segment_Image_Data inputImageData,
		Segmentation_Parameters segParameters, Coefficient_Pointers CoefPtrs, Point3D* centers, size_t no_of_centers);

#ifdef __cplusplus
}
#endif