#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdbool.h>
#include <stddef.h>
#include "common_functions.h"
#include "segmentation.h"
#include "imageInterpolation.h"

	// Structures

	/// <summary>
	/// The method segment 3d image by 3d curve with lagrangean approach (explicit scheme)
	/// </summary>
	/// <param name="inputImage3D">Input mage data</param>
	/// <param name="pSegmentationParams">Parameters needed for segmentation</param>
	/// <param name="pOutputPathPtr">Destination path, where resulting segmentation should be stored</param>
	/// <param name="resultSegmentation">Resulting Curve3D segmentation curve</param>
	/// <returns>Returns true in case of sucessfully finished segmentation process. Otherwise returns false</returns>
	bool lagrangeanExplicit3DCurveSegmentation(Image_Data inputImage3D,
		const Lagrangean3DSegmentationParameters* pSegmentationParams, unsigned char* pOutputPathPtr, Curve3D* resultSegmentation);

	/// <summary>
	/// The function move a 3D curve with lagrangean approach (semi-implicit scheme)
	/// </summary>
	/// <param name="inputImage3D">Input mage data</param>
	/// <param name="pSegmentationParams">Parameters needed for the evolution</param>
	/// <param name="pOutputPathPtr">Destination path, where final curve should be strored</param>
	/// <param name="resultSegmentation">Resulting 3D curve</param>
	/// <returns>Returns true in case of sucessfully finished evolution process. Otherwise returns false</returns>
	bool lagrangeanSemiImplicit3DCurveSegmentation(Image_Data inputImage3D, const Lagrangean3DSegmentationParameters* pSegmentationParams, unsigned char* pOutputPathPtr, Curve3D* pResultSegmentation);

#ifdef __cplusplus
}
#endif
