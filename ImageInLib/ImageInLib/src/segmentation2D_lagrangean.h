#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdbool.h>
#include <stddef.h>
#include "common_functions.h"
#include "segmentation.h"

	// STRUCTs

	/// <summary>
	/// The method segments 2d image by 2d open curve with lagrangean approach (explicit scheme)
	/// </summary>
	/// <param name="inputImage2D">Input mage data</param>
	/// <param name="pSegmentationParams">Parameters needed for segmentation</param>
	/// <param name="pOutputPathPtr">Destination path, where resulting segmentation should be strored</param>
	/// <param name="resultSegmentation">Resulting Curve2D segmentation curve</param>
	/// <returns>Returns true in case of sucessfully finished segmentatio process. Otherwise returns false</returns>
	bool lagrangeanExplicitOpen2DCurveSegmentation(Image_Data2D inputImage2D, 
		const Lagrangean2DSegmentationParameters* pSegmentationParams, unsigned char* pOutputPathPtr, Curve2D* resultSegmentation);

#ifdef __cplusplus
}
#endif
