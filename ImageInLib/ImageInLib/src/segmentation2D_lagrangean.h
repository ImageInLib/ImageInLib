#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdbool.h>
#include <stddef.h>
#include "common_functions.h"
#include "segmentation.h"

	// STRUCTs

	/* lagrangean2dOpenCurveSegmentation segments 2d image by 2d open curve with lagrangean approach (explicit scheme)
	*/
	bool lagrangeanExplicitOpen2DCurveSegmentation(Image_Data inputImage2D, const Lagrangean2DSegmentationParameters* pSegmentationParams, 
		unsigned char* pOutputPathPtr, Curve2D* resultSegmentation);

#ifdef __cplusplus
}
#endif
