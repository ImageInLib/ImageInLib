#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdbool.h>
#include <stddef.h>
#include "common_functions.h"


	// STRUCTs

//Structure that holds the parameters used during SUBSURF segmentation process.
	typedef struct
	{
		size_t numberOfTimeStep;	// Number of current time step
		Point2D* pInitialCondition;	// Initial curve
		size_t numberOfCurvePoints;// Number of initial segmentation curve points
		bool isCurveClosed;			// The flag defines, if the curve is open (with fixed the 1st and last point) or closed
	} Lagrangean2DSegmentationParameters;

	/* lagrangean2dOpenCurveSegmentation segments 2d image by 2d open curve with lagrangean approach (explicit scheme)
	*/
	bool lagrangeanExplicitOpen2DCurveSegmentation(Image_Data inputImage2D, const Lagrangean2DSegmentationParameters* pSegmentationParams, 
		unsigned char* pOutputPathPtr, Curve2D* resultSegmentation);

#ifdef __cplusplus
}
#endif
