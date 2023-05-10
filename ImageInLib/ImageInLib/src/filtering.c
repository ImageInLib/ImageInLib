#include "filtering.h"
#include "heat_equation.h"
#include "non_linear_heat_equation.h"
#include "filter_params.h"
#include "common_functions.h"

void filterImage(Image_Data inputImageData, Filter_Parameters filterParameters, const FilterMethod method)
{
	switch (method)
	{
	case LINEAR_HEATEQUATION_EXPLICIT:
		heatExplicitScheme(inputImageData, filterParameters);
		break;
	case LINEAR_HEATEQUATION_IMPLICIT:
		heatImplicitScheme(inputImageData, filterParameters);
		break;
	case NONLINEAR_HEATEQUATION_EXPLICIT:
		nonLinearHeatExplicitScheme(inputImageData, filterParameters);
		break;
	case NONLINEAR_HEATEQUATION_IMPLICIT:
		nonLinearHeatImplicitScheme(inputImageData, filterParameters);
		break;
	case MEAN_CURVATURE_FILTER:
		meanCurvatureTimeStep(inputImageData, filterParameters);
		break;
	case GEODESIC_MEAN_CURVATURE_FILTER:
		geodesicMeanCurvatureTimeStep(inputImageData, filterParameters);
		break;
	default:
		break;
	}
}