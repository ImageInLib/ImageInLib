#include "segmentation3D_common.h"
#include <math.h>

dataType chooseGamma(dataType value, dataType dist1, dataType dist2)
{
	if (fabs(value) > dist1)
	{
		return 0.;
	}
	else if (fabs(value) < dist2)
	{
		return 1.;
	}
	else
	{
		return 0.5;
	}
}
