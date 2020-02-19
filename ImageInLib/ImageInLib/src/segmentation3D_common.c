#include "segmentation3D_common.h"
#include <math.h>

dataType approx_U_pq(dataType U_p, dataType U_q, dataType A_pq)
{
	if (A_pq > 0)
	{
		return U_p;
	}
	else if (A_pq < 0)
	{
		return U_q;
	}
}

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
