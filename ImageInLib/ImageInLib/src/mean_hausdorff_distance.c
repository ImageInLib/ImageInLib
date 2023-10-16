//==============================================================================
#include "mean_hausdorff_distance.h"
//==============================================================================
#include <float.h>
#include <math.h>
//==============================================================================
// MHD Function
dataType mean_hausdorff(dataType ** curveA_Pointer, dataType ** curveB_Pointer, dataType curveIntensity, size_t height, size_t length, size_t width)
{
	dataType mhd = 0., dist, dz, dx, dy;
	double shortest;
	//==============================================================================
	size_t k, i, j, xd;
	size_t pts_a = 0, pts_b = 0;
	const size_t dim2D = (length)*(width);
	const size_t dim3D = (height)*dim2D;
	Point3D * surface_points_a = (Point3D *)malloc(sizeof(Point3D) * dim3D);
	Point3D * surface_points_b = (Point3D *)malloc(sizeof(Point3D) * dim3D);
	//==============================================================================
	// Save the zero surface points
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				// 1D Conversion of row and column
				xd = x_new(i, j, length);
				if (curveA_Pointer[k][xd] == curveIntensity)
				{
					surface_points_a[pts_a].x = (dataType)i;
					surface_points_a[pts_a].y = (dataType)j;
					surface_points_a[pts_a].z = (dataType)k;
					pts_a++;
				}
				if (curveB_Pointer[k][xd] == curveIntensity)
				{
					surface_points_b[pts_b].x = (dataType)i;
					surface_points_b[pts_b].y = (dataType)j;
					surface_points_b[pts_b].z = (dataType)k;
					pts_b++;
				}
			}
		}
	}
	//==============================================================================
	double *dist1 = (double *)malloc(pts_a * sizeof(double));
	//==============================================================================
	for (i = 0; i < pts_a; i++)
	{
		shortest = DBL_MAX;
		dist1[i] = 0;
		for (j = 0; j < pts_b; j++)
		{
			dz = surface_points_a[i].z - surface_points_b[j].z;
			dx = surface_points_a[i].x - surface_points_b[j].x;
			dy = surface_points_a[i].y - surface_points_b[j].y;
			dist = dx * dx + dy * dy + dz * dz;
			dist = (dataType)sqrt(dist);
			if (dist < shortest)
			{
				shortest = dist;
			}
		}
		dist1[i] = (shortest) / (double)pts_a;
	}
	//==============================================================================
	for (i = 0; i < pts_a; i++)
	{
		mhd += (dataType)dist1[i];
	}
	//==============================================================================
	// Free Pointers
	free(dist1);
	free(surface_points_a);
	free(surface_points_b);
	dist1 = NULL;
	surface_points_a = NULL;
	surface_points_b = NULL;
	//==============================================================================
	return mhd;
	//==============================================================================
}
