/*
* Author: Markjoe Olunna UBA
* Purpose: ImageInLife project - 4D Image Segmentation Methods
* Language:  C
*/
#include <stdio.h>
#include <stdlib.h>
#include "common_math.h"
#include "trajectories.h"
#include "common_functions.h"
#include "generate_3d_shapes.h"

bool ballsOnHelix(dataType ** image3DPtr, dataType radius, size_t xDim, size_t yDim,
	size_t zDim, size_t pitch, dataType bRadius, dataType fillValue)
{
	if (image3DPtr == NULL)
		return false;

	double t, pi;
	pi = 2 * M_PI * pitch;// * pitch
	Point3D ballCenter;

	for (t = 0; t <= pi; t += 0.6)
	{
		ballCenter.x = (dataType)((xDim/2) + radius * cos(t));
		ballCenter.y = (dataType)((yDim/2) + radius * sin(t));
		ballCenter.z = (dataType)((zDim / pi) * t);
		fillBall3D(image3DPtr, zDim, xDim, yDim, bRadius, ballCenter, fillValue);
	}
	return true;
}
