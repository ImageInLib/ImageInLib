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
#include "generate_3D_shapes.h"


bool ballsOnCircle(dataType ** image3DPtr, dataType radius, size_t xDim, size_t yDim,
	size_t zDim, dataType xCordctr, dataType yCordctr, dataType bRadius, dataType fillValue)
{
	if (image3DPtr == NULL)
		return false;

	dataType pi, t;
	pi = 2 * M_PI;
	Point3D ballCenter;

	for (t = 0; t <= pi; t += 0.06)
	{
		ballCenter.x = (dataType)(xCordctr + radius * cos(t));
		ballCenter.y = (dataType)(yCordctr + radius * sin(t));
		ballCenter.z = (dataType)(zDim / 2);
		fillBall3D(image3DPtr, zDim, xDim, yDim, bRadius, ballCenter, fillValue);
	}
	return true;
}