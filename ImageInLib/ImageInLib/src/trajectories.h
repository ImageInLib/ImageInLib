#ifdef __cplusplus
extern "C" {
#endif
	/*
* Author: Markjoe Olunna UBA
* Purpose: ImageInLife project - 4D Image Segmentation Methods
* Language:  C
*/
#pragma once
#include <stdbool.h>
#include <stddef.h>
#include "common_functions.h"

/* ballsOnHelix is a function that positions balls on a given helix
* pitch of a helix is the height of one complete helix turn, measured parallel to the axis of the helix. That is, if
pitch is set to 5, five complete helix turns will be implemented.
* radius is the radius of the helix turn(0 < radius <= (min{xDim, yDim})/2).
* fillValue is the value the voxel will be set to
* bRadius is the radius of the ball
*xDim, yDim, zDim are the dimensions of x, y and z respectively
*/
	bool ballsOnHelix(dataType ** image3DPtr, dataType radius, size_t xDim, size_t yDim,
		size_t zDim, size_t pitch, dataType bRadius, dataType fillValue);

	/* ballsOnCircle is a function that positions balls on a given circle
	* radius is the radius of the circle(0 < radius <= (min{xDim, yDim})/2).
	* fillValue is the value the voxel will be set to
	* bRadius is the radius of the ball
	*xDim, yDim, zDim are the dimensions of x, y and z respectively
	* xCordctr, yCordctr are the x and y coordinates respectively (of the circle centre)
	*/
	bool ballsOnCircle(dataType ** image3DPtr, dataType radius, size_t xDim, size_t yDim,
		size_t zDim, dataType xCordctr, dataType yCordctr, dataType bRadius, dataType fillValue);

#ifdef __cplusplus
}
#endif