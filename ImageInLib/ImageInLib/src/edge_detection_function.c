/*
* Author: Markjoe Olunna UBA
* Purpose: ImageInLife project - 4D Image Segmentation Methods
* Language:  C
*/
#include <stdio.h>
#include <stdlib.h>
#include "edgedetection.h"
#include "data_initialization.h"
#include "thresholding.h"
#include "common_functions.h"

bool edgeDetection3dFunctionD(dataType ** image3DPtr, dataType ** edge3DPtr, const size_t xDim,
	const size_t yDim, const size_t zDim, const dataType bgroundvalue, const dataType fgroundvalue)
{
	size_t i, k;//loop counter for z dimension
	const size_t dim2D = xDim * yDim;

	//checks if the memory was allocated
	if (image3DPtr == NULL)
		return false;
	if (edge3DPtr == NULL)
		return false;

	//thresholding of voxel values
	thresholding3dFunctionD(image3DPtr, xDim, yDim, zDim, (bgroundvalue + fgroundvalue) / 2, bgroundvalue, fgroundvalue);

	// initialize the edge array to background value
	initialize3dArrayD(edge3DPtr, xDim, yDim, zDim, bgroundvalue);

	// filling the 2d array with number=value
	for (k = 0; k < zDim; k++)
	{
		for (i = 0; i < dim2D; i++) {
			if (image3DPtr[k][i] == fgroundvalue)//while on the object
			{
				if (k == 0)
				{
					if (i == 0)//if at the top left corner,
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue)
							|| (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == (xDim - 1)) //if at the top right coner,
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue)
							|| (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == ((yDim - 1) * xDim)) //if at the bottom left corner,
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue)
							|| (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == (dim2D - 1)) //if at the bottom right corner,
					{// check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue)
							|| (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i > 0) && (i < (xDim - 1)))//if at first row (excluding left and right corner),
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i + xDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i >((yDim - 1) * xDim)) && (i < (dim2D - 1))) //if at last row
					{//(excluding left and right corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i != 0) && (i != ((yDim - 1) * xDim)) && ((i % xDim) == 0))//if at first column
					{//(excluding top and bottom left corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) ||
							(image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i != (xDim - 1)) && (i != (dim2D - 1)) && ((i % xDim) == (xDim - 1)))//if at last column
					{//(excluding top and bottom right corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) ||
							(image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else {//otherwise, check if pixel neighbors has background value

						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) ||
							(image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
				}

				else if (k == (zDim - 1))
				{
					if (i == 0)//if at the top left corner,
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue)
							|| (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == (xDim - 1)) //if at the top right coner,
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue)
							|| (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == ((yDim - 1) * xDim)) //if at the bottom left corner,
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue)
							|| (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == (dim2D - 1)) //if at the bottom right corner,
					{// check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue)
							|| (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i > 0) && (i < (xDim - 1)))//if at first row (excluding left and right corner),
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i + xDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i >((yDim - 1) * xDim)) && (i < (dim2D - 1))) //if at last row
					{//(excluding left and right corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i != 0) && (i != ((yDim - 1) * xDim)) && ((i % xDim) == 0))//if at first column
					{//(excluding top and bottom left corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) ||
							(image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i != (xDim - 1)) && (i != (dim2D - 1)) && ((i % xDim) == (xDim - 1)))//if at last column
					{//(excluding top and bottom right corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) ||
							(image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else {//otherwise, check if pixel neighbors has background value

						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
				}

				else
				{
					if (i == 0)//if at the top left corner,
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == (xDim - 1)) //if at the top right coner,
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == ((yDim - 1) * xDim)) //if at the bottom left corner,
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == (dim2D - 1)) //if at the bottom right corner,
					{// check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i > 0) && (i < (xDim - 1)))//if at first row (excluding left and right corner),
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i + xDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i >((yDim - 1) * xDim)) && (i < (dim2D - 1))) //if at last row
					{//(excluding left and right corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i != 0) && (i != ((yDim - 1) * xDim)) && ((i % xDim) == 0))//if at first column
					{//(excluding top and bottom left corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) ||
							(image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i != (xDim - 1)) && (i != (dim2D - 1)) && ((i % xDim) == (xDim - 1)))//if at last column
					{//(excluding top and bottom right corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) ||
							(image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else {//otherwise, check if pixel neighbors has background value

						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
				}
			}

		}

	}

	return true;
}


bool edgeDetection3dFunctionUC(unsigned char ** image3DPtr, unsigned char ** edge3DPtr, const size_t xDim, 
	const size_t yDim, 	const size_t zDim, const unsigned char bgroundvalue, const unsigned char fgroundvalue)
{
	size_t i, k;//loop counter for z dimension
	const size_t dim2D = xDim * yDim;

	//checks if the memory was allocated
	if (image3DPtr == NULL)
		return false;
	if (edge3DPtr == NULL)
		return false;

	// initialize the edge array to background value
	initialize3dArrayUC(edge3DPtr, xDim, yDim, zDim, bgroundvalue);

	// filling the 2d array with number=value
	for (k = 0; k < zDim; k++)
	{
		for (i = 0; i < dim2D; i++) {
			if (image3DPtr[k][i] == fgroundvalue)//while on the object
			{ 
				if (k == 0)
				{
					if (i == 0)//if at the top left corner, 
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue)
							|| (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == (yDim - 1)) //if at the top right coner, 
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue)
							|| (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == ((xDim - 1) * yDim)) //if at the bottom left corner, 
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i - yDim] == bgroundvalue)
							|| (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == (dim2D - 1)) //if at the bottom right corner,
					{// check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i - yDim] == bgroundvalue)
							|| (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i > 0) && (i < (yDim - 1)))//if at first row (excluding left and right corner),
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i + yDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i >((xDim - 1) * yDim)) && (i < (dim2D - 1))) //if at last row 
					{//(excluding left and right corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i != 0) && (i != ((xDim - 1) * yDim)) && ((i % yDim) == 0))//if at first column 
					{//(excluding top and bottom left corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i != (yDim - 1)) && (i != (dim2D - 1)) && ((i % yDim) == (yDim - 1)))//if at last column 
					{//(excluding top and bottom right corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else {//otherwise, check if pixel neighbors has background value

						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
				}
				
				else if (k == (zDim - 1))
				{
					if (i == 0)//if at the top left corner, 
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue)
							|| (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == (yDim - 1)) //if at the top right coner, 
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue)
							|| (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == ((xDim - 1) * yDim)) //if at the bottom left corner, 
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i - yDim] == bgroundvalue)
							|| (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == (dim2D - 1)) //if at the bottom right corner,
					{// check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i - yDim] == bgroundvalue)
							|| (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i > 0) && (i < (yDim - 1)))//if at first row (excluding left and right corner),
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i + yDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i >((xDim - 1) * yDim)) && (i < (dim2D - 1))) //if at last row 
					{//(excluding left and right corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i != 0) && (i != ((xDim - 1) * yDim)) && ((i % yDim) == 0))//if at first column 
					{//(excluding top and bottom left corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i != (yDim - 1)) && (i != (dim2D - 1)) && ((i % yDim) == (yDim - 1)))//if at last column 
					{//(excluding top and bottom right corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else {//otherwise, check if pixel neighbors has background value

						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
				}

				else
				{
					if (i == 0)//if at the top left corner, 
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == (yDim - 1)) //if at the top right coner, 
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == ((xDim - 1) * yDim)) //if at the bottom left corner, 
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i - yDim] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == (dim2D - 1)) //if at the bottom right corner,
					{// check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i - yDim] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i > 0) && (i < (yDim - 1)))//if at first row (excluding left and right corner),
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i + yDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i >((xDim - 1) * yDim)) && (i < (dim2D - 1))) //if at last row 
					{//(excluding left and right corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i != 0) && (i != ((xDim - 1) * yDim)) && ((i % yDim) == 0))//if at first column 
					{//(excluding top and bottom left corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i != (yDim - 1)) && (i != (dim2D - 1)) && ((i % yDim) == (yDim - 1)))//if at last column 
					{//(excluding top and bottom right corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else {//otherwise, check if pixel neighbors has background value

						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
				}
			}

		}

	}

	return true;
}

//function for edgedetection in 2D.
bool edgeDetection2dFunctionUC(unsigned char * image2DPtr, unsigned char * edge2DPtr, const size_t xDim, 
	const size_t yDim, 	const unsigned char bgroundvalue, const unsigned char fgroundvalue)
{
	const size_t dim2D = xDim * yDim;
	size_t i;

	if (image2DPtr == NULL)
		return false;
	if (edge2DPtr == NULL)
		return false;

	// Edge detection in 2D
	for (i = 0; i < dim2D; i++) {
		if (image2DPtr[i] == fgroundvalue)//while on the object
		{
			if (i == 0)//if at the top left corner, 
			{//check if pixel neighbors has background value
				if ((image2DPtr[i + 1] == bgroundvalue) || (image2DPtr[i + yDim] == bgroundvalue))
					edge2DPtr[i] = fgroundvalue;
				else
					edge2DPtr[i] = bgroundvalue;
			}
			else if (i == (yDim - 1)) //if at the top right coner, 
			{//check if pixel neighbors has background value
				if ((image2DPtr[i - 1] == bgroundvalue) || (image2DPtr[i + yDim] == bgroundvalue))
					edge2DPtr[i] = fgroundvalue;
				else
					edge2DPtr[i] = bgroundvalue;
			}
			else if (i == ((xDim - 1) * yDim)) //if at the bottom left corner, 
			{//check if pixel neighbors has background value
				if ((image2DPtr[i + 1] == bgroundvalue) || (image2DPtr[i - yDim] == bgroundvalue))
					edge2DPtr[i] = fgroundvalue;
				else
					edge2DPtr[i] = bgroundvalue;
			}
			else if (i == (dim2D - 1)) //if at the bottom right corner,
			{// check if pixel neighbors has background value
				if ((image2DPtr[i - 1] == bgroundvalue) || (image2DPtr[i - yDim] == bgroundvalue))
					edge2DPtr[i] = fgroundvalue;
				else
					edge2DPtr[i] = bgroundvalue;
			}
			else if ((i > 0) && (i < (yDim - 1)))//if at first row (excluding left and right corner),
			{//check if pixel neighbors has background value
				if ((image2DPtr[i - 1] == bgroundvalue) || (image2DPtr[i + 1] == bgroundvalue) ||
					(image2DPtr[i + yDim] == bgroundvalue))
					edge2DPtr[i] = fgroundvalue;
				else
					edge2DPtr[i] = bgroundvalue;
			}
			else if ((i >((xDim - 1) * yDim)) && (i < (dim2D - 1))) //if at last row 
			{//(excluding left and right corner), check if pixel neighbors has background value
				if ((image2DPtr[i - 1] == bgroundvalue) || (image2DPtr[i + 1] == bgroundvalue) ||
					(image2DPtr[i - yDim] == bgroundvalue))
					edge2DPtr[i] = fgroundvalue;
				else
					edge2DPtr[i] = bgroundvalue;
			}
			else if ((i != 0) && (i != ((xDim - 1) * yDim)) && ((i % yDim) == 0))//if at first column 
			{//(excluding top and bottom left corner), check if pixel neighbors has background value
				if ((image2DPtr[i + 1] == bgroundvalue) || (image2DPtr[i + yDim] == bgroundvalue) ||
					(image2DPtr[i - yDim] == bgroundvalue))
					edge2DPtr[i] = fgroundvalue;
				else
					edge2DPtr[i] = bgroundvalue;
			}
			else if ((i != (yDim - 1)) && (i != (dim2D - 1)) && ((i % yDim) == (yDim - 1)))//if at last column 
			{//(excluding top and bottom right corner), check if pixel neighbors has background value
				if ((image2DPtr[i - 1] == bgroundvalue) || (image2DPtr[i + yDim] == bgroundvalue) ||
					(image2DPtr[i - yDim] == bgroundvalue))
					edge2DPtr[i] = fgroundvalue;
				else
					edge2DPtr[i] = bgroundvalue;
			}
			else {//otherwise, check if pixel neighbors has background value

				if ((image2DPtr[i - 1] == bgroundvalue) || (image2DPtr[i + 1] == bgroundvalue) ||
					(image2DPtr[i - yDim] == bgroundvalue) || (image2DPtr[i + yDim] == bgroundvalue))
					edge2DPtr[i] = fgroundvalue;
				else
					edge2DPtr[i] = bgroundvalue;
			}

		}

	}


	return true;
}