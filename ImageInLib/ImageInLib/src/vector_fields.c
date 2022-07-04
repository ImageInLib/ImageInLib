#include "vector_fields.h"
#include <math.h>
//==============================================================================
/*
* Function generating 3D vector field
*/
void generate3DVector(Point3D ** vectorPtr, Image_Data inputPtr, Vector_Parameters vectorVariables, enum FiniteDifference fieldDirection, dataType coeff)
{
	size_t k, x;
	const size_t height = inputPtr.height, _width = inputPtr.length, p = vectorVariables.p;
	dataType h = vectorVariables.h;

	// Mid center Pointers Pointers for the 3 Dimensions
	Vector_Direction zCent = { NULL }, xCent = { NULL }, yCent = { NULL };
	// Calculate the Mid center before
	for (k = p; k <= height + p; k++)
	{
		for (x = p; x <= _width + p; x++)
		{
			// 2D to 3D representation for i, j
			// Assuming dim i = dim j ...
			int i = (int)floor((double)x / _width + 0.5);
			int j = (int)floor((double)(x - i) / _width + 0.5);
			// Fill the central pointers
			xCent.fieldPtr[k][x] = (inputPtr.imageDataPtr[k][x] + inputPtr.imageDataPtr[k][x_new(i - 1, j - 1, _width)] + inputPtr.imageDataPtr[k][x_new(i, j - 1, _width)] + inputPtr.imageDataPtr[k][x_new(i - 1, j, _width)]) / 4;
			yCent.fieldPtr[k][x] = (inputPtr.imageDataPtr[k][x] + inputPtr.imageDataPtr[k - 1][x_new(i - 1, j, _width)] + inputPtr.imageDataPtr[k - 1][x] + inputPtr.imageDataPtr[k][x_new(i - 1, j, _width)]) / 4;
			zCent.fieldPtr[k][x] = (inputPtr.imageDataPtr[k][x] + inputPtr.imageDataPtr[k - 1][x_new(i, j - 1, _width)] + inputPtr.imageDataPtr[k - 1][x] + inputPtr.imageDataPtr[k][x_new(i, j - 1, _width)]) / 4;
		}
	}
	// Coordinate Direction - x Direction, y Direction, z Direction
	dataType zD, xD, yD;

	for (k = p; k <= height + p; k++)
	{
		for (x = p; x <= _width + p; x++)
		{
			// 2D to 3D representation for i, j
			// Assuming dim i = dim j ...
			int i = (int)floor((double)x / _width + 0.5);
			int j = (int)floor((double)(x - i) / _width + 0.5);
			// Checks the direction passed and fill values for that
			if (fieldDirection == FINITE_FORWARD)
			{
				// East
				xD = (inputPtr.imageDataPtr[k][x_new(i, j + 1, _width)] - inputPtr.imageDataPtr[k][x]) / h;
				yD = (xCent.fieldPtr[k][x_new(i + 1, j + 1, _width)] - xCent.fieldPtr[k][x_new(i, j + 1, _width)]) / h;// Y Direction
				zD = (zCent.fieldPtr[k + 1][x_new(i, j + 1, _width)] - zCent.fieldPtr[k][x_new(i, j + 1, _width)]) / h;// Z Direction
				vectorPtr[k][x].x = gradientFunction(xD*xD + yD * yD + zD * zD, coeff);
				// South
				xD = (inputPtr.imageDataPtr[k][x_new(i + 1, j, _width)] - inputPtr.imageDataPtr[k][x]) / h; // Y Direction
				yD = (xCent.fieldPtr[k][x_new(i + 1, j + 1, _width)] - xCent.fieldPtr[k][x_new(i + 1, j, _width)]) / h;// X Direction
				zD = (yCent.fieldPtr[k + 1][x_new(i + 1, j, _width)] - yCent.fieldPtr[k][x_new(i + 1, j, _width)]) / h;// Z Direction
				vectorPtr[k][x].y = gradientFunction(xD*xD + yD * yD + zD * zD, coeff);
				// Top
				xD = (inputPtr.imageDataPtr[k + 1][x] - inputPtr.imageDataPtr[k][x]) / h; // Z Direction
				yD = (yCent.fieldPtr[k + 1][x_new(i + 1, j, _width)] - yCent.fieldPtr[k + 1][x]) / h;// Y Direction
				zD = (zCent.fieldPtr[k + 1][x_new(i, j + 1, _width)] - zCent.fieldPtr[k + 1][x]) / h;// X Direction
				vectorPtr[k][x].z = gradientFunction(xD*xD + yD * yD + zD * zD, coeff);
			}
			else if (fieldDirection == FINITE_BACKWARD)
			{
				// West
				xD = (inputPtr.imageDataPtr[k][x] - inputPtr.imageDataPtr[k][x_new(i, j - 1, _width)]) / h;
				yD = (xCent.fieldPtr[k][x_new(i + 1, j + 1, _width)] - xCent.fieldPtr[k][x]) / h;// Y Direction
				zD = (zCent.fieldPtr[k + 1][x] - zCent.fieldPtr[k][x]) / h;// Z Direction
				vectorPtr[k][x].x = gradientFunction(xD*xD + yD * yD + zD * zD, coeff);
				// North
				xD = (inputPtr.imageDataPtr[k][x] - inputPtr.imageDataPtr[k][x_new(i - 1, j, _width)]) / h; // Y Direction
				yD = (xCent.fieldPtr[k][x_new(i, j + 1, _width)] - xCent.fieldPtr[k][x]) / h;// X Direction
				zD = (yCent.fieldPtr[k + 1][x] - yCent.fieldPtr[k][x]) / h;// Z Direction
				vectorPtr[k][x].y = gradientFunction(xD*xD + yD * yD + zD * zD, coeff);
				// Bottom
				xD = (inputPtr.imageDataPtr[k][x] - inputPtr.imageDataPtr[k - 1][x]) / h; // Z Direction
				yD = (yCent.fieldPtr[k][x_new(i + 1, j, _width)] - yCent.fieldPtr[k][x]) / h;// Y Direction
				zD = (zCent.fieldPtr[k][x_new(i, j + 1, _width)] - zCent.fieldPtr[k][x]) / h;// X Direction
				vectorPtr[k][x].z = gradientFunction(xD*xD + yD * yD + zD * zD, coeff);
			}
			else if (fieldDirection == FINITE_CENTRAL)
			{
				vectorPtr[k][x].x = 0;
				vectorPtr[k][x].y = 0;
				vectorPtr[k][x].z = 0;
			}
		}
	}
}
//==============================================================================