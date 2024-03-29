#include "interpolations.h"

/*
* Function to calculate Tri-linear Interpolation
*/
dataType trilinearInterpolation(dataType x, dataType x1, dataType x2, dataType y, dataType y1, dataType y2, dataType c000, dataType c001, dataType c010, dataType c011, dataType c100, dataType c101, dataType c110, dataType c111, dataType z, dataType z1, dataType z2)
{
	dataType c1 = bilinearInterpolation(x, x1, x2, c001, c011, c101, c111, y, y1, y2);
	dataType c2 = bilinearInterpolation(x, x1, x2, c000, c010, c100, c110, y, y1, y2);

	dataType c = linearInterpolation(z, z1, z2, c1, c2); // Interpolated point

	return c;
}
/*
* Function to calculate linear values R1, R2 inputs for Bilinear
*/
dataType linearInterpolation(dataType x, dataType x1, dataType x2, dataType q00, dataType q01)
{
	return ((x2 - x) / (x2 - x1))*q00 + ((x - x1) / (x2 - x1))*q01;// Interpolated point
}
/*
* Function to calculate Bilinear Interpolation
*/
dataType bilinearInterpolation(dataType x, dataType x1, dataType x2, dataType q11, dataType q12, dataType q21, dataType q22, dataType y, dataType y1, dataType y2)
{
	dataType r1 = linearInterpolation(x, x1, x2, q11, q21);
	dataType r2 = linearInterpolation(x, x1, x2, q12, q22);
	dataType p = linearInterpolation(y, y1, y2, r1, r2); // Interpolated point
	return p;
}