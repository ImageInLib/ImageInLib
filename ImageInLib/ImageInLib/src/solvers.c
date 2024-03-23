#include <memory.h>
#include <corecrt_malloc.h>

#include "solvers.h"

bool sherman_morris(SchemeData* pscheme_data, const size_t number_of_points)
{
	if (pscheme_data == NULL)
	{
		return false;
	}

	//bb, aa, cc, cc[n], bb[1], u, ps
	//double *a, double * b, double * c, double alpha, double beta, double * x, double * ps, int N

	int i = 0;
	double fact = 0, gamma = 0;
	const double beta = pscheme_data[1].b;
	const double alpha = pscheme_data[number_of_points].c;

	for (size_t i = 0; i < number_of_points + 2; i++)
	{
		pscheme_data[i].bb = pscheme_data[i].a; // toto bolo popletene pomenovanie
	}

	gamma = -pscheme_data[1].a;  //popletene pomenovanie
	pscheme_data[1].bb -= gamma;
	pscheme_data[number_of_points].bb -= alpha * beta / gamma;

	for (size_t i = 0; i < number_of_points + 2; i++)
	{
		pscheme_data[i].thomas_a = pscheme_data[i].b;
		pscheme_data[i].thomas_b = pscheme_data[i].bb;
		pscheme_data[i].thomas_c = pscheme_data[i].c;
		pscheme_data[i].thomas_ps = pscheme_data[i].ps;
	}
	thomas(pscheme_data, number_of_points);

	for (size_t i = 0; i < number_of_points + 2; i++)
	{
		pscheme_data[i].sol = pscheme_data[i].thomas_x;
		pscheme_data[i].thomas_ps = 0;
		pscheme_data[i].thomas_b = pscheme_data[i].bb;
	}

	pscheme_data[1].thomas_ps = gamma;
	pscheme_data[number_of_points].thomas_ps = alpha;

	thomas(pscheme_data, number_of_points);

	fact = (pscheme_data[1].sol + beta * pscheme_data[number_of_points].sol / gamma) /
		(1.0 + pscheme_data[1].thomas_x + beta * pscheme_data[number_of_points].thomas_x / gamma);

	for (size_t i = 1; i < number_of_points + 1; i++)
	{
		pscheme_data[i].sol -= fact * pscheme_data[i].thomas_x;
	}

	pscheme_data[0].sol = pscheme_data[number_of_points].sol;
	pscheme_data[number_of_points + 1].sol = pscheme_data[1].sol;

	return true;
}

bool thomas(SchemeData* pscheme_data, const size_t number_of_points)
{
	if (pscheme_data == NULL)
	{
		return false;
	}

	double m = 0;

	//forward elimination

	for (size_t i = 2; i < number_of_points + 1; i++)
	{
		m = pscheme_data[i].thomas_a / pscheme_data[i - 1].thomas_b;
		pscheme_data[i].thomas_b -= m * pscheme_data[i - 1].c;
		pscheme_data[i].thomas_ps -= m * pscheme_data[i - 1].thomas_ps;
	}

	//backward substitution

	pscheme_data[number_of_points].thomas_x = pscheme_data[number_of_points].thomas_ps / pscheme_data[number_of_points].thomas_b;

	for (size_t i = number_of_points - 1; i >= 1; i--)
	{
		pscheme_data[i].thomas_x = (pscheme_data[i].thomas_ps - pscheme_data[i].thomas_c * pscheme_data[i + 1].thomas_x) / pscheme_data[i].thomas_b;
	}

	return true;
}
