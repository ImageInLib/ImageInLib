#include "pca.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//==============================================================================
dataType * vector(const size_t n)
{
	dataType *v;

	v = (dataType *)malloc((unsigned)n * sizeof(dataType));
	if (!v) erhand("Allocation failure in vector().");
	return v - 1;
}
//==============================================================================
dataType ** matrix(const size_t m, const size_t n)
{
	int i;
	dataType **mat;

	/* Allocate pointers to rows. */
	mat = (dataType **)malloc((unsigned)(n) * sizeof(dataType*));
	if (!mat) erhand("Allocation failure 1 in matrix().");
	mat -= 1;

	/* Allocate rows and set pointers to them. */
	for (i = 1; i <= n; i++)
	{
		mat[i] = (dataType *)malloc((unsigned)(m) * sizeof(dataType));
		if (!mat[i]) erhand("Allocation failure 2 in matrix().");
		mat[i] -= 1;
	}

	/* Return pointer to array of pointers to rows. */
	return mat;
}
//==============================================================================
/* Householder reduction of matrix a to tridiagonal form.
   Algorithm: Martin et al., Num. Math. 11, 181-195, 1968.
   Ref: Smith et al., Matrix Eigensystem Routines -- EISPACK Guide
		Springer-Verlag, 1976, pp. 489-494.
		W H Press et al., Numerical Recipes in C, Cambridge U P,
		1988, pp. 373-374.
*/
void triDecomp(dataType ** a, const size_t n, dataType * d, dataType * e)
{
	size_t l, k, j, i;
	dataType scale, hh, h, g, f;

	for (i = n; i >= 2; i--)
	{
		l = i - 1;
		h = scale = 0.0;
		if (l > 1)
		{
			for (k = 1; k <= l; k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)
				e[i] = a[i][l];
			else
			{
				for (k = 1; k <= l; k++)
				{
					a[i][k] /= scale;
					h += a[i][k] * a[i][k];
				}
				f = a[i][l];
				g = f > 0 ? -sqrt(h) : sqrt(h);
				e[i] = scale * g;
				h -= f * g;
				a[i][l] = f - g;
				f = 0.0;
				for (j = 1; j <= l; j++)
				{
					a[j][i] = a[i][j] / h;
					g = 0.0;
					for (k = 1; k <= j; k++)
						g += a[j][k] * a[i][k];
					for (k = j + 1; k <= l; k++)
						g += a[k][j] * a[i][k];
					e[j] = g / h;
					f += e[j] * a[i][j];
				}
				hh = f / (h + h);
				for (j = 1; j <= l; j++)
				{
					f = a[i][j];
					e[j] = g = e[j] - hh * f;
					for (k = 1; k <= j; k++)
						a[j][k] -= (f * e[k] + g * a[i][k]);
				}
			}
		}
		else
			e[i] = a[i][l];
		d[i] = h;
	}
	d[1] = 0.0;
	e[1] = 0.0;
	for (i = 1; i <= n; i++)
	{
		l = i - 1;
		if (d[i])
		{
			for (j = 1; j <= l; j++)
			{
				g = 0.0;
				for (k = 1; k <= l; k++)
					g += a[i][k] * a[k][j];
				for (k = 1; k <= l; k++)
					a[k][j] -= g * a[k][i];
			}
		}
		d[i] = a[i][i];
		a[i][i] = 1.0;
		for (j = 1; j <= l; j++)
			a[j][i] = a[i][j] = 0.0;
	}
}
//==============================================================================
// dataTyperidiagonal QL (Q- Orthogonal, R -Upper triangular) implicit Algorithm
void triDian(dataType d[], dataType e[], const size_t n, dataType ** z)
{
	int m, l, iter, i, k;
	dataType s, r, p, g, f, dd, c, b;
	void erhand();

	for (i = 2; i <= n; i++)
		e[i - 1] = e[i];
	e[n] = 0.0;
	for (l = 1; l <= n; l++)
	{
		iter = 0;
		do
		{
			for (m = l; m <= n - 1; m++)
			{
				dd = fabs(d[m]) + fabs(d[m + 1]);
				if (fabs(e[m]) + dd == dd) break;
			}
			if (m != l)
			{
				if (iter++ == 30) erhand("No convergence in TLQI.");
				g = (d[l + 1] - d[l]) / (2.0 * e[l]);
				r = sqrt((g * g) + 1.0);
				g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
				s = c = 1.0;
				p = 0.0;
				for (i = m - 1; i >= l; i--)
				{
					f = s * e[i];
					b = c * e[i];
					if (fabs(f) >= fabs(g))
					{
						c = g / f;
						r = sqrt((c * c) + 1.0);
						e[i + 1] = f * r;
						c *= (s = 1.0 / r);
					}
					else
					{
						s = f / g;
						r = sqrt((s * s) + 1.0);
						e[i + 1] = g * r;
						s *= (c = 1.0 / r);
					}
					g = d[i + 1] - p;
					r = (d[i] - g) * s + 2.0 * c * b;
					p = s * r;
					d[i + 1] = g + p;
					g = c * r - b;
					for (k = 1; k <= n; k++)
					{
						f = z[k][i + 1];
						z[k][i + 1] = s * z[k][i] + c * f;
						z[k][i] = c * z[k][i] - s * f;
					}
				}
				d[l] = d[l] - p;
				e[l] = g;
				e[m] = 0.0;
			}
		} while (m != l);
	}
}
//==============================================================================
void erhand(char err_msg[])
{
	fprintf(stderr, "Run-time error:\n");
	fprintf(stderr, "%s\n", err_msg);
	fprintf(stderr, "Exiting to system.\n");
	exit(1);
}
//==============================================================================
void free_vector(dataType * v, const size_t n)
{
	free((char*)(v + 1));
}
//==============================================================================
void free_matrix(dataType ** mat, const size_t n, const size_t m)
{
	size_t i;

	for (i = n; i >= 1; i--)
	{
		free((char*)(mat[i] + 1));
	}
	free((char*)(mat + 1));
}
void eigenSort(dataType * eigval, dataType ** eigvectors, const size_t dim)
{
	size_t k, i, j;
	dataType p;
	for (i = 1; i < dim; i++)
	{
		k = i;
		p = eigval[i];
		for (j = i + 1; j <= dim; j++)
		{
			if (eigval[j] >= p)
			{
				k = j;
				p = eigval[j];
			}
		}
		if (k != i)
		{
			eigval[k] = eigval[i];
			eigval[i] = p;
			for (j = 1; j <= dim; j++)
			{
				p = eigvectors[j][i];
				eigvectors[j][i] = eigvectors[j][k];
				eigvectors[j][k] = p;
			}
		}
	}
}
void covcol(dataType ** dta, const size_t n, const size_t m, dataType ** symmat)
{
	dataType *mean, *vector();
	int i, j, j1, j2;

	/* Allocate storage for mean vector */
	mean = vector(m);
	/* Determine mean of column vectors of input data matrix */
	for (j = 1; j <= m; j++)
	{
		mean[j] = 0.0;
		for (i = 1; i <= n; i++)
		{
			mean[j] += dta[i][j];
		}
		mean[j] /= (dataType)n;
	}
	/* Center the column vectors. */
	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= m; j++)
		{
			dta[i][j] -= mean[j];
		}
	}
	/* Calculate the m * m covariance matrix. */
	for (j1 = 1; j1 <= m; j1++)
	{
		for (j2 = j1; j2 <= m; j2++)
		{
			symmat[j1][j2] = 0.0;
			for (i = 1; i <= n; i++)
			{
				symmat[j1][j2] += dta[i][j1] * dta[i][j2];
			}
			symmat[j2][j1] = symmat[j1][j2];
		}
	}
	// Free
	free_vector(mean, m);
}
//==============================================================================