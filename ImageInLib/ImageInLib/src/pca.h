#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef PCA
#define PCA
	/*********************************/
	/* Principal Components Analysis */
	/*********************************/

	/* Principal Components Analysis or the Karhunen-Loeve expansion is a
	   classical method for dimensionality reduction or exploratory data
	   analysis.  One reference among many is: F. Murtagh and A. Heck,
	   Multivariate Data Analysis, Kluwer Academic, Dordrecht, 1987.*/

	   //==============================================================================
#include "common_functions.h"
//==============================================================================
// Macro
#define SIGN(a, b) ( (b) < 0 ? -fabs(a) : fabs(a) ) // Changes the sign
//==============================================================================
	dataType *vector(const size_t n); // Allocation of vector storage
	dataType **matrix(const size_t m, const size_t n); // Allocation of float matrix storage
	void triDecomp(dataType **a, const size_t n, dataType *d, dataType *e); // Reduce a real, symmetric matrix to a symmetric, tridiagonal matrix.
	void triDian(dataType d[], dataType e[], const size_t n, dataType **z); // Tridiagonal QL algorithm -- Implicit
	void erhand(char err_msg[]);// Error handler
	void free_vector(dataType *v, const size_t n); //Deallocate vector storage allocated by vector()
	void free_matrix(dataType **mat, const size_t n, const size_t m); // Deallocate float matrix storage
	void covcol(dataType ** data, const size_t n, const size_t m, dataType ** symmat);
	//==============================================================================
	void eigenSort(dataType *eigval, dataType **eigvectors, const size_t dim); // Sorts pointers in descending order of eigenvalues
	//==============================================================================
#endif // !PCA

#ifdef __cplusplus
}
#endif