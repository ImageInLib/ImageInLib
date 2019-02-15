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
dataType *vector(int n); // Allocation of vector storage
dataType **matrix(int m, int n); // Allocation of float matrix storage
void triDecomp(dataType **a, int n, dataType *d, dataType *e); // Reduce a real, symmetric matrix to a symmetric, tridiagonal matrix.
void triDian(dataType d[], dataType e[], int n, dataType **z); // Tridiagonal QL algorithm -- Implicit
void erhand(char err_msg[]);// Error handler
void free_vector(dataType *v, int n); //Deallocate vector storage allocated by vector()
void free_matrix(dataType **mat, int n, int m); // Deallocate float matrix storage
void covcol(dataType ** data, int n, int m, dataType ** symmat);
//==============================================================================
void eigenSort(dataType *eigval, dataType **eigvectors, int dim); // Sorts pointers in descending order of eigenvalues
//==============================================================================
#endif // !PCA
