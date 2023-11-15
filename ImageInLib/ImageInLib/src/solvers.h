#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef SOLVERS
#define SOLVERS

#include "common_functions.h"

    /// <summary>
    /// Sherman-Morris formula for solving of tridiagonal cyclic systems of linear equations
    /// </summary>
    /// <param name="a">upper diagonal</param>
    /// <param name="b">central diagonal</param>
    /// <param name="c">lower diagonal</param>
    /// <param name="alpha">upeer right corner value</param>
    /// <param name="beta">bottom left corner value</param>
    /// <param name="x">vector of unknown values</param>
    /// <param name="ps">righ hand side values</param>
    /// <param name="N">matrix dimension (it is square matrix)</param>
     bool sherman_morris(dataType* a, dataType* b, dataType* c, dataType alpha, dataType beta, dataType* x, dataType* ps, const size_t N);
 
    /// <summary>
    /// Thomas algorithm for solving of tridiagonal system of linear equations
    /// </summary>
    /// <param name="a">upper diagonal</param>
    /// <param name="b">central diagonal</param>
    /// <param name="c">lower diagonal</param>
    /// <param name="x">vector of unknown values</param>
    /// <param name="ps">right hand side values</param>
    /// <param name="N">matrix dimension (it is square matrix)</param>
    bool thomas(dataType* a, dataType* b, dataType* c, dataType* x, dataType* ps, const size_t N);
 

#endif // !SOLVERS

#ifdef __cplusplus
}
#endif