#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef SOLVERS
#define SOLVERS

#include "common_functions.h"

    typedef struct SchemeData
    {
        double curvature;
        double alfa;
        double beta;
        double beta_ps_expl;
        double ps;
        double a;
        double b;
        double c;
        double sol;
        double f;
        double m;

        //helper attributes - solvers
        //sherman
        double bb;

        //thomas
        double thomas_a;
        double thomas_ps;
        double thomas_b;
        double thomas_c;
        double thomas_x;
    } SchemeData;


    ///// <summary>
    ///// Sherman-Morris formula for solving of tridiagonal cyclic systems of linear equations
    ///// </summary>
    ///// <param name="a">upper diagonal</param>
    ///// <param name="b">central diagonal</param>
    ///// <param name="c">lower diagonal</param>
    ///// <param name="alpha">upeer right corner value</param>
    ///// <param name="beta">bottom left corner value</param>
    ///// <param name="x">vector of unknown values</param>
    ///// <param name="ps">righ hand side values</param>
    ///// <param name="N">matrix dimension (it is square matrix)</param>
    // bool sherman_morris(dataType* a, dataType* b, dataType* c, dataType alpha, dataType beta, dataType* x, dataType* ps, const size_t N);
 
    ///// <summary>
    ///// Thomas algorithm for solving of tridiagonal system of linear equations
    ///// </summary>
    ///// <param name="a">upper diagonal</param>
    ///// <param name="b">central diagonal</param>
    ///// <param name="c">lower diagonal</param>
    ///// <param name="x">vector of unknown values</param>
    ///// <param name="ps">right hand side values</param>
    ///// <param name="N">matrix dimension (it is square matrix)</param>
    //bool thomas(dataType* a, dataType* b, dataType* c, dataType* x, dataType* ps, const size_t N);
 

    /// <summary>
    /// Solves tri-diagonal cyclic system of egautions given by coefficients in pscheme_data
    /// </summary>
    /// <param name="pscheme_data">pointer to coefficients  of 3 diagonal system of equations</param>
    /// <param name="number_of_points">dimension of system of equations</param>
    /// <returns>true if the solver succeeds, otherwise false</returns>
    bool sherman_morris(SchemeData* pscheme_data, const size_t number_of_points);
    /// <summary>
    /// Solves tridiagonal system of equations by thomas algorithm
    /// </summary>
    /// <param name="pscheme_data">pointer to coefficients  of 3 diagonal system of equations</param>
    /// <param name="number_of_points">dimension of system of equations</param>
    /// <returns>true if the solver succeeds, otherwise false</returns>
    bool thomas(SchemeData* pscheme_data, const size_t number_of_points);


#endif // !SOLVERS

#ifdef __cplusplus
}
#endif