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

    /// <summary>
    /// Solves tridiagonal system of equations by thomas algorithm, without shermann-morisson updates
    /// </summary>
    /// <param name="pscheme_data">pointer to coefficients  of 3 diagonal system of equations</param>
    /// <param name="number_of_points">dimension of system of equations</param>
    /// <returns>true if the solver succeeds, otherwise false</returns>
    bool calculate_by_thomas(SchemeData* pscheme_data, const size_t number_of_points);

    typedef struct SchemeData3D
    {
        double k1;
        double k2;
        double u;
        double v;

        double alfa;
        double ps;

        double a;
        double b;
        double c;
        double sol;
        double m;

        ////helper attributes - solvers
        ////sherman
        //double bb;

        //thomas
        double thomas_a;
        double thomas_b;
        double thomas_c;
        double thomas_ps;
        double thomas_x;

    } SchemeData3D;

    /// <summary>
    /// Solves tridiagonal system of equations by thomas algorithm
    /// </summary>
    /// <param name="pscheme_data">pointer to coefficients  of system of equations</param>
    /// <param name="number_of_points">dimension of system of equations</param>
    /// <returns>true if the solver succeeds, otherwise false</returns>
    bool thomas3D(SchemeData3D* pscheme_data, const size_t number_of_points);

    /// <summary>
    /// Solves tridiagonal system of equations by thomas algorithm, without shermann-morisson updates
    /// </summary>
    /// <param name="pscheme_data">pointer to coefficients  of system of equations</param>
    /// <param name="number_of_points">dimension of system of equations</param>
    /// <returns>true if the solver succeeds, otherwise false</returns>
    bool calculate_by_thomas3D(SchemeData3D* pscheme_data, const size_t number_of_points);

    /// <summary>
    /// Solves tri-diagonal cyclic system of equations given by coefficients in pscheme_data
    /// </summary>
    /// <param name="pscheme_data">pointer to coefficients  of system of equations</param>
    /// <param name="number_of_points">dimension of system of equations</param>
    /// <returns>true if the solver succeeds, otherwise false</returns>
    bool sherman_morris3D(SchemeData3D* pscheme_data, const size_t number_of_points);

#endif // !SOLVERS

#ifdef __cplusplus
}
#endif