#ifndef SOLVERS_H
#define SOLVERS_H

#include <utility>

#include "config.h"
#include "matrix.h"
#include "real.h"


/**
 *  @brief Computes the right hand side for the discrete 2D poisson equation.
 *  @param geoConfig Information for the discretized grid
 *  @param delt Time step size
 *  @param F, G Intermediate velocity computed by @compIntermediateVelocity()
 *  @return The iteration number and achieved accuracy of the solution
 */
Matrix RHS_Poisson(const Config::geo geoConfig, const REAL delt, const Matrix& F, const Matrix& G);

/**
 *  @brief Solves the discrete 2D poisson equation using the SOR method.
 *  @param P, rhs delta P = rhs
 *  @param geoConfig Information for the discretized grid
 *  @param solverConfig Contains parameter omega, residual error eps and maximal iteration number
 *  @return The iteration number and achieved accuracy of the solution
 */
std::pair<unsigned, REAL> SOR_Poisson(Matrix& P, const Matrix& rhs, const Config::geo geoConfig, 
				      const Config::solver solverConfig);

#endif
