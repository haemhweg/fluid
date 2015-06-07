#ifndef SOLVERS_H
#define SOLVERS_H

#include <utility>

#include "config.h"
#include "matrix.h"
#include "real.h"
#include "geometry.h"

/**
 *  @brief Solves the discrete 2D poisson equation using the SOR method.
 *  @param P, rhs delta P = rhs
 *  @param geoConfig Information for the discretized grid
 *  @param solverConfig Contains parameter omega, residual error eps and maximal iteration number
 *  @return The iteration number and achieved accuracy of the solution
 */
std::pair<unsigned, REAL> SOR_Poisson(const Config::geo& geoConfig, const Config::solver solverConfig, 
				      const Geometry& geometry, Matrix& P, const Matrix& rhs);

#endif
