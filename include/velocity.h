#ifndef VELOCITY_H
#define VELOCITY_H

#include <utility>

#include "real.h"
#include "config.h"

class Matrix;

std::pair<Matrix,Matrix> compIntermediateVelocity(const Config::geo geoConfig, const Config::constants constantsConfig,
						  const REAL delt, const Matrix& U, const Matrix& V);

/**
 *  @brief Computes the velocity for the next time step using the given intermediate values F,G and the new pressure P
 *  @param geoConfig, delt Grid and time dicretization
 *  @param U,V Contains the new velocity afterwards
 *  @param F,G Intermediate values for the velocity, see @compIntermediateVelocity()
 *  @param P New pressure from @SOR_Poisson()
 */
void compNewVelocity(const Config::geo geoConfig, const REAL delt, Matrix& U, Matrix& V,
		     const Matrix& F, const Matrix& G, const Matrix& P);

#endif
