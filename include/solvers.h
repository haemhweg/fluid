#ifndef SOLVERS_H
#define SOLVERS_H

#include <utility>

#include "config.h"
#include "matrix.h"
#include "real.h"


std::pair<unsigned, REAL> SOR_Poisson(Matrix& P, const Matrix& rhs, const Config::geo geoConfig, 
				      const Config::solver solverConfig);

#endif
