#include "real.h"
#include "config.h"
#include "matrix.h"

REAL dU2dx(Config::geo geoConfig, Config::constants constantsConfig, Matrix U, size_t i, size_t j);
REAL dUVdy(Config::geo geoConfig, Config::constants constantsConfig, Matrix U, Matrix V, size_t i, size_t j);
REAL d2Udx2(Config::geo geoConfig, Config::constants constantsConfig, Matrix U, size_t i, size_t j);
REAL d2Udy2(Config::geo geoConfig, Config::constants constantsConfig, Matrix U, size_t i, size_t j);

REAL dUVdx(Config::geo geoConfig, Config::constants constantsConfig, Matrix U, Matrix V, size_t i, size_t j);
REAL dV2dy(Config::geo geoConfig, Config::constants constantsConfig, Matrix V, size_t i, size_t j);
REAL d2Vdx2(Config::geo geoConfig, Config::constants constantsConfig, Matrix V, size_t i, size_t j);
REAL d2Vdy2(Config::geo geoConfig, Config::constants constantsConfig, Matrix V, size_t i, size_t j);

REAL dpdx(Config::geo geoConfig, Config::constants constantsConfig, Matrix P, size_t i, size_t j);
REAL dpdy(Config::geo geoConfig, Config::constants constantsConfig, Matrix P, size_t i, size_t j);