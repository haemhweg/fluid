#include "real.h"
#include "config.h"
#include "matrix.h"

REAL dU2dx(Config::geo geoConfig, Config::constants constantsConfig, Matrix U, unsigned i, unsigned j);
REAL dUVdy(Config::geo geoConfig, Config::constants constantsConfig, Matrix U, Matrix V, unsigned i, unsigned j);
REAL d2Udx2(Config::geo geoConfig, Config::constants constantsConfig, Matrix U, unsigned i, unsigned j);
REAL d2Udy2(Config::geo geoConfig, Config::constants constantsConfig, Matrix U, unsigned i, unsigned j);

REAL dUVdx(Config::geo geoConfig, Config::constants constantsConfig, Matrix U, Matrix V, unsigned i, unsigned j);
REAL dV2dy(Config::geo geoConfig, Config::constants constantsConfig, Matrix V, unsigned i, unsigned j);
REAL d2Vdx2(Config::geo geoConfig, Config::constants constantsConfig, Matrix V, unsigned i, unsigned j);
REAL d2Vdy2(Config::geo geoConfig, Config::constants constantsConfig, Matrix V, unsigned i, unsigned j);

REAL dpdx(Config::geo geoConfig, Config::constants constantsConfig, Matrix P, unsigned i, unsigned j);
REAL dpdy(Config::geo geoConfig, Config::constants constantsConfig, Matrix P, unsigned i, unsigned j);
