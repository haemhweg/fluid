#ifndef PROBLEMS
#define PROBLEMS

#include <functional>
#include <vector>

#include "matrix.h"
#include "config.h"

enum CELL { FLUID=0, BLOCK, B_N, B_O, B_S, B_W, B_NW, B_SW, B_SO, B_NO };

typedef std::function<void(const unsigned, const unsigned, Matrix&, Matrix&)> special_boundary_fct;

typedef std::function<std::vector<CELL>(const Config::geo& geoConfig)> init_geometry_fct;

void bc_DRIVEN_CAVITY(const unsigned imax, const unsigned jmax, Matrix& U, Matrix&); 

void bc_STEP(const unsigned, const unsigned jmax, Matrix& U, Matrix&);

void bc_KARMAN(const unsigned, const unsigned jmax, Matrix& U, Matrix&);

std::vector<CELL> geometry_DRIVEN_CAVITY(const Config::geo& geoConfig);

std::vector<CELL> geometry_STEP(const Config::geo& geoConfig);

std::vector<CELL> geometry_KARMAN(const Config::geo& geoConfig);

#endif
