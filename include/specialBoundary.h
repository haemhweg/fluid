#ifndef SPECIALBOUNDARY
#define SPECIALBOUNDARY

#include <functional>
#include <vector>

#include "geometry.h"
#include "matrix.h"

using special_boundary_fct = std::function<void(const unsigned, const unsigned, Matrix&, Matrix&)>;

void bc_DRIVEN_CAVITY(const unsigned imax, const unsigned jmax, Matrix& U, Matrix&); 

void bc_STEP(const unsigned, const unsigned jmax, Matrix& U, Matrix&);

void bc_KARMAN(const unsigned, const unsigned jmax, Matrix& U, Matrix&);

std::vector<CELL> geometry_DRIVEN_CAVITY(const Config::geo& geoConfig);

std::vector<CELL> geometry_STEP(const Config::geo& geoConfig);

std::vector<CELL> geometry_KARMAN(const Config::geo& geoConfig);

#endif
