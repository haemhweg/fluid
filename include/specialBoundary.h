#ifndef SPECIALBOUNDARY
#define SPECIALBOUNDARY

#include <functional>
#include <vector>

#include "geometry.h"
#include "matrix.h"

typedef std::function<void(const unsigned,const unsigned,Matrix&,Matrix&)> special_boundary;

//void bc_NONE(const unsigned, const unsigned, Matrix&, Matrix&) { }

void bc_DRIVEN_CAVITY(const unsigned imax, const unsigned jmax, Matrix& U, Matrix&);

std::vector<CELL> geometry_DEFAULT(const Config::geo& geoConfig);

std::vector<CELL> geometry_STEP(const Config::geo& geoConfig);
#endif
