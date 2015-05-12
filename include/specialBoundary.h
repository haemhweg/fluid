#ifndef SPECIALBOUNDARY
#define SPECIALBOUNDARY

#include <functional>
#include <vector>

#include "geometry.h"
#include "matrix.h"

using special_boundary = std::function<void(const unsigned,const unsigned,Matrix&,Matrix&)>;

//void bc_NONE(const unsigned, const unsigned, Matrix&, Matrix&) { }

void bc_DRIVEN_CAVITY(const unsigned imax, const unsigned jmax, Matrix& U, Matrix&);

std::vector<CELL> geometry_DEFAULT(const unsigned imax, const unsigned jmax);

#endif
