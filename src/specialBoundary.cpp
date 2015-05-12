#include <vector>

#include "specialBoundary.h"
#include "geometry.h"
#include "real.h"
#include "matrix.h"

void bc_DRIVEN_CAVITY(const unsigned imax, const unsigned jmax, Matrix& U, Matrix&)
{
  for(unsigned i=1; i<imax+1; ++i){
    U.at(i,jmax+1) = REAL(2.) - U.at(i,jmax);
  }
}

std::vector<CELL> geometry_DEFAULT(const unsigned imax, const unsigned jmax)
{
  std::vector<CELL> cells((imax+2)*(jmax+2), FLUID);

  for(unsigned i=1; i<imax+1; ++i){
    cells[i*(jmax+2)] = B_N;
    cells[(i+1)*(jmax+2)-1] = B_S;
  }
  for(unsigned j=1; j<jmax+1; ++j){
    cells[j] = B_O;
    cells[(imax+1)*(jmax+2)+j] = B_W;
  }

  cells[0] = cells[jmax+1] = cells[(imax+1)*(jmax+2)] = cells[(imax+2)*(jmax+2)-1] = BLOCK;

  return cells;
}
