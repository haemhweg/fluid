#include "specialBoundary.h"
#include "real.h"
#include "matrix.h"

void bc_DRIVEN_CAVITY(const unsigned imax, const unsigned jmax, Matrix& U, Matrix&)
{
  for(unsigned i=1; i<imax+1; ++i){
    U.at(i,jmax+1) = REAL(2.) - U.at(i,jmax);
  }
}
