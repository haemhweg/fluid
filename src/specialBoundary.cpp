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

std::vector<CELL> geometry_DRIVEN_CAVITY(const Config::geo& geoConfig)
{
  const unsigned imax = geoConfig.imax, jmax = geoConfig.jmax;
  auto at = [imax] (const unsigned i, const unsigned j) -> unsigned { return j*(imax+2)+i; };
  
  std::vector<CELL> cells((imax+2)*(jmax+2), FLUID);

  for(unsigned i=1; i<imax+1; ++i){
    cells[at(i,0)] = B_N;
    cells[at(i,jmax+1)] = B_S;
  }
  for(unsigned j=1; j<jmax+1; ++j){
    cells[at(0,j)] = B_O;
    cells[at(imax+1,j)] = B_W;
  }

  cells[at(0,0)] = cells[at(imax+1,0)] = cells[at(0,jmax+1)] = cells[at(imax+1,jmax+1)] = BLOCK;

  return cells;
}

void bc_STEP(const unsigned, const unsigned jmax, Matrix& U, Matrix&)
{
  for(unsigned j=0; j<jmax/2+1; ++j){
    U.at(0,j) = 0;
  }
  for(unsigned j=jmax/2+1; j<jmax+1; ++j){
    U.at(0,j) = REAL(1.);
  }
}

std::vector<CELL> geometry_STEP(const Config::geo& geoConfig)
{  
  const unsigned imax = geoConfig.imax, jmax = geoConfig.jmax;
  auto at = [imax] (const unsigned i, const unsigned j) -> unsigned { return j*(imax+2)+i; };

  std::vector<CELL> cells(geometry_DRIVEN_CAVITY(geoConfig));

  for(unsigned j=1; j<jmax/2; ++j){
    cells[at(jmax/2,j)] = B_O;
    cells[at(j,jmax/2)] = B_N;
  }

  for(unsigned j=0; j<jmax/2; ++j){
    for(unsigned i=0; i<jmax/2; ++ i){
      cells[at(i,j)] = BLOCK;
    }    
  }

  cells[at(jmax/2,jmax/2)] = B_NO;
  cells[at(0,jmax/2)] = BLOCK;
  cells[at(jmax/2,0)] = BLOCK;
 
  return cells;
}

void bc_KARMAN(const unsigned, const unsigned jmax, Matrix& U, Matrix&)
{
  for(unsigned j=1; j<jmax+1; ++j){
    U.at(0,j) = REAL(1.);
  }  
}

std::vector<CELL> geometry_KARMAN(const Config::geo& geoConfig)
{  
  const unsigned imax = geoConfig.imax, jmax = geoConfig.jmax;
  const unsigned offset = 2*jmax/5+(jmax%5)/2;
  auto at = [imax] (const unsigned i, const unsigned j) -> unsigned { return j*(imax+2)+i; };
  unsigned j;

  std::vector<CELL> cells(geometry_DRIVEN_CAVITY(geoConfig));

  // Untere 2 Kanten des Balken
  cells[at(offset,offset)] = B_SW;
  cells[at(offset+1,offset)] = B_SO;

  cells[at(offset,offset+1)] = B_NW;
  cells[at(offset+1,offset+1)] = BLOCK;
  cells[at(offset+2,offset+1)] = B_O;

  // Innerer Teil des Balken
  for(j=2; j<jmax/5-1; ++j){
    cells[at(offset+j-1,offset+j)] = B_W;
    cells[at(offset+j,offset+j)] = BLOCK;
    cells[at(offset+j+1,offset+j)] = B_O;    
  }
  
  // Obersten 2 Kanten des Balken
  cells[at(offset+jmax/5-2,offset+jmax/5-1)] = B_W;
  cells[at(offset+jmax/5-1,offset+jmax/5-1)] = BLOCK;
  cells[at(offset+jmax/5,offset+jmax/5-1)] = B_SO;
  
  cells[at(offset+jmax/5-1,offset+jmax/5)] = B_NW;
  cells[at(offset+jmax/5,offset+jmax/5)] = B_NO;
 
  return cells;
}
