#include <vector>
#include <iostream>
#include <cassert>

#include "specialBoundary.h"
#include "geometry.h"
#include "real.h"
#include "field.h"

void bc_DRIVEN_CAVITY(const unsigned imax, const unsigned jmax, Field2D& U, Field2D&)
{
  for(unsigned i=1; i<imax+1; ++i){
    U.at(i,jmax+1) = REAL(2.) - U.at(i,jmax);
  }
}

std::vector<CELL> geometry_DRIVEN_CAVITY(const Config::geo& geoConfig)
{
  const unsigned imax = geoConfig.imax, jmax = geoConfig.jmax;
  
  std::vector<CELL> cells((imax+2)*(jmax+2), FLUID);
  
  for(unsigned i=1; i<imax+1; ++i){
    cells[i] = B_N;
    cells[(jmax+1)*(imax+2)+i] = B_S;
  }
  for(unsigned j=1; j<jmax+1; ++j){
    cells[j*(imax+2)] = B_O;
    cells[(j+1)*(imax+2)-1] = B_W;
  }

  cells[0] = cells[(imax+2)*(jmax+2)-1] = cells[imax+1] = cells[(jmax+1)*(imax+2)] = BLOCK;

  return cells;
}

void bc_STEP(const unsigned, const unsigned jmax, Field2D& U, Field2D&)
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
  const unsigned jmax = geoConfig.jmax;
  const unsigned imax = geoConfig.imax;
  assert(imax>jmax/2);

  std::vector<CELL> cells(geometry_DRIVEN_CAVITY(geoConfig));

  for (int j = jmax+1; j >=0 ; --j){
    for (unsigned i = 0; i < imax+2; ++i){
      std::cout << int(cells[j*(imax+2)+i]) << "  ";
    }
    std::cout << std::endl;
  }

  std::cout << std::endl;

  for(unsigned i=1; i<jmax/2; ++i){
    cells[(jmax/2)*(imax+2)+i] = B_N;
  }
  for(unsigned j=1; j<jmax/2; ++j){
    cells[j*(imax+2)+jmax/2] = B_O;
  }

  for(unsigned j=0; j<jmax/2; ++j){
    for(unsigned i=0; i<jmax/2; ++ i){
      cells[j*(imax+2)+i] = BLOCK;
    }    
  }

  cells[(jmax/2)*(imax+2)+jmax/2] = B_NO;
  cells[(jmax/2)*(imax+2)] = BLOCK;
  cells[jmax/2] = BLOCK;
 
  return cells;
}

void bc_KARMAN(const unsigned, const unsigned jmax, Field2D& U, Field2D&)
{
  for(unsigned j=1; j<jmax+1; ++j){
    U.at(0,j) = REAL(1.);
  }  
}

std::vector<CELL> geometry_KARMAN(const Config::geo& geoConfig)
{  
  const unsigned jmax = geoConfig.jmax;
  const unsigned offset = 2*jmax/5+(jmax%5)/2;
  unsigned i, j;

  std::vector<CELL> cells(geometry_DRIVEN_CAVITY(geoConfig));

  // Untere 2 Kanten des Balken
  cells[(offset)*(jmax+2)+offset+1] = B_SW;
  cells[(offset+1)*(jmax+2)+offset+1] = B_SO;

  cells[(offset)*(jmax+2)+offset+2] = B_NW;
  cells[(offset+1)*(jmax+2)+offset+2] = BLOCK;
  cells[(offset+2)*(jmax+2)+offset+2] = B_O;

  // Innerer Teil des Balken
  for(j=3, i=1; j<jmax/5-1; ++j, ++i){
    cells[(offset+i)*(jmax+2)+offset+j] = B_W;
    cells[(offset+i+1)*(jmax+2)+offset+j] = BLOCK;
    cells[(offset+i+2)*(jmax+2)+offset+j] = B_O;    
  }
  
  // Obersten 2 Kanten des Balken
  cells[(offset+i)*(jmax+2)+offset+jmax/5-1] = B_W;
  cells[(offset+i+1)*(jmax+2)+offset+jmax/5-1] = BLOCK;
  cells[(offset+i+2)*(jmax+2)+offset+jmax/5-1] = B_SO;
  
  cells[(offset+i+1)*(jmax+2)+offset+jmax/5] = B_NW;
  cells[(offset+i+2)*(jmax+2)+offset+jmax/5] = B_NO;
 
  return cells;
}
