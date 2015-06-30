#include <vector>
#include "mpi.h"

#include "specialBoundary.h"
#include "geometry.h"
#include "real.h"
#include "matrix.h"
#include "parallel.h"

void bc_DRIVEN_CAVITY(const unsigned imax, const unsigned jmax, Matrix& U, Matrix&)
{
  for(unsigned i=1; i<imax+1; ++i){
    U.at(i,jmax+1) = REAL(2.) - U.at(i,jmax);
  }
}

std::vector<CELL> geometry_DRIVEN_CAVITY(const MPI_Comm& comm_grid, const Config::geo& geoConfig)
{
  const unsigned imax = geoConfig.imax, jmax = geoConfig.jmax;
  const auto coords = get_MPI_Cart_coords(comm_grid, 2);
  const auto dims = get_MPI_Dims_create(MPI_COMM_WORLD, 2);
  
  std::vector<CELL> cells((imax+2)*(jmax+2), FLUID);
  
  // for(unsigned i=1; i<imax+1; ++i){
  //   cells[i*(jmax+2)] = B_N;
  //   cells[(i+1)*(jmax+2)-1] = B_S;
  // }

  // for(unsigned j=1; j<jmax+1; ++j){
  //   cells[j] = B_O;
  //   cells[(imax+1)*(jmax+2)+j] = B_W;
  // }


  if(coords[0]==0){
    for(unsigned i=1; i<imax+1; ++i){
      cells[i*(jmax+2)] = B_N;
    }
  }else{
    for(unsigned i=1; i<imax+1; ++i){
      cells[i*(jmax+2)] = BLOCK;
    }
  }
  if(coords[0]==dims[0]-1){
    for(unsigned i=1; i<imax+1; ++i){
      cells[(i+1)*(jmax+2)-1] = B_S;
    }
  }else{
    for(unsigned i=1; i<imax+1; ++i){
      cells[(i+1)*(jmax+2)-1] = BLOCK;
    }
  }

  if(coords[1]==0){
    for(unsigned j=1; j<jmax+1; ++j){
      cells[j] = B_O;
    }
  }else{
    for(unsigned j=1; j<jmax+1; ++j){
      cells[j] = BLOCK;
    }
  }
  if(coords[1]==dims[1]-1){
    for(unsigned j=1; j<jmax+1; ++j){
      cells[(imax+1)*(jmax+2)+j] = B_W;
    }
  }else{
    for(unsigned j=1; j<jmax+1; ++j){
      cells[(imax+1)*(jmax+2)+j] = BLOCK;
    }
  }

  cells[0] = cells[jmax+1] = cells[(imax+1)*(jmax+2)] = cells[(imax+2)*(jmax+2)-1] = BLOCK;

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

std::vector<CELL> geometry_STEP(const MPI_Comm& comm_grid, const Config::geo& geoConfig)
{  
  const unsigned jmax = geoConfig.jmax;

  std::vector<CELL> cells(geometry_DRIVEN_CAVITY(comm_grid, geoConfig));

  for(unsigned j=1; j<jmax/2; ++j){
    cells[j*(jmax+2)+jmax/2] = B_N;
    cells[(jmax/2)*(jmax+2)+j] = B_O;
  }

  for(unsigned j=0; j<jmax/2; ++j){
    for(unsigned i=0; i<jmax/2; ++ i){
      cells[j+i*(jmax+2)] = BLOCK;
    }    
  }

  cells[(jmax/2)*(jmax+2)+jmax/2] = B_NO;
  cells[jmax/2] = BLOCK;
  cells[(jmax/2)*(jmax+2)] = BLOCK;
 
  return cells;
}

void bc_KARMAN(const unsigned, const unsigned jmax, Matrix& U, Matrix&)
{
  for(unsigned j=1; j<jmax+1; ++j){
    U.at(0,j) = REAL(1.);
  }  
}

std::vector<CELL> geometry_KARMAN(const MPI_Comm& comm_grid, const Config::geo& geoConfig)
{  
  const unsigned jmax = geoConfig.jmax;
  const unsigned offset = 2*jmax/5+(jmax%5)/2;
  unsigned i, j;

  std::vector<CELL> cells(geometry_DRIVEN_CAVITY(comm_grid, geoConfig));

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
