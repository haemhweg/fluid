#include <vector>
#include <functional>
#include <utility>
#include <iostream>
#include "mpi.h"

#include "geometry.h"
#include "parallel.h"

Geometry::Geometry(const MPI_Comm& comm_grid, const Config::geo& geoConfig, const init_geometry_fct& geoInit)
  : cells(geoInit(geoConfig)), imax(geoConfig.imax), jmax(geoConfig.jmax), cells_boundary(), cells_fluid()
{
  auto coords = get_MPI_Cart_coords(comm_grid, 2);
  auto dims = get_MPI_Dims_create(comm_grid, 2);
  for(unsigned i=1; i<imax+1; ++i){
    for(unsigned j=1; j<jmax+1; ++j){ 
      if(at(i,j)==FLUID){
	cells_fluid.emplace_back(std::make_pair(i,j));
      }
      else if(at(i,j)!=BLOCK && 
	      ( coords[0]==0 || coords[1]==0 || coords[0]==dims[0]-1 || coords[1]==dims[1]-1 )){
	cells_boundary.emplace_back(std::make_pair(i,j));
      }
    }
  }
}

void Geometry::print()
{
  for (int j = jmax+1; j >=0 ; --j)
    {
      for (unsigned i = 0; i < imax+2; ++i)
	{
	  std::cout << int(at(i,j)) << "  ";
	}
      std::cout << std::endl;
    }
}
