#include <utility>
#include <iostream>
#include <cmath>
#include "mpi.h"

#include "solvers.h"
#include "real.h"
#include "config.h"
#include "matrix.h"
#include "differences.h"
#include "geometry.h"
#include "parallel.h"

// computes the residual
REAL comp_res(const MPI_Comm & comm_grid, const Config::geo& geoConfig, const Geometry& geometry, const Matrix& P, const Matrix& rhs)
{
  REAL res_loc = 0, res_glob = 0;
  REAL delx = geoConfig.delx;
  REAL dely = geoConfig.dely;

  const auto& fluid = geometry.get_fluid();
  const unsigned nfluid_loc = fluid.size();
  unsigned nfluid_glob = 0;

  MPI_Allreduce((void*)&nfluid_loc, (void*)&nfluid_glob, 1, MPI_UNSIGNED, MPI_SUM, comm_grid);

  for(const auto& cell : fluid){
    unsigned i = cell.first, j = cell.second;  
    res_loc += std::pow(d2fdx(delx, P, i, j) + d2fdy(dely, P, i, j) - rhs.at(i,j), 2) / double(nfluid_glob);  
  }

  MPI_Allreduce((void*)&res_loc, (void*)&res_glob, 1, MPI_DOUBLE, MPI_SUM, comm_grid);

  return std::sqrt(res_glob);
}

void updatePressureBoundary(const MPI_Comm& comm_grid, const Config::geo& geoConfig, const Geometry& geometry, Matrix& P)
{  
  REAL delx2 = geoConfig.delx*geoConfig.delx;
  REAL dely2 = geoConfig.dely*geoConfig.dely;

  const auto& boundary = geometry.get_boundary();
  const auto coords = get_MPI_Cart_coords(comm_grid, 2);
  const auto dims = get_MPI_Dims_create(MPI_COMM_WORLD, 2);

  // set zero dirichlet bc if P is at the boundary
  if(coords[0]==0){
    for(unsigned j=1; j<geoConfig.jmax+1; ++j){
      P.at(0,j) = P.at(1,j);
    }    
  }else if(coords[0]==dims[0]-1){
    for(unsigned j=1; j<geoConfig.jmax+1; ++j){
      P.at(geoConfig.imax+1,j) = P.at(geoConfig.imax,j); 
    }
  }
  
  if(coords[1]==0){
    for(unsigned i=1; i<geoConfig.imax+1; ++i){
      P.at(i,0) = P.at(i,1);
    }
  }else if(coords[1]==dims[1]-1){
    for(unsigned i=1; i<geoConfig.imax+1; ++i){
      P.at(i,geoConfig.jmax+1) = P.at(i,geoConfig.jmax); 
    }
  }

  // Hindernis
  // for(const auto& cell : boundary){
  //   unsigned i = cell.first, j = cell.second;
  //   switch (geometry.at(i, j))
  //     {
  // 	// Randkanten
  //     case B_N:
  // 	P.at(i,j) = P.at(i,j+1);
  // 	break;
  //     case B_W:
  // 	P.at(i,j) = P.at(i-1,j);
  // 	break;
  //     case B_S:
  // 	P.at(i,j) = P.at(i,j-1);
  // 	break;
  //     case B_O:
  // 	P.at(i,j) = P.at(i+1,j);
  // 	break;

  // 	// Randecken
  //     case B_NO:
  // 	P.at(i,j) = (delx2*P.at(i,j+1)+dely2*P.at(i+1,j))/(delx2+dely2);
  // 	break;
  //     case B_SO:
  // 	P.at(i,j) = (delx2*P.at(i,j-1)+dely2*P.at(i+1,j))/(delx2+dely2);
  // 	break;
  //     case B_SW:
  // 	P.at(i,j) = (delx2*P.at(i,j-1)+dely2*P.at(i-1,j))/(delx2+dely2);
  // 	break;
  //     case B_NW:
  // 	P.at(i,j) = (delx2*P.at(i,j+1)+dely2*P.at(i-1,j))/(delx2+dely2);
  // 	break;      
  //     }
  // }
}

// Solves the discrete poisson equation using the SOR method
std::pair<unsigned, REAL> SOR_Poisson(const MPI_Comm& comm_grid, const Config::geo& geoConfig, const Config::solver solverConfig, 
				      const Geometry& geometry, Matrix& P, const Matrix& rhs)
{
  // comupte initial residual  
  REAL res = comp_res(comm_grid, geoConfig, geometry, P, rhs);
  
  const REAL omega = solverConfig.omega;
  const REAL delx2 = geoConfig.delx*geoConfig.delx;
  const REAL dely2 = geoConfig.dely*geoConfig.dely;
  const unsigned imax = geoConfig.imax;
  const unsigned jmax = geoConfig.jmax;

  const auto& fluid = geometry.get_fluid();

  unsigned it = 0, rank = get_MPI_Comm_rank(comm_grid);

  for(it = 1; it<solverConfig.itmax && res>solverConfig.eps; ++it) { 
    updatePressureBoundary(comm_grid, geoConfig, geometry, P); 

    // compute next iteration
    for(unsigned i=1; i<imax+1; ++i){
      for(unsigned j=1; j<jmax+1; ++j){
	P.at(i,j) = (1-omega)*P.at(i,j) + omega/(2*(1/delx2 + 1/dely2))
	  *((P.at(i+1,j)+P.at(i-1,j))/delx2 + (P.at(i,j+1)+P.at(i,j-1))/dely2
	    - rhs.at(i,j));
      }
    }

    // for(const auto& cell : fluid){
    //   unsigned i=cell.first, j=cell.second;

    //   P.at(i,j) = (1-omega)*P.at(i,j) + omega/(2*(1/delx2 + 1/dely2))
    // 	*((P.at(i+1,j)+P.at(i-1,j))/delx2 + (P.at(i,j+1)+P.at(i,j-1))/dely2
    // 	  - rhs.at(i,j));
    // }

    //MPI_Barrier(comm_grid);
    MPI_Matrix_exchange(comm_grid, P);
    //MPI_Barrier(comm_grid);   

    // compute residual
    res = comp_res(comm_grid, geoConfig, geometry, P, rhs);
    
    if(rank==0)
      std::cout << it << ": " << res << std::endl;
  }

  return std::make_pair(it, res);
}
