#include <utility>
#include <iostream>
#include <cmath>

#include "solvers.h"
#include "real.h"
#include "config.h"
#include "matrix.h"
#include "differences.h"
#include "geometry.h"



// computes the residual
REAL comp_res(const Config::geo& geoConfig, const Geometry& geometry, const Matrix& P, const Matrix& rhs)
{
  REAL res = 0;
  const REAL delx = geoConfig.delx;
  const REAL dely = geoConfig.dely;

  const auto& fluid = geometry.get_fluid();
  const unsigned nfluid = fluid.size();

  unsigned i=0, j=0;
  for(const auto& cell : fluid){
    i = cell.first; j = cell.second;  
    res += std::pow(d2fdx(delx, P, i, j) + d2fdy(dely, P, i, j) - rhs.at(i,j), 2) / double(nfluid);  
  }

  return std::sqrt(res);
}

void updatePressureBoundary(const Config::geo& geoConfig, const Geometry& geometry, Matrix& P)
{  
  const REAL delx2 = geoConfig.delx*geoConfig.delx;
  const REAL dely2 = geoConfig.dely*geoConfig.dely;
  const unsigned imax = geoConfig.imax;
  unsigned jmax = geoConfig.jmax, i=0, j=0;

  const auto& boundary = geometry.get_boundary();

  // set zero dirichlet bc on the square for P
  for(j=1; j<jmax+1; ++j){
    P.at(0,j) = P.at(1,j);
    P.at(imax+1,j) = P.at(imax,j); 
  }
  for(i=1; i<imax+1; ++i){
    P.at(i,0) = P.at(i,1);
    P.at(i,jmax+1) = P.at(i,jmax); 
  }

  // Hindernis
  for(const auto& cell : boundary){
    i = cell.first; j = cell.second;
    switch (geometry.at(i, j))
      {
  	// Randkanten
      case B_N:
	P.at(i,j) = P.at(i,j+1);
  	break;
      case B_W:
	P.at(i,j) = P.at(i-1,j);
  	break;
      case B_S:
	P.at(i,j) = P.at(i,j-1);
  	break;
      case B_O:
	P.at(i,j) = P.at(i+1,j);
  	break;

  	// Randkanten
      case B_NO:
	P.at(i,j) = (delx2*P.at(i,j+1)+dely2*P.at(i+1,j))/(delx2+dely2);
  	break;
      case B_SO:
	P.at(i,j) = (delx2*P.at(i,j-1)+dely2*P.at(i+1,j))/(delx2+dely2);
  	break;
      case B_SW:
	P.at(i,j) = (delx2*P.at(i,j-1)+dely2*P.at(i-1,j))/(delx2+dely2);
  	break;
      case B_NW:
	P.at(i,j) = (delx2*P.at(i,j+1)+dely2*P.at(i-1,j))/(delx2+dely2);
  	break;      
      }
  }
}

// Solves the discrete poisson equation using the SOR method
std::pair<unsigned, REAL> SOR_Poisson(const Config::geo& geoConfig, const Config::solver solverConfig, const Geometry& geometry, 
				      Matrix& P, const Matrix& rhs)
{
  // comupte initial residual  
  REAL res = comp_res(geoConfig, geometry, P, rhs);
  
  const REAL omega = solverConfig.omega;
  const REAL delx2 = geoConfig.delx*geoConfig.delx;
  const REAL dely2 = geoConfig.dely*geoConfig.dely;

  const auto& fluid = geometry.get_fluid();

  unsigned it = 0, i=0, j=0;

  for(it = 1; it<solverConfig.itmax && res>solverConfig.eps; ++it) { 
    updatePressureBoundary(geoConfig, geometry, P);    
   
    // compute next iteration
    for(const auto& cell : fluid){
      i=cell.first; j=cell.second;

      P.at(i,j) = (1-omega)*P.at(i,j) + omega/(2*(1/delx2 + 1/dely2))
	*((P.at(i+1,j)+P.at(i-1,j))/delx2 + (P.at(i,j+1)+P.at(i,j-1))/dely2
	  - rhs.at(i,j));
    }

    // compute residual
    res = comp_res(geoConfig, geometry, P, rhs);
  }

  return std::make_pair(it, res);
}
