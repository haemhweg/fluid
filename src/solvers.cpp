#include <utility>
#include <iostream>
#include <cmath>

#include "solvers.h"
#include "real.h"
#include "config.h"
#include "matrix.h"
#include "differences.h"


// Computes the right hand side for the poisson problem
Matrix RHS_Poisson(const Config::geo geoConfig, const REAL delt, const Matrix& F, const Matrix& G)
{
  Matrix rhs(geoConfig.imax+1, geoConfig.jmax+1);

  for(unsigned i=1; i<geoConfig.imax+1; ++i){
    for(unsigned j=1; j<geoConfig.jmax+1; ++j){
      rhs.at(i,j) = (  (F.at(i,j) - F.at(i-1,j))/geoConfig.delx 
		     + (G.at(i,j) - G.at(i,j-1))/geoConfig.dely) / delt;
    }
  }

  return rhs;
}


// computes the residual
REAL comp_res(const Config::geo geoConfig, const Matrix& P, const Matrix& rhs)
{
  REAL res = 0;
  REAL delx = geoConfig.delx;
  REAL dely = geoConfig.dely;

  unsigned imax = geoConfig.imax;
  unsigned jmax = geoConfig.jmax;


  for(unsigned i=1; i<imax+1; ++i) {
    for(unsigned j=1; j<jmax+1; ++j) {
      res += std::pow(d2fdx(delx, P, i, j) + d2fdy(dely, P, i, j) - rhs.at(i,j), 2) / double(imax*jmax);
    }
  }

  return std::sqrt(res);
}

// Solves the discrete poisson equation using the SOR method
std::pair<unsigned, REAL> SOR_Poisson(const Config::geo geoConfig, const Config::solver solverConfig, 
				      Matrix& P, const Matrix& rhs)
{
  // comupte initial residual  
  REAL res = comp_res(geoConfig, P, rhs);
  
  REAL omega = solverConfig.omega;
  REAL delx2 = geoConfig.delx*geoConfig.delx;
  REAL dely2 = geoConfig.dely*geoConfig.dely;

  unsigned it = 0;

  for(it = 1; it<solverConfig.itmax && res>solverConfig.eps; ++it) { 
    // set zero dirichlet bc for P
    for(unsigned j=1; j<geoConfig.jmax+1; ++j){
      P.at(0,j) = P.at(1,j);
      P.at(geoConfig.imax+1,j) = P.at(geoConfig.imax,j); 
    }
    for(unsigned i=1; i<geoConfig.imax+1; ++i){
      P.at(i,0) = P.at(i,1);
      P.at(i,geoConfig.jmax+1) = P.at(i,geoConfig.jmax); 
    }
   
    // compute next iteration
    for(unsigned i=1; i<geoConfig.imax+1; ++i) {
      for(unsigned j=1; j<geoConfig.jmax+1; ++j) {
	P.at(i,j) = (1-omega)*P.at(i,j) + omega/(2*(1/delx2 + 1/dely2))
	  *((P.at(i+1,j)+P.at(i-1,j))/delx2 + (P.at(i,j+1)+P.at(i,j-1))/dely2
	    - rhs.at(i,j));
      }
    }

    // compute residual
    res = comp_res(geoConfig, P, rhs);

    //std::cout << res << std::endl;
  }

  return std::make_pair(it, res);
}
