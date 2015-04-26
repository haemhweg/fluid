#include <utility>
#include <iostream>
#include <cmath>

#include "solvers.h"
#include "real.h"
#include "config.h"
#include "matrix.h"


// computes the residual
REAL comp_res(const Matrix& P, const Matrix& rhs, const Config::geo geoConfig)
{
  REAL res = 0;
  REAL delx2 = geoConfig.delx*geoConfig.delx;
  REAL dely2 = geoConfig.dely*geoConfig.dely;

  unsigned imax = geoConfig.imax;
  unsigned jmax = geoConfig.jmax;


  for(unsigned i=1; i<imax+1; ++i) {
    for(unsigned j=1; j<jmax+1; ++j) {
      res += std::pow((P.at(i+1,j)-2*P.at(i,j)+P.at(i-1,j))/delx2 + (P.at(i,j+1)-2*P.at(i,j)+P.at(i,j-1))/dely2
		      - rhs.at(i,j), 2) / double(imax*jmax);
    }
  }

  return std::sqrt(res);
}

std::pair<unsigned, REAL> SOR_Poisson(Matrix& P, const Matrix& rhs, const Config::geo geoConfig, 
				    const Config::solver solverConfig)
{
  // comupte initial residual  
  REAL res = comp_res(P, rhs, geoConfig);
  
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
	  *((P.at(i+1,j)+P.at(i-1,j))/delx2 + (P.at(i,j+1)+P.at(i,j-1))/delx2
	    - rhs.at(i,j));
      }
    }

    // compute residual
    res = comp_res(P, rhs, geoConfig);

    std::cout << "Schritt " << it << ": " << res << std::endl;
  }

  return std::make_pair(it, res);
}