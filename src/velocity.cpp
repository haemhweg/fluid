
#include "velocity.h"
#include "config.h"
#include "matrix.h"
#include "real.h"
#include "differences.h"

void compIntermediateVelocity(const Config::geo geoConfig, const Config::constants constantsConfig,
			      const REAL delt, const Matrix& U, const Matrix& V, Matrix& F, Matrix& G)
{
  const REAL delx = geoConfig.delx;
  const REAL dely = geoConfig.dely;
  const REAL alpha = constantsConfig.alpha;
  const REAL Re = constantsConfig.Re;
  const REAL gx = constantsConfig.GX;
  const REAL gy = constantsConfig.GY;

  
  // Set inner values of F
  for(unsigned i=1; i<geoConfig.imax; ++i){
    for(unsigned j=1; j<geoConfig.jmax+1; ++j){
      F.at(i,j) = U.at(i,j) + delt*((d2fdx(delx, U, i, j) + d2fdy(dely, U, i, j))/Re - df2dx(delx, alpha, U, i, j)
				    - dfgdy(dely, alpha, U, V, i, j) + gx);
    }
  }
  // Set inner values of G
  for(unsigned i=1; i<geoConfig.imax+1; ++i){
    for(unsigned j=1; j<geoConfig.jmax; ++j){
      G.at(i,j) = V.at(i,j) + delt*((d2fdx(delx, V, i, j) + d2fdy(dely, V, i, j))/Re - df2dy(dely, alpha, V, i, j)
				    - dfgdx(delx, alpha, U, V, i, j) + gy);
    }
  }

  // Set boundary values for F,G
  for(unsigned j=1; j<geoConfig.jmax; ++j){
    F.at(0,j) = U.at(0,j);
    F.at(geoConfig.imax,j) = U.at(geoConfig.imax,j);
  }
  for(unsigned i=1; i<geoConfig.imax; ++i){
    G.at(i,0) = V.at(i,0);
    G.at(i,geoConfig.jmax) = V.at(i,geoConfig.imax);
  }
  
}

void compNewVelocity(const Config::geo geoConfig, const REAL delt, Matrix& U, Matrix& V,
		     const Matrix& F, const Matrix& G, const Matrix& P)
{

  const REAL delx = geoConfig.delx;
  const REAL dely = geoConfig.dely;
  
  // Set inner values of U
  for(unsigned i=1; i<geoConfig.imax; ++i){
    for(unsigned j=1; j<geoConfig.jmax+1; ++j){
      U.at(i,j) = F.at(i,j) - delt/delx * (P.at(i+1,j) - P.at(i,j));
    }
  }

  // Set inner values of V
  for(unsigned i=1; i<geoConfig.imax+1; ++i){
    for(unsigned j=1; j<geoConfig.jmax; ++j){
      V.at(i,j) = G.at(i,j) - delt/dely * (P.at(i,j+1) - P.at(i,j));
    }
  }

  // Set boundary on vertical boundary
  for(unsigned j=1; j<geoConfig.jmax+1; ++j){
    // Invoke zero dirichlet bc on U
    U.at(0,j) = U.at(geoConfig.imax,j) = 0;

    // Invoke zero dirichlet bc on V by linear interpolation
    V.at(0,j) = -V.at(1,j);
    V.at(geoConfig.imax+1,j) = -V.at(geoConfig.imax,j);
  }
  // Set boundary on horizontal boundary
  for(unsigned i=0; i<geoConfig.imax+1; ++i){
    // Invoke zero dirichlet bc on U
    V.at(i,0) = V.at(i,geoConfig.jmax) = 0;

    // Invoke zero dirichlet bc on V by linear interpolation
    U.at(i,0) = -U.at(i,1);
    U.at(i,geoConfig.jmax+1) = -U.at(i,geoConfig.jmax);
  }

}
