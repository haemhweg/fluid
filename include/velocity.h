#ifndef VELOCITY_H
#define VELOCITY_H

#include <utility>
#include <string>

#include "real.h"
#include "config.h"
#include "matrix.h"
#include "specialBoundary.h"
#include "geometry.h"


class Velocity
{

  Matrix F, G;
  Matrix U, V;

  const Config::geo geoConfig;
  const Config::constants constantsConfig;
  const Config::boundaryCondition bc;

  const special_boundary& updateSPBoundary;

  const Geometry& geometry;
  
  /**
   *
   */
  void updateBoundary();
  
  /**
   *
   */
  void updateIntermidiate(const REAL delt);
  
public:
  
  Velocity(const Config::geo geo_, const Config::constants constants_, const Config::boundaryCondition bc_, 
	   const Geometry& geometry_, const special_boundary& bc_sp_) 
    : F(geo_.imax+1, geo_.jmax+1), G(geo_.imax+1, geo_.jmax+1), U(geo_.imax+2, geo_.jmax+2, constants_.UI), 
      V(geo_.imax+2, geo_.jmax+2, constants_.VI), geoConfig(geo_), constantsConfig(constants_), 
      bc(bc_), updateSPBoundary(bc_sp_), geometry(geometry_) { }

  /**
   *  @brief Computes the right hand side for the discrete 2D poisson equation.
   *  @param geoConfig Information for the discretized grid
   *  @param delt Time step size
   *  @param F, G Intermediate velocity computed by @compIntermediateVelocity()
   *  @return The iteration number and achieved accuracy of the solution
   */
  const Matrix getDivergenceIntermidiate(const REAL delt);    
  
  /**
   *  @brief Computes the velocity for the next time step using the given intermediate values F,G and the new pressure P
   *  @param geoConfig, delt Grid and time dicretization
   *  @param U,V Contains the new velocity afterwards
   *  @param F,G Intermediate values for the velocity, see @compIntermediateVelocity()
   *  @param P New pressure from @SOR_Poisson()
   */
  void update(const REAL delt, const Matrix& P);

  REAL getMaxU() const { return U.getMax(); }
  REAL getMaxV() const { return V.getMax(); }

  void writeVTK(const unsigned step) 
  { writeVectorFieldVTK("Velocity"+std::to_string(step)+".vtk", "Velocity", U, V, geoConfig.delx, geoConfig.dely); }
};

#endif
