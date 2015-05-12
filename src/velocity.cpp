#include <utility>
#include <algorithm>
#include <iostream>
#include <string>

#include "velocity.h"
#include "config.h"
#include "matrix.h"
#include "real.h"
#include "differences.h"
#include "geometry.h"


void Velocity::updateBoundary()
{
  unsigned jmax = geoConfig.jmax;
  unsigned imax = geoConfig.imax;

  // Setze Randwerte bzgl. der Hindernisse
  auto boundary = geometry.get_boundary();

  for(const auto& cell : boundary){
    unsigned i = cell.first, j = cell.second;
    switch (geometry.at(i, j))
      {
  	// Randkanten
      case B_N:
  	U.at(i,j) = -U.at(i,j+1);
  	U.at(i-1,j) = -U.at(i-1,j+1);
  	V.at(i,j) = 0;
  	break;
      case B_W:
  	U.at(i,j) = 0;
  	V.at(i,j) = -V.at(i-1,j);
  	V.at(i,j-1) = -V.at(i-1,j-1);
  	break;
      case B_S:
  	U.at(i,j) = -U.at(i,j-1);
  	U.at(i-1,j) = -U.at(i-1,j-1);
  	V.at(i,j-1) = 0;
  	break;
      case B_O:
  	U.at(i-1,j) = 0;
  	V.at(i,j) = -V.at(i+1,j);
  	V.at(i,j-1) = -V.at(i+1,j-1);
  	break;

  	// Randkanten
      case B_NO:
  	U.at(i,j) = 0;
  	U.at(i-1,j) = -U.at(i-1,j+1);
  	V.at(i,j) = 0;
  	V.at(i,j-1) = -V.at(i+1,j-1);
  	break;
      case B_SO:
  	U.at(i,j) = 0;
  	U.at(i-1,j) = -U.at(i-1,j-1);
  	V.at(i,j) = -V.at(i+1,j);
  	V.at(i,j-1) = 0;
  	break;
      case B_SW:
  	U.at(i,j) = -U.at(i,j-1);
  	U.at(i-1,j) = 0;
  	V.at(i,j) = -V.at(i-1,j);
  	V.at(i,j-1) = 0;
  	break;
      case B_NW:
  	U.at(i,j) = -U.at(i,j+1);
  	U.at(i-1,j) = 0;
  	V.at(i,j) = 0;
  	V.at(i,j-1) = -V.at(i-1,j-1);
  	break;      
      }
  }

  // Setze Randwerte bzgl. den Randbedingungen
  // top boundary, i.e. j = jmax,jmax+1
  switch (bc.wt)
    {
    case NO_SLIP:
      for(unsigned i=1; i<imax+1; ++i){
	U.at(i,jmax+1) = -U.at(i,jmax);
	V.at(i,jmax) = 0;
      }
      break;
    case FREE_SLIP:
      for(unsigned i=1; i<imax+1; ++i){
	U.at(i,jmax+1) = U.at(i,jmax);
	V.at(i,jmax) = 0;
      }
      break;
    case OUTFLOW:
      for(unsigned i=1; i<imax+1; ++i){
	U.at(i,jmax+1) = U.at(i,jmax);
	V.at(i,jmax) = V.at(i,jmax-1);
      }
      break;
    }

  // right boundary, i.e. i = imax,imax+1
  switch (bc.wr)
    {
    case NO_SLIP:
      for(unsigned j=1; j<jmax+1; ++j){
	U.at(imax,j) = 0;
	V.at(imax+1,j) = -V.at(imax,j);
      }
      break;
    case FREE_SLIP:
      for(unsigned j=1; j<jmax+1; ++j){
	U.at(imax,j) = 0;
	V.at(imax+1,j) = V.at(imax,j);	
      }
      break;
    case OUTFLOW:
      for(unsigned j=1; j<jmax+1; ++j){
	U.at(imax,j) = U.at(imax-1,j);
	V.at(imax+1,j) = V.at(imax,j);
      }
      break;
    }

  // bottom boundary, i.e. j=0
  switch (bc.wb)
    {
    case NO_SLIP:
      for(unsigned i=1; i<imax+1; ++i){
	U.at(i,0) = -U.at(i,1);
	V.at(i,0) = 0;
      }
      break;
    case FREE_SLIP:
      for(unsigned i=1; i<imax+1; ++i){
	U.at(i,0) = U.at(i,1);
	V.at(i,0) = 0;
      }
      break;
    case OUTFLOW:
      for(unsigned i=1; i<imax+1; ++i){
	U.at(i,0) = U.at(i,1);
	V.at(i,0) = V.at(i,1);
      }
      break;
    }

  // left boundary, i.e. i=0
  switch (bc.wl)
    {
    case NO_SLIP:
      for(unsigned j=1; j<jmax+1; ++j){
	U.at(0,j) = 0;
	V.at(0,j) = -V.at(1,j);
      }
      break;
    case FREE_SLIP:
      for(unsigned j=1; j<jmax+1; ++j){
	U.at(0,j) = U.at(1,j);
	V.at(0,j) = 0;
      }
      break;
    case OUTFLOW:
      for(unsigned j=1; j<jmax+1; ++j){
	U.at(0,j) = U.at(1,j);
	V.at(0,j) = V.at(1,j);
      }
      break;
    }
}

void Velocity::updateIntermidiate(const REAL delt)
{
  const REAL delx = geoConfig.delx;
  const REAL dely = geoConfig.dely;
  const REAL alpha = constantsConfig.alpha;
  const REAL Re = constantsConfig.Re;
  const REAL gx = constantsConfig.GX;
  const REAL gy = constantsConfig.GY;

  auto fluid = geometry.get_fluid();

  for(const auto& cell : fluid){
    unsigned i = cell.first, j = cell.second;

    if(i<geoConfig.imax){
      F.at(i,j) = U.at(i,j) + delt*((d2fdx(delx, U, i, j) + d2fdy(dely, U, i, j))/Re - df2dx(delx, alpha, U, i, j)
				    - dfgdy(dely, alpha, U, V, i, j) + gx);
    }
    if(j<geoConfig.jmax){
      G.at(i,j) = V.at(i,j) + delt*((d2fdx(delx, V, i, j) + d2fdy(dely, V, i, j))/Re - df2dy(dely, alpha, V, i, j)
				    - dfgdx(delx, alpha, U, V, i, j) + gy);
    }
  }

  // Set boundary values of F and G
  auto boundary = geometry.get_boundary();
  for(const auto& cell : boundary){
    unsigned i = cell.first, j = cell.second;
    switch (geometry.at(i,j))
      {
	// Randkanten
      case B_N:
	G.at(i,j) = V.at(i,j);
	break;
      case B_O:
	F.at(i,j) = U.at(i,j);
	break;
      case B_S:
	G.at(i,j-1) = V.at(i,j-1);
	break;
      case B_W:
	F.at(i-1,j) = U.at(i-1,j);
	break;

	// Randecken
      case B_NO:
	F.at(i,j) = U.at(i,j);
	G.at(i,j) = V.at(i,j);
	break;
      case B_SO:
	F.at(i,j) = U.at(i,j);
	G.at(i,j-1) = V.at(i,j-1);
	break;
      case B_SW:
	F.at(i-1,j) = U.at(i-1,j);
	G.at(i,j-1) = V.at(i,j-1);
	break;
      case B_NW:
	F.at(i-1,j) = U.at(i-1,j);
	G.at(i,j) = V.at(i,j);
	break;
      }
  }
}


const Matrix Velocity::getDivergenceIntermidiate(const REAL delt)
{
  Matrix divergence{geoConfig.imax+1, geoConfig.jmax+1};

  updateBoundary();

  updateSPBoundary(geoConfig.imax, geoConfig.jmax, U, V);
  
  updateIntermidiate(delt);
  
  for(unsigned i=1; i<geoConfig.imax+1; ++i){
    for(unsigned j=1; j<geoConfig.jmax+1; ++j){
      divergence.at(i,j) = (  (F.at(i,j) - F.at(i-1,j))/geoConfig.delx 
		       + (G.at(i,j) - G.at(i,j-1))/geoConfig.dely) / delt;
    }
  }

  return divergence;
}

void Velocity::update(const REAL delt, const Matrix& P)
{
  const REAL delx = geoConfig.delx;
  const REAL dely = geoConfig.dely;
  
  auto fluid = geometry.get_fluid();

  for(const auto cells : fluid){
    unsigned i = cells.first, j = cells.second;

    if(i<geoConfig.imax){
      U.at(i,j) = F.at(i,j) - delt/delx * (P.at(i+1,j) - P.at(i,j));
    }
    if(j<geoConfig.imax){
      V.at(i,j) = G.at(i,j) - delt/dely * (P.at(i,j+1) - P.at(i,j));
    }
  }
}
