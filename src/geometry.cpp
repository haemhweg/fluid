#include <vector>
#include <functional>
#include <utility>
#include <iostream>

#include "geometry.h"

Geometry::Geometry(const Config::geo& geoConfig, const init_geometry_fct& geoInit)
  : cells(geoInit(geoConfig)), imax(geoConfig.imax), jmax(geoConfig.jmax), cells_boundary(), cells_fluid()
{
  for(unsigned j=1; j<jmax+1; ++j){ 
    for(unsigned i=1; i<imax+1; ++i){
      if(at(i,j)==FLUID){
	cells_fluid.emplace_back(std::make_pair(i,j));
      }
      else if(at(i,j)!=BLOCK){
	cells_boundary.emplace_back(std::make_pair(i,j));
      }
    }
  }
}

void Geometry::print()
{
  for (int j = jmax+1; j >=0 ; --j){
    for (unsigned i = 0; i < imax+2; ++i){
      std::cout << int(at(i,j)) << "  ";
    }
    std::cout << std::endl;
  }
}
