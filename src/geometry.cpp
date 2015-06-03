#include <vector>
#include <functional>
#include <utility>
#include <iostream>

#include "geometry.h"


Geometry::Geometry(const Config::geo& geoConfig, const std::function<std::vector<CELL>(const Config::geo&)> fill_geometry)
  : cells(fill_geometry(geoConfig)), imax(geoConfig.imax), jmax(geoConfig.jmax), cells_boundary()
{
  for(unsigned i=0; i<imax+2; ++i){
    for(unsigned j=0; j<jmax+2; ++j){ 
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
  for (int j = jmax+1; j >=0 ; --j)
    {
      for (unsigned i = 0; i < imax+2; ++i)
	{
	  std::cout << int(at(i,j)) << "\t";
	}
      std::cout << std::endl;
    }
}
