#include <vector>
#include <functional>
#include <utility>
#include <iostream>

#include "geometry.h"


Geometry::Geometry(const unsigned imax_, const unsigned jmax_, 
		   const std::function<std::vector<CELL>(unsigned,unsigned)> fill_geometry)
  : cells(fill_geometry(imax_,jmax_)), imax(imax_), jmax(jmax_), cells_boundary()
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
  for (unsigned i = 0; i < imax+2; i++)
    {
      for (unsigned j = 0; j < jmax+2; j++)
	{
	  std::cout << int(at(i,j)) << "\t";
	}
      std::cout << std::endl;
    }
}
