#include <vector>
#include <functional>
#include <utility>
#include <iostream>
#include <fstream>

#include "geometry.h"

Geometry::Geometry(const Config::geo& geoConfig, const init_geometry_fct& geoInit)
  : cells(geoInit(geoConfig)), imax_(geoConfig.imax), jmax_(geoConfig.jmax), cells_boundary(), cells_fluid()
{
  for(unsigned j=1; j<jmax_+1; ++j){ 
    for(unsigned i=1; i<imax_+1; ++i){
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
  for (int j = jmax_+1; j >=0 ; --j){
    for (unsigned i = 0; i < imax_+2; ++i){
      std::cout << int(at(i,j)) << "  ";
    }
    std::cout << std::endl;
  }
}

void Geometry::writeVTK(const std::string & filename, const std::string& descr, const double dx, const double dy) const
{
  std::ofstream fs(filename);

  fs << "# vtk DataFile Version 3.0\n"
     << "Geometry\n"
     << "ASCII\n"
     << "DATASET RECTILINEAR_GRID\n"
     << "DIMENSIONS " << imax_+2 << " " << jmax_+2 << " 1\n"
     << "X_COORDINATES " << imax_+2 << " double\n";

  for(unsigned i=0; i<imax_+2; ++i) {
    fs << dx*i << " ";
  }

  fs << "\nY_COORDINATES " << jmax_+2 << " double\n";
  
  for(unsigned j=0; j<jmax_+2; ++j) {
    fs << dy*j << " ";
  }

  fs << "\nZ_COORDINATES 1 double\n"
     << "0.0\n"
     << "POINT_DATA " << (imax_+2)*(jmax_+2) << "\n"
     << "SCALARS " << descr << " double 1\n"
     << "LOOKUP_TABLE default\n";

  
  for(const auto& cell : cells) {
    if(cell == BLOCK){
      fs << 0. << "\n";
    }else{
      fs << 1. << "\n";
    }
  }

  fs << std::endl;    
}
