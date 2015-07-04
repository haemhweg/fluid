#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <functional>

#include "config.h"
#include "specialBoundary.h"

class Geometry
{
  std::vector<CELL> cells;
  unsigned imax_, jmax_;
  
  std::vector<std::pair<unsigned,unsigned>> cells_boundary;
  std::vector<std::pair<unsigned,unsigned>> cells_fluid;

 public:

  Geometry(const Config::geo& geoConfig, const init_geometry_fct& geoInit);

  CELL at(unsigned i, unsigned j) const { return cells[j*(imax_+2)+i]; }
  CELL& at(unsigned i, unsigned j) { return cells[j*(imax_+2)+i]; }

  const std::vector<std::pair<unsigned,unsigned>>& get_boundary() const { return cells_boundary; }
  const std::vector<std::pair<unsigned,unsigned>>& get_fluid() const { return cells_fluid; }

  void print();
  void writeVTK(const std::string & filename, const std::string& descr, 
		const double dx=0.1, const double dy=0.1) const;
};

#endif
