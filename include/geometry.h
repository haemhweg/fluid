#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <functional>

#include "config.h"
#include "specialBoundary.h"

class Geometry
{
  std::vector<CELL> cells;
  unsigned imax, jmax;
  
  std::vector<std::pair<unsigned,unsigned>> cells_boundary;
  std::vector<std::pair<unsigned,unsigned>> cells_fluid;

 public:

  Geometry(const Config::geo& geoConfig, const init_geometry_fct& geoInit);

  CELL at(unsigned i, unsigned j) const { return cells[i*(jmax+2)+j]; }
  CELL& at(unsigned i, unsigned j) { return cells[i*(jmax+2)+j]; }

  const std::vector<std::pair<unsigned,unsigned>>& get_boundary() const { return cells_boundary; }
  const std::vector<std::pair<unsigned,unsigned>>& get_fluid() const { return cells_fluid; }

  void print();
};

#endif
