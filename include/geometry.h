#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <functional>

#include "config.h"

enum CELL { FLUID=0, BLOCK, B_N, B_O, B_S, B_W, B_NW, B_SW, B_SO, B_NO };

class Geometry
{
  std::vector<CELL> cells;
  unsigned imax, jmax;
  
  std::vector<std::pair<unsigned,unsigned>> cells_boundary;
  std::vector<std::pair<unsigned,unsigned>> cells_fluid;

 public:

  Geometry(const Config::geo& geoConfig, const std::function<std::vector<CELL>(const Config::geo&)> fill_geometry);

  CELL at(unsigned i, unsigned j) const { return cells[i*(jmax+2)+j]; }
  CELL& at(unsigned i, unsigned j) { return cells[i*(jmax+2)+j]; }

  const std::vector<std::pair<unsigned,unsigned>>& get_boundary() const { return cells_boundary; }
  const std::vector<std::pair<unsigned,unsigned>>& get_fluid() const { return cells_fluid; }

  void print();
};

#endif
