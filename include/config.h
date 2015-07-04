#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <string>

#include "real.h"

enum  BCType { NO_SLIP=1, FREE_SLIP=2, OUTFLOW=3 };

struct Config
{

  struct geo
  {
    REAL xlength;
    REAL ylength;
    unsigned imax;
    unsigned jmax;
    REAL delx;
    REAL dely;
  };
  struct time
  {
    REAL t_end;
    REAL tau;
    REAL del_vec;
  };
  struct solver
  {
    unsigned itmax;
    REAL eps;
    REAL omega;
  };
  struct constants
  {
    REAL Re;
    REAL alpha;
    REAL GX, GY;
    REAL UI, VI, PI;
  };
  struct boundaryCondition
  {
    BCType wl, wr, wt, wb;
    std::string problem="default";
  };


  geo _geo=geo();
  time _time=time();
  solver _solver=solver();
  constants _constants=constants();
  boundaryCondition _bc=boundaryCondition();

  Config(const std::string& filename);

  void print(const std::string& problem);
};

#endif
