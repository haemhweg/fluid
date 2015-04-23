#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <string>

#include "real.h"

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
    REAL delt;
    REAL tau;
    REAL del_vec;
  };
  struct solver
  {
    int itmax;
    REAL eps;
    REAL omega;
    REAL alpha;
  };
  struct constants
  {
    REAL Re;
    REAL GX, GY;
    REAL UI, VI, PI;
  };

  geo _geo=geo();
  time _time=time();
  solver _solver=solver();
  constants _constants=constants();

  Config(const std::string& filename);
};

#endif
