#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <string>

#include "real.h"

enum  BCType { NO_SLIP=1, FREE_SLIP=2, OUTFLOW=3 };
enum  TracingType { PATHLINES=1, STREAKLINES=2 };

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
  struct tracing
  {
    unsigned int N;
    REAL x1, x2, y1, y2;
    TracingType tr;
    REAL delt_write, delt_inject;
  };


  geo _geo=geo();
  time _time=time();
  solver _solver=solver();
  constants _constants=constants();
  boundaryCondition _bc=boundaryCondition();
  tracing _tracing=tracing();

  Config(const std::string& filename);
};

#endif
