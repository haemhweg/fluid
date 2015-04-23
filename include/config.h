#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "real.h"

struct Config
{
  /**
   * Nur ein Vorschlag. Eigentlich muss man hier den Config irgendwie enkapsulieren.
   */
  struct geo
  {
    REAL xlength=10;
    REAL ylength=10;
    unsigned imax=50;
    unsigned jmax=50;
    REAL delx=0.2;
    REAL dely=0.2;
  };
  struct time
  {
    REAL t_end=2.;
    REAL delt=0.02;
    REAL tau=0.5;
    REAL del_vec=2.0;
  };
  struct solver
  {
    unsigned itmax=100;
    REAL eps=0.001;
    REAL omg=1.7;
    REAL alpha=0.5;
  };
  struct constants
  {
    REAL Re=10;
    REAL GX=0.0, GY=0.0;
    REAL UI=0.0, VI=0.0, PI=0.0;
  };

  geo _geo=geo();
  time _time=time();
  solver _pressure=solver();
  constants _constants=constants();
};

#endif
