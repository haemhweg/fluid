#include <iostream>
#include <algorithm>
#include <string>
#include <cmath>

#include "real.h"
#include "matrix.h"
#include "config.h"

#include "differences.h"

#include "solvers.h"

#include "velocity.h"

#include "specialBoundary.h"

REAL compDelt(Config::geo geoConfig, Config::time timeConfig, Config::constants constantsConfig, 
	      const Velocity& Velocity) 
{	
  REAL vals[3];

  REAL maxU = Velocity.getMaxU();
  REAL maxV = Velocity.getMaxV();

  if(maxU!=0 && maxV!=0) {
    vals[0] = constantsConfig.Re / 2 / (1 / geoConfig.delx / geoConfig.delx + 1 / geoConfig.dely / geoConfig.dely);
    vals[1] = geoConfig.delx / maxU;
    vals[2] = geoConfig.dely / maxV;  

    return timeConfig.tau * *std::min_element(vals, vals + 3);
  }
  else return timeConfig.tau * constantsConfig.Re / 2 / (  1 / geoConfig.delx / geoConfig.delx 
							 + 1 / geoConfig.dely / geoConfig.dely);
}

int main()
{  
  Config conf{"config"};

  Geometry geometry(conf._geo, geometry_STEP);

  Matrix Pressure{conf._geo.imax + 2, conf._geo.jmax + 2, conf._constants.PI};
  Velocity Velocity{conf._geo, conf._constants, conf._bc, geometry, bc_DRIVEN_CAVITY};

  REAL t=0;
  REAL delt=0;
  unsigned step=0;
  
  while(t<conf._time.t_end ){
    delt = compDelt(conf._geo, conf._time, conf._constants, Velocity);

    auto RHS = Velocity.getDivergenceIntermidiate(delt);

    auto it_res = SOR_Poisson(conf._geo, conf._solver, Pressure, RHS);

    std::cout << "Schritt " << step  << ": delt = " << delt 
  	      << ", Iterationen: " << it_res.first << ", Residuum: " << it_res.second << std::endl;

    Velocity.update(delt, Pressure);

    Velocity.writeVTK(step);
    Pressure.writeVTK("Pressure"+std::to_string(step)+".vtk", "Pressure", conf._geo.delx, conf._geo.dely);

    t += delt;    
    ++step;
  }
}
