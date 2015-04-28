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

void initUVP(Config::geo geoConfig, Config::constants constantsConfig, Matrix & U, Matrix & V, Matrix & P)
{
  U = Matrix(geoConfig.imax + 2, geoConfig.jmax + 2, constantsConfig.UI);
  V = Matrix(geoConfig.imax + 2, geoConfig.jmax + 2, constantsConfig.VI);
  P = Matrix(geoConfig.imax + 2, geoConfig.jmax + 2, constantsConfig.PI);
}

REAL compDelt(Config::geo geoConfig, Config::time timeConfig, Config::constants constantsConfig, 
	      const Matrix& U, const Matrix& V) 
{	
  REAL vals[3];

  if(U.getMax()!=0 && V.getMax()!=0) {
    vals[0] = constantsConfig.Re / 2 / (1 / geoConfig.delx / geoConfig.delx + 1 / geoConfig.dely / geoConfig.dely);
    vals[1] = geoConfig.delx / U.getMax();
    vals[2] = geoConfig.dely / V.getMax();  

    return timeConfig.tau * *std::min_element(vals, vals + 3);
  }
  else return timeConfig.tau * constantsConfig.Re / 2 / (  1 / geoConfig.delx / geoConfig.delx 
							 + 1 / geoConfig.dely / geoConfig.dely);
}

void outputVTK(const Config::geo geoConfig, const Matrix& U, const Matrix& V, const Matrix& P, const unsigned n)
{
  std::string filenameUV = "UV" + std::to_string(n) + ".vtk";
  std::string filenameP = "P" + std::to_string(n) + ".vtk";
  writeVectorFieldVTK(filenameUV, "Velocity", U, V, geoConfig.delx, geoConfig.dely);
  P.writeVTK(filenameP, "Pressure", geoConfig.delx, geoConfig.dely);
}


int main()
{
  Matrix U{};
  Matrix V{};
  Matrix P{};
  
  Config conf{"config"};

  REAL t=0;
  REAL delt=0;
  unsigned n=0;

  initUVP(conf._geo, conf._constants, U, V, P);
  
  while(t<conf._time.t_end ){
    delt = compDelt(conf._geo, conf._time, conf._constants, U, V);

    auto FG = compIntermediateVelocity(conf._geo, conf._constants, delt, U, V);

    auto RHS = RHS_Poisson(conf._geo, delt, FG.first, FG.second);

    auto it_res = SOR_Poisson(conf._geo, conf._solver, P, RHS);

    std::cout << "Schritt " << n  << ": delt = " << delt 
	      << ", Iterationen: " << it_res.first << ", Residuum: " << it_res.second << std::endl;
    
    compNewVelocity(conf._geo, delt, U, V, FG.first, FG.second, P);

    outputVTK(conf._geo, U, V, P, n++);

    t += delt;    
  }
}
