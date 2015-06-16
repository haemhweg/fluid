#include <iostream>
#include <algorithm>
#include <string>
#include <cmath>

#include "mpi.h"
#include "real.h"
#include "matrix.h"
#include "config.h"
#include "differences.h"
#include "solvers.h"
#include "velocity.h"
#include "specialBoundary.h"
#include "parallel.h"

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

int main(int argc, char* argv[])
{  
  MPI_Comm comm_grid;
  Config conf{"config_DRIVEN_CAVITY"};

  init_MPI_Grid(conf._geo, comm_grid);

  MPI_Finalize();
  
  // std::string cfg_PROBLEM;
  // special_boundary_fct bc_sp_PROBLEM;
  // init_geometry_fct initGeometry_PROBLEM;

  // auto toUpper = [] (std::string& str) -> void { std::transform(str.begin(), str.end(),str.begin(), ::toupper); };
  
  // if(argc<2){
  //   cfg_PROBLEM = "config_DIRVEN_CAVATY";
  //   bc_sp_PROBLEM = bc_DRIVEN_CAVITY;
  //   initGeometry_PROBLEM = geometry_DRIVEN_CAVITY;
  // }
  // else{
  //   std::string problem(argv[1]);
  //   toUpper(problem);
    
  //   if(problem==std::string("DRIVEN_CAVITY")){
  //     cfg_PROBLEM = "config_DIRVEN_CAVATY";
  //     bc_sp_PROBLEM = bc_DRIVEN_CAVITY;
  //     initGeometry_PROBLEM = geometry_DRIVEN_CAVITY;      
  //   }else if(problem==std::string("STEP")){
  //     cfg_PROBLEM = "config_STEP";
  //     bc_sp_PROBLEM = bc_STEP;
  //     initGeometry_PROBLEM = geometry_STEP; 	
  //   }else if(problem==std::string("KARMAN")){
  //     cfg_PROBLEM = "config_KARMAN";
  //     bc_sp_PROBLEM = bc_KARMAN;
  //     initGeometry_PROBLEM = geometry_KARMAN;   	  
  //   }else{
  //     std::cout << "Case >" << argv[1] << "< not implemented." << std::endl;
  //     return 1;
  //   }
  // }
      
  // Config conf{cfg_PROBLEM};
  
  // Geometry geometry(conf._geo, initGeometry_PROBLEM);
  
  // Matrix Pressure{conf._geo.imax + 2, conf._geo.jmax + 2, conf._constants.PI};
  // Matrix Div_velocity{conf._geo.imax + 1, conf._geo.jmax + 1, 0};
  // Velocity Velocity{conf._geo, conf._constants, conf._bc, geometry, bc_sp_PROBLEM};
  
  // REAL t=0;
  // REAL delt=0;
  // REAL next_output=conf._time.del_vec;
  // unsigned step=1;

  // system("rm *.vtk");

  // Velocity.writeVTK(0);
  // Pressure.writeVTK("Pressure0.vtk", "Pressure", conf._geo.delx, conf._geo.dely);
    
  // while(t<conf._time.t_end ){
  //   delt = compDelt(conf._geo, conf._time, conf._constants, Velocity);

  //   Velocity.setDivergenceIntermidiate(delt, Div_velocity);

  //   auto it_res = SOR_Poisson(conf._geo, conf._solver, geometry, Pressure, Div_velocity);

  //   Velocity.update(delt, Pressure);

  //   if(t>next_output){      
  //     std::cout << "Ausgabe " << step  << ": delt = " << delt 
  // 		<< ", Iterationen: " << it_res.first << ", Residuum: " << it_res.second << std::endl;

  //     Velocity.writeVTK(step);
  //     Pressure.writeVTK("Pressure"+std::to_string(step)+".vtk", "Pressure", conf._geo.delx, conf._geo.dely);
      
  //     next_output += conf._time.del_vec;   
  //     ++step;
  //   }

  //   t += delt; 
  // }

  // Velocity.writeVTK(step);
  // Pressure.writeVTK("Pressure"+std::to_string(step)+".vtk", "Pressure", conf._geo.delx, conf._geo.dely);
}
