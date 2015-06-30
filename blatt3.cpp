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
  MPI_Init(0, 0);
  init_MPI_Grid(comm_grid);

  int rank = get_MPI_Comm_rank(comm_grid);

  Config conf{"config_DRIVEN_CAVITY", comm_grid};

  double F_0 = rank, G_0 = rank;
  Matrix F{conf._geo.imax + 2, conf._geo.jmax + 2, F_0};
  Matrix G{conf._geo.imax + 2, conf._geo.jmax + 2, G_0};

  // MPI_Barrier(comm_grid);
  // Matrix_exchange(comm_grid, Pressure);
  // MPI_Barrier(comm_grid);
  if(rank==0){
    system("rm *.vtk");
  }
  
  MPI_Barrier(comm_grid);
  MPI_VectorFieldVTK(comm_grid, "test.vtk", "test", F, G, conf._geo.delx, conf._geo.dely);
  int size = F.getCols()*F.getRows();
  if(rank==0){
    MPI_Status status;
    Matrix F1{conf._geo.imax + 2, conf._geo.jmax + 2, F_0};
    Matrix G1{conf._geo.imax + 2, conf._geo.jmax + 2, G_0};
    Matrix F2{conf._geo.imax + 2, conf._geo.jmax + 2, F_0};
    Matrix G2{conf._geo.imax + 2, conf._geo.jmax + 2, G_0};
    Matrix F3{conf._geo.imax + 2, conf._geo.jmax + 2, F_0};
    Matrix G3{conf._geo.imax + 2, conf._geo.jmax + 2, G_0};
    Matrix F_{8, 8, F_0};
    Matrix G_{8, 8, G_0};

    int rank_src = get_MPI_Cart_rank(comm_grid, std::array<int,2>{0,1}.data());
    MPI_Recv((void*)F1.begin(), size, MPI_DOUBLE, rank_src, rank_src, comm_grid, &status);
    MPI_Recv((void*)G1.begin(), size, MPI_DOUBLE, rank_src, size*rank_src, comm_grid, &status);

    rank_src = get_MPI_Cart_rank(comm_grid, std::array<int,2>{1,0}.data());
    MPI_Recv((void*)F2.begin(), size, MPI_DOUBLE, rank_src, rank_src, comm_grid, &status);
    MPI_Recv((void*)G2.begin(), size, MPI_DOUBLE, rank_src, size*rank_src, comm_grid, &status);

    rank_src = get_MPI_Cart_rank(comm_grid, std::array<int,2>{1,1}.data());
    MPI_Recv((void*)F3.begin(), size, MPI_DOUBLE, rank_src, rank_src, comm_grid, &status);
    MPI_Recv((void*)G3.begin(), size, MPI_DOUBLE, rank_src, size*rank_src, comm_grid, &status);

    for(unsigned i=0; i<F_.getRows(); ++i){
      for(unsigned j=0; j<F_.getCols(); ++j){
	if(i<4 && j<4){
	  F_.at(i,j) = F.at(i,j);
	  G_.at(i,j) = G.at(i,j);
	}else if(i>=4 && j<4){
	  F_.at(i,j) = F2.at(i+1,j);
	  G_.at(i,j) = G2.at(i+1,j);
	}else if(i<4 && j>=4){
	  F_.at(i,j) = F1.at(i,j+1);
	  G_.at(i,j) = G1.at(i,j+1);
	}else if(i>=4 && j>=4){
	  F_.at(i,j) = F3.at(i+1,j+1);
	  G_.at(i,j) = G3.at(i+1,j+1);
	}
      }
    }
    writeVectorFieldVTK("test2.vtk", "test", F_, G_, conf._geo.delx, conf._geo.dely);
  }
  else{
    MPI_Send((void*)F.begin(), size, MPI_DOUBLE, 0, rank, comm_grid);
    MPI_Send((void*)G.begin(), size, MPI_DOUBLE, 0, size*rank, comm_grid);
  }
  MPI_Barrier(comm_grid);
  

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
