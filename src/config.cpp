#include <fstream>
#include <string>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <cstdlib>
#include "mpi.h"

#include "config.h"
#include "parallel.h"

Config::Config(const std::string& filename, MPI_Comm comm_grid)
{
  std::ifstream file(filename);

  std::string input;

  std::vector<bool> is_set(21, false);


  while(file.good()) {
    file >> input;

    // convert input to lowercase
    std::transform(std::begin(input), std::end(input), std::begin(input), ::tolower);

    // geometry
    if(input == "xlength") {

      if(!is_set[0]) {
	is_set[0] = true;
	file >> _geo.xlength;
      }      
    } else if(input == "ylength") {
      
      if(!is_set[1]) {
	is_set[1] = true;
	file >> _geo.ylength;
      }
    }else if(input == "imax") {
      
      if(!is_set[2]) {
	is_set[2] = true;
	file >> _geo.imax;
      }
    }else if(input == "jmax") {
      
      if(!is_set[3]) {
	is_set[3] = true;
	file >> _geo.jmax;
      }
    }
    // time
    else if(input == "t_end") {
      
      if(!is_set[4]) {
	is_set[4] = true;
	file >> _time.t_end;
      }
    }else if(input == "tau") {
      
      if(!is_set[5]) {
	is_set[5] = true;
	file >> _time.tau;
      }
    }else if(input == "del_vec") {
      
      if(!is_set[6]) {
	is_set[6] = true;
	file >> _time.del_vec;
      }
    }
    // solver
    else if(input == "itmax") {
      
      if(!is_set[7]) {
	is_set[7] = true;
	file >> _solver.itmax;
      }
    }else if(input == "eps") {
      
      if(!is_set[8]) {
	is_set[8] = true;
	file >> _solver.eps;
      }
    }else if(input == "omega") {
      
      if(!is_set[9]) {
	is_set[9] = true;
	file >> _solver.omega;
      }
    }    
    // constants
    else if(input == "alpha") {
      
      if(!is_set[10]) {
	is_set[10] = true;
	file >> _constants.alpha;
      }
    }else if(input == "re") {
      
      if(!is_set[11]) {
	is_set[11] = true;
	file >> _constants.Re;
      }
    }else if(input == "gx") {
      
      if(!is_set[12]) {
	is_set[12] = true;
	file >> _constants.GX;
      }
    }else if(input == "gy") {
      
      if(!is_set[13]) {
	is_set[13] = true;
	file >> _constants.GY;
      }
    }else if(input == "ui") {
      
      if(!is_set[14]) {
	is_set[14] = true;
	file >> _constants.UI;
      }
    }else if(input == "vi") {
      
      if(!is_set[15]) {
	is_set[15] = true;
	file >> _constants.VI;
      }
    }else if(input == "pi") {
      
      if(!is_set[16]) {
	is_set[16] = true;
	file >> _constants.PI;
      }
    }else if(input == "wl") {
      
      if(!is_set[17]) {
	int bc_val;
	is_set[17] = true;
	file >> bc_val;
	_bc.wl = static_cast<BCType>(bc_val);
      }
    }else if(input == "wr") {
      
      if(!is_set[18]) {
	int bc_val;
	is_set[18] = true;
	file >> bc_val;
	_bc.wr = static_cast<BCType>(bc_val);
      }
    }else if(input == "wt") {
      
      if(!is_set[19]) {
	int bc_val;
	is_set[19] = true;
	file >> bc_val;
	_bc.wt = static_cast<BCType>(bc_val);
      }
    }else if(input == "wb") {
      
      if(!is_set[20]) {
	int bc_val;
	is_set[20] = true;
	file >> bc_val;
	_bc.wb = static_cast<BCType>(bc_val);
      }
    }
  }

  if(std::equal(std::begin(is_set), std::end(is_set), std::begin(std::vector<bool>(21,true)))) {
    std::array<int,2> coords = get_MPI_Cart_coords(comm_grid, 2);      
    std::array<int,2> dims = get_MPI_Dims_create(MPI_COMM_WORLD, 2);   

    unsigned imax = _geo.imax, jmax = _geo.jmax;

    if(coords[0]==dims[0]-1 && dims[0]>1){
      _geo.imax -= (dims[0]-1)*(_geo.imax/dims[0]);
      _geo.xlength -= (dims[0]-1)*(_geo.xlength/dims[0]);
    }else{
      _geo.imax /= dims[0];
      _geo.xlength /= dims[0];
    }
    if(coords[1]==dims[1]-1 && dims[1]>1){
      _geo.jmax -= (dims[1]-1)*(_geo.jmax/dims[1]);
      _geo.ylength -= (dims[1]-1)*(_geo.ylength/dims[1]);
    }else{
      _geo.jmax /= dims[1];
      _geo.ylength /= dims[1];
    }
      
    _geo.delx = REAL(_geo.xlength)/_geo.imax;
    _geo.dely = REAL(_geo.ylength)/_geo.jmax;
  }
  else{
    std::cout << "Wrong input!" << std::endl;
    std::exit(EXIT_FAILURE);

  }  
}
