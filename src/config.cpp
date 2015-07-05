#include <fstream>
#include <string>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <cstdlib>

#include "config.h"

Config::Config(const std::string& filename)
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
    _geo.delx = REAL(_geo.xlength)/_geo.imax;
    _geo.dely = REAL(_geo.ylength)/_geo.jmax;
  }
  else{
    std::cout << "Wrong input!" << std::endl;
    std::exit(EXIT_FAILURE);

  }  
}


void Config::print(const std::string& example)
{
  std::cout << "Solving the Navier-Stokes equation for the " << example << " example.\n";

  std::cout << "\nGeometry parameters:"<< "\n"
	    << "\txlength = " << _geo.xlength << "\n"
	    << "\tylength = " << _geo.ylength << "\n"
	    << "\timax = " << _geo.imax << "\n"
	    << "\tjmax = " << _geo.jmax << "\n";

  std::cout << "\nTimestep parameters:"<< "\n"
	    << "\tt_end = " << _time.t_end << "\n"
	    << "\ttau = " << _time.tau  << "\n"
	    << "\tdel_vec = " << _time.del_vec << "\n";

  std::cout << "\nSOR parameters:" << "\n"
	    << "\teps = " << _solver.eps << "\n"
	    << "\titmax = " << _solver.itmax << "\n"
	    << "\tomega = " << _solver.omega << "\n";

  std::cout << "\nConstants:" << "\n"
	    << "\talpha = " << _constants.alpha << "\n"
	    << "\tG = (" << _constants.GX << ", " << _constants.GY << ")" << "\n"
	    << "\tRe = " << _constants.Re << "\n"
	    << "\tInit velocity = (" << _constants.UI << ", " << _constants.VI << ")" << "\n"
	    << "\tInit pressure = " << _constants.PI << "\n";

  std::cout << "\nBoundary Conditions (without special boundary):" << "\n";
  std::cout << "\tTop = ";
  switch(_bc.wt){
  case NO_SLIP:
    std::cout << "NO_SLIP" << "\n";
    break;
  case FREE_SLIP:
    std::cout << "FREE_SLIP" << "\n";
    break;
  case OUTFLOW:
    std::cout << "OUTFLOW" << "\n";
    break;
  }
  std::cout << "\tRight = ";
  switch(_bc.wr){
  case NO_SLIP:
    std::cout << "NO_SLIP" << "\n";
    break;
  case FREE_SLIP:
    std::cout << "FREE_SLIP" << "\n";
    break;
  case OUTFLOW:
    std::cout << "OUTFLOW" << "\n";
    break;
  }
  std::cout << "\tBottom = ";
  switch(_bc.wb){
  case NO_SLIP:
    std::cout << "NO_SLIP" << "\n";
    break;
  case FREE_SLIP:
    std::cout << "FREE_SLIP" << "\n";
    break;
  case OUTFLOW:
    std::cout << "OUTFLOW" << "\n";
    break;
  }
  std::cout << "\tLeft = ";
  switch(_bc.wl){
  case NO_SLIP:
    std::cout << "NO_SLIP" << "\n";
    break;
  case FREE_SLIP:
    std::cout << "FREE_SLIP" << "\n";
    break;
  case OUTFLOW:
    std::cout << "OUTFLOW" << "\n";
    break;
  }
  std::cout << std::endl;    
}
