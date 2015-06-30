#ifndef PROBLEMS
#define PROBLEMS

#include <functional>
#include <vector>

#include "matrix.h"
#include "config.h"

enum CELL { FLUID=0, BLOCK, B_N, B_O, B_S, B_W, B_NW, B_SW, B_SO, B_NO };

using special_boundary_fct = std::function<void(const unsigned, const unsigned, Matrix&, Matrix&)>;

using init_geometry_fct = std::function<std::vector<CELL>(const MPI_Comm&, const Config::geo&)>;

void bc_DRIVEN_CAVITY(const unsigned imax, const unsigned jmax, Matrix& U, Matrix&); 

void bc_STEP(const unsigned, const unsigned jmax, Matrix& U, Matrix&);

void bc_KARMAN(const unsigned, const unsigned jmax, Matrix& U, Matrix&);

std::vector<CELL> geometry_DRIVEN_CAVITY(const MPI_Comm& comm_grid, const Config::geo& geoConfig);

std::vector<CELL> geometry_STEP(const MPI_Comm& comm_grid, const Config::geo& geoConfig);

std::vector<CELL> geometry_KARMAN(const MPI_Comm& comm_grid, const Config::geo& geoConfig);

#endif
