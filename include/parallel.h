#ifndef PARALLEL_H
#define PARALLEL_H

#include "mpi.h"
#include "config.h"

int get_MPI_Comm_size(MPI_Comm comm);
int get_MPI_Comm_rank(MPI_Comm comm);


int init_MPI_Grid(const Config::geo& geoConfig, MPI_Comm comm_grid);

#endif
