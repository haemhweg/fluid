#include <algorithm>
#include <iostream>

#include "parallel.h"
#include "mpi.h"


int get_MPI_Comm_size(MPI_Comm comm)
{
  int size;

  MPI_Comm_size(comm, &size);

  return size;
}
int get_MPI_Comm_rank(MPI_Comm comm)
{
  int rank;

  MPI_Comm_rank(comm, &rank);

  return rank;
}


int init_MPI_Grid(const Config::geo& geoConfig, MPI_Comm comm_grid)
{
  int nprocs = get_MPI_Comm_size(MPI_COMM_WORLD);
  int rank = get_MPI_Comm_rank(MPI_COMM_WORLD);
  int dims[2] = {0,0};
  int periods[2] = {0,0};

  unsigned imax = geoConfig.imax, jmax = geoConfig.jmax;

  MPI_Dims_create(nprocs, 2, dims);

  if((imax>jmax && dims[0]<dims[1]) || (imax<jmax && dims[0]>dims[1])){
    std::swap(dims[0], dims[1]);
  }

  if(rank==0){
    std::cout << "Proposal from MPI for " << nprocs << " Nodes and 2 Dims : " 
	      << dims[0] << " " << dims[1] << std::endl;
  }

  return MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm_grid);
}
