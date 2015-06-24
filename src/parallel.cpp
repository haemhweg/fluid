#include <algorithm>
#include <iostream>
#include <array>
#include "mpi.h"

#include "parallel.h"
#include "matrix.h"
#include "config.h"


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
int get_MPI_Cart_rank(MPI_Comm comm, const int* coords)
{
  int rank;

  MPI_Cart_rank(comm, coords, &rank);

  return rank;
}
std::array<int,2> get_MPI_Dims_create(MPI_Comm comm, int ndims)
{
  int nprocs = get_MPI_Comm_size(comm);
  int dims[2] ={0,0};

  MPI_Dims_create(nprocs, ndims, dims);

  return std::array<int,2>{dims[0],dims[1]};
}
std::array<int,2> get_MPI_Cart_coords(MPI_Comm comm_grid, int ndims)
{
  int rank = get_MPI_Comm_rank(comm_grid);
  int coords[2] ={0,0};

  MPI_Cart_coords(comm_grid, rank, ndims, coords);

  return std::array<int,2>{coords[0],coords[1]};
}


int init_MPI_Grid(MPI_Comm& comm_grid)
{
  int rank = get_MPI_Comm_rank(MPI_COMM_WORLD);
  int nprocs = get_MPI_Comm_size(MPI_COMM_WORLD);
  std::array<int,2> dims = get_MPI_Dims_create(MPI_COMM_WORLD, 2);
  int periods[2] = {0,0};

  if(rank==0){
    std::cout << "Proposal from MPI for " << nprocs << " Nodes and 2 Dims : " 
	      << dims[0] << " " << dims[1] << std::endl;
  }

  return MPI_Cart_create(MPI_COMM_WORLD, 2, dims.data(), periods, 0, &comm_grid);
}


void Matrix_exchange(const MPI_Comm& comm_grid, Matrix& M)
{
  std::array<int,2> dims = get_MPI_Dims_create(MPI_COMM_WORLD,2);
  std::array<int,2> coords_src, coords_dest, mycoords = get_MPI_Cart_coords(comm_grid,2);
  std::vector<double> buffer_send(std::max(M.getRows(),M.getCols()), 0), buffer_recv(std::max(M.getRows(),M.getCols()), 0);
  unsigned n = M.getCols(), m = M.getRows();
  int rank = get_MPI_Comm_rank(comm_grid);

  MPI_Request request;
  MPI_Status status;

  // send right
  if(mycoords[0]<dims[0]-1){
    if(rank==0)
      std::cout << "send RIGHT" << std::endl;
    coords_dest[0] = mycoords[0]+1; coords_dest[1] = mycoords[1];
    for(unsigned i=0; i<m; ++i){
      buffer_send[i] = M.at(n-2,i);
    }
    MPI_Send((void*)buffer_send.data(), m, MPI_DOUBLE, get_MPI_Cart_rank(comm_grid, coords_dest.data()), 0, comm_grid);
  }
  // recive left
  if(mycoords[0]>0){
    if(rank==0)
      std::cout << "recive LEFT" << std::endl;
    coords_src[0] = mycoords[0]-1; coords_src[1] = mycoords[1];
    MPI_Recv((void*)buffer_recv.data(), m, MPI_DOUBLE, get_MPI_Cart_rank(comm_grid, coords_src.data()), 0, comm_grid, &status);
    for(unsigned i=0; i<m; ++i){
      M.at(0,i) = buffer_recv[i];
    }    
  }

  // send left
  if(mycoords[0]>0){
    if(rank==0)
      std::cout << "send LEFT" << std::endl;
    coords_dest[0] = mycoords[0]-1; coords_dest[1] = mycoords[1];
    for(unsigned i=0; i<m; ++i){
      buffer_send[i] = M.at(1,i);
    }
    MPI_Send((void*)buffer_send.data(), m, MPI_DOUBLE, get_MPI_Cart_rank(comm_grid, coords_dest.data()), 0, comm_grid);
  }
  // recive right
  if(mycoords[0]<dims[0]-1){
    if(rank==0)
      std::cout << "recive RIGHT" << std::endl;
    coords_src[0] = mycoords[0]+1; coords_src[1] = mycoords[1];
    MPI_Recv((void*)buffer_recv.data(), m, MPI_DOUBLE, get_MPI_Cart_rank(comm_grid, coords_src.data()), 0, comm_grid, &status);
    for(unsigned i=0; i<m; ++i){
      M.at(n-1,i) = buffer_recv[i];
    }    
  }



  // send up
  if(mycoords[1]<dims[1]-1){
    if(rank==0)
      std::cout << "send UP" << std::endl;
    coords_dest[0] = mycoords[0]; coords_dest[1] = mycoords[1]+1;
    for(unsigned i=0; i<n; ++i){
      buffer_send[i] = M.at(i,m-2);
    }
    MPI_Send((void*)buffer_send.data(), n, MPI_DOUBLE, get_MPI_Cart_rank(comm_grid, coords_dest.data()), 0, comm_grid);
  }
  // recive down
  if(mycoords[1]>0){
    if(rank==0)
      std::cout << "recive DOWN" << std::endl;
    coords_src[0] = mycoords[0]; coords_src[1] = mycoords[1]-1;
    MPI_Recv((void*)buffer_recv.data(), n, MPI_DOUBLE, get_MPI_Cart_rank(comm_grid, coords_src.data()), 0, comm_grid, &status);
    for(unsigned i=0; i<n; ++i){
      M.at(i,0) = buffer_recv[i];
    }    
  }

  // send down
  if(mycoords[1]>0){
    if(rank==0)
      std::cout << "send DOWN" << std::endl;
    coords_dest[0] = mycoords[0]; coords_dest[1] = mycoords[1]-1;
    for(unsigned i=0; i<n; ++i){
      buffer_send[i] = M.at(i,0);
    }
    MPI_Send((void*)buffer_send.data(), n, MPI_DOUBLE, get_MPI_Cart_rank(comm_grid, coords_dest.data()), 0, comm_grid);
  }
  // recive up
  if(mycoords[1]<dims[1]-1){
    if(rank==0)
      std::cout << "recive UP" << std::endl;
    coords_src[0] = mycoords[0]; coords_src[1] = mycoords[1]+1;
    MPI_Recv((void*)buffer_recv.data(), n, MPI_DOUBLE, get_MPI_Cart_rank(comm_grid, coords_src.data()), 0, comm_grid, &status);
    for(unsigned i=0; i<n; ++i){
      M.at(i,m-1) = buffer_recv[i];
    }    
  }
}
