#include <algorithm>
#include <iostream>
#include <vector>
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
std::vector<int> get_MPI_Dims_create(MPI_Comm comm, int ndims)
{
  int nprocs = get_MPI_Comm_size(comm);
  int dims[2] ={0,0};

  MPI_Dims_create(nprocs, ndims, dims);

  return std::vector<int>{dims[0],dims[1]};
}
std::vector<int> get_MPI_Cart_coords(MPI_Comm comm_grid, int ndims)
{
  int rank = get_MPI_Comm_rank(comm_grid);
  int coords[2] ={0,0};

  MPI_Cart_coords(comm_grid, rank, ndims, coords);

  return std::vector<int>{coords[0],coords[1]};
}


int init_MPI_Grid(MPI_Comm& comm_grid)
{
  int rank = get_MPI_Comm_rank(MPI_COMM_WORLD);
  int nprocs = get_MPI_Comm_size(MPI_COMM_WORLD);
  std::vector<int> dims = get_MPI_Dims_create(MPI_COMM_WORLD, 2);
  int periods[2] = {0,0};

  if(rank==0){
    std::cout << "Proposal from MPI for " << nprocs << " Nodes and 2 Dims : " 
	      << dims[0] << " " << dims[1] << std::endl;
  }

  return MPI_Cart_create(MPI_COMM_WORLD, 2, dims.data(), periods, 0, &comm_grid);
}


void Matrix_exchange(const MPI_Comm& comm_grid, Matrix& M)
{
  std::vector<int> dims = get_MPI_Dims_create(MPI_COMM_WORLD,2);
  std::vector<int> coords_src(2), coords_dest(2), mycoords = get_MPI_Cart_coords(comm_grid,2);
  std::vector<double> buffer_send(std::max(M.getRows(),M.getCols()), 0), buffer_recv(std::max(M.getRows(),M.getCols()), 0);
  unsigned n = M.getCols(), m = M.getRows();
  int rank = get_MPI_Comm_rank(comm_grid), tagLR=1, tagRL=2, tagUD=3, tagDU=4;

  MPI_Request request;
  MPI_Status status;

  // send right
  if(mycoords[0]<dims[0]-1){
    // if(rank==0)
    //   std::cout << "send RIGHT" << std::endl;
    coords_dest[0] = mycoords[0]+1; coords_dest[1] = mycoords[1];
    for(unsigned i=0; i<m; ++i){
      buffer_send[i] = M.at(n-2,i);
    }
    MPI_Send((void*)buffer_send.data(), buffer_send.size(), MPI_DOUBLE, get_MPI_Cart_rank(comm_grid, coords_dest.data()), tagRL, comm_grid);
  }
  // recive left
  if(mycoords[0]>0){
    // if(rank==0)
    //   std::cout << "recive LEFT" << std::endl;
    coords_src[0] = mycoords[0]-1; coords_src[1] = mycoords[1];
    MPI_Recv((void*)buffer_recv.data(), buffer_recv.size(), MPI_DOUBLE, get_MPI_Cart_rank(comm_grid, coords_src.data()), tagRL, comm_grid, &status);
    for(unsigned i=0; i<m; ++i){
      M.at(0,i) = buffer_recv[i];
    }    
  }

  // send left
  if(mycoords[0]>0){
    // if(rank==0)
    //   std::cout << "send LEFT" << std::endl;
    coords_dest[0] = mycoords[0]-1; coords_dest[1] = mycoords[1];
    for(unsigned i=0; i<m; ++i){
      buffer_send[i] = M.at(1,i);
    }
    MPI_Send((void*)buffer_send.data(), m, MPI_DOUBLE, get_MPI_Cart_rank(comm_grid, coords_dest.data()), tagLR, comm_grid);
  }
  // recive right
  if(mycoords[0]<dims[0]-1){
    // if(rank==0)
    //   std::cout << "recive RIGHT" << std::endl;
    coords_src[0] = mycoords[0]+1; coords_src[1] = mycoords[1];
    MPI_Recv((void*)buffer_recv.data(), m, MPI_DOUBLE, get_MPI_Cart_rank(comm_grid, coords_src.data()), tagLR, comm_grid, &status);
    for(unsigned i=0; i<m; ++i){
      M.at(n-1,i) = buffer_recv[i];
    }    
  }



  // send up
  if(mycoords[1]<dims[1]-1){
    coords_dest[0] = mycoords[0]; coords_dest[1] = mycoords[1]+1;
    for(unsigned i=0; i<n; ++i){
      buffer_send[i] = M.at(i,m-2);
    }
    MPI_Send((void*)buffer_send.data(), n, MPI_DOUBLE, get_MPI_Cart_rank(comm_grid, coords_dest.data()), tagUD, comm_grid);
  }
  // recive down
  if(mycoords[1]>0){
    coords_src[0] = mycoords[0]; coords_src[1] = mycoords[1]-1;
    MPI_Recv((void*)buffer_recv.data(), n, MPI_DOUBLE, get_MPI_Cart_rank(comm_grid, coords_src.data()), tagUD, comm_grid, &status);
    for(unsigned i=0; i<n; ++i){
      M.at(i,0) = buffer_recv[i];
    }    
  }

  // send down
  if(mycoords[1]>0){
    coords_dest[0] = mycoords[0]; coords_dest[1] = mycoords[1]-1;
    for(unsigned i=0; i<n; ++i){
      buffer_send[i] = M.at(i,1);
    }
    MPI_Send((void*)buffer_send.data(), n, MPI_DOUBLE, get_MPI_Cart_rank(comm_grid, coords_dest.data()), tagDU, comm_grid);
  }
  // recive up
  if(mycoords[1]<dims[1]-1){
    coords_src[0] = mycoords[0]; coords_src[1] = mycoords[1]+1;
    MPI_Recv((void*)buffer_recv.data(), n, MPI_DOUBLE, get_MPI_Cart_rank(comm_grid, coords_src.data()), tagDU, comm_grid, &status);
    for(unsigned i=0; i<n; ++i){
      M.at(i,m-1) = buffer_recv[i];
    }    
  }
}

void MPI_Matrix_print(const MPI_Comm& comm_grid, const Matrix& M)
{
  int rank = get_MPI_Comm_rank(comm);
  

  if(rank==0){
    int nprocs = get_MPI_Comm_size(comm);
    Matrix M_tmp(M.getRows(), M.getCols());
    for(int i=1; i<nprocs; ++i){
      MPI_Recv((void*)M_tmp.begin(), M.getRows()*M.getCols(), MPI_DOUBLE, 0, i, comm);
      
    }    
  }else{
    MPI_Send((void*)M.begin(), M.getRows()*M.getCols(), MPI_DOUBLE, 0, rank, comm);
  }
}
