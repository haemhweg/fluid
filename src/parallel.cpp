#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <iterator>
#include <cassert>
#include "mpi.h"

#include "parallel.h"
#include "matrix.h"
#include "config.h"


int get_MPI_Comm_size(const MPI_Comm& comm)
{
  int size;

  MPI_Comm_size(comm, &size);

  return size;
}
int get_MPI_Comm_rank(const MPI_Comm& comm)
{
  int rank;

  MPI_Comm_rank(comm, &rank);

  return rank;
}
int get_MPI_Cart_rank(const MPI_Comm& comm, const int* coords)
{
  int rank;

  MPI_Cart_rank(comm, coords, &rank);

  return rank;
}
std::array<int,2> get_MPI_Dims_create(const MPI_Comm& comm, int ndims)
{
  int nprocs = get_MPI_Comm_size(comm);
  int dims[2] ={0,0};

  MPI_Dims_create(nprocs, ndims, dims);

  return std::array<int,2>{dims[0],dims[1]};
}
std::array<int,2> get_MPI_Cart_coords(const MPI_Comm& comm_grid, int ndims)
{
  int rank = get_MPI_Comm_rank(comm_grid);
  int coords[2] ={0,0};

  MPI_Cart_coords(comm_grid, rank, ndims, coords);

  return std::array<int,2>{coords[0],coords[1]};
}
std::array<int,2> get_MPI_Cart_dims(const MPI_Comm& comm_grid){
  int dims[2];

  MPI_Cartdim_get(comm_grid, dims);

  return std::array<int,2>{dims[0],dims[1]};
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
  int tagLR=1, tagRL=2, tagUD=3, tagDU=4;

  MPI_Status status;

  // send right
  if(mycoords[0]<dims[0]-1){
    // if(rank==0)
    //   std::cout << "send RIGHT" << std::endl;
    coords_dest[0] = mycoords[0]+1; coords_dest[1] = mycoords[1];
    for(unsigned i=0; i<m; ++i){
      buffer_send[i] = M.at(n-2,i);
    }
    MPI_Send((void*)buffer_send.data(), m, MPI_DOUBLE, get_MPI_Cart_rank(comm_grid, coords_dest.data()), tagRL, comm_grid);
  }
  // recive left
  if(mycoords[0]>0){
    // if(rank==0)
    //   std::cout << "recive LEFT" << std::endl;
    coords_src[0] = mycoords[0]-1; coords_src[1] = mycoords[1];
    MPI_Recv((void*)buffer_recv.data(), m, MPI_DOUBLE, get_MPI_Cart_rank(comm_grid, coords_src.data()), tagRL, comm_grid, &status);
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


void MPI_VectorFieldVTK(const MPI_Comm& comm_grid, const std::string& filename, const std::string& descr, 
			const Matrix& U, const Matrix& V, const double dx, const double dy)
{
  assert(U.getRows() == V.getRows() && U.getCols() == V.getCols());
  const int rank = get_MPI_Comm_rank(comm_grid), nprocs = get_MPI_Comm_size(comm_grid);
  const std::array<int,2> coords = get_MPI_Cart_coords(comm_grid, 2);
  const std::array<int,2> dims = get_MPI_Dims_create(MPI_COMM_WORLD, 2);
  unsigned M_loc=U.getRows(), N_loc=U.getCols(), M=0, N=0, M_max;

  MPI_Request request;
  MPI_Status status, status_dims;
  

  if(rank==0){  
    unsigned M_tmp, N_tmp;
    for(int i=0; i<dims[0]; ++i){
      if(coords[1]==0 && coords[0]==i){
	M += M_loc-1;
      }else{
	MPI_Recv((void*)&M_tmp, 1, MPI_UNSIGNED, get_MPI_Cart_rank(comm_grid, std::array<int,2>{i,0}.data()), 1, comm_grid,
		 &status);
	M += M_tmp-1;
      }
    }
    for(int j=0; j<dims[1]; ++j){
      if(coords[0]==0 && coords[1]==j){
    	N += N_loc-1;
      }else{
    	MPI_Recv((void*)&N_tmp, 1, MPI_UNSIGNED, get_MPI_Cart_rank(comm_grid, std::array<int,2>{0,j}.data()), 2, comm_grid,
    		 &status);
	N += N_tmp-1;
      }
    }
  }else{
    if(coords[1]==0){
      MPI_Send((void*)&M_loc, 1, MPI_UNSIGNED, 0, 1, comm_grid);
    }
    if(coords[0]==0){
      MPI_Send((void*)&N_loc, 1, MPI_UNSIGNED, 0, 2, comm_grid);
    }
  }
  MPI_Barrier(comm_grid);
  MPI_Reduce((void*)&M_loc, (void*)&M_max, 1, MPI_UNSIGNED, MPI_MAX, 0, comm_grid);

  if(rank==0){
    std::ofstream fs(filename);

    fs << std::fixed<< std::setprecision(6);

    fs << "# vtk DataFile Version 3.0\n"
       << "Vector Field\n"
       << "ASCII\n"
       << "DATASET RECTILINEAR_GRID\n"
       << "DIMENSIONS " << M << " " << N << " 1\n"
       << "X_COORDINATES " << M << " double\n";

    for(unsigned i=0; i<M; ++i) {
      fs << dx*i << " ";
    }

    fs << "\nY_COORDINATES " << N << " double\n";
  
    for(unsigned j=0; j<N; ++j) {
      fs << dy*j << " ";
    }

    fs << "\nZ_COORDINATES 1 double\n"
       << "0.0\n"
       << "POINT_DATA " << M*N << "\n"
       << "VECTORS " << descr << " double\n";

    fs << std::scientific;
  
    std::array<int,2> coords_src;
    unsigned first_N, first_M, M_tmp;

    std::vector<double> U_tmp(M_max), V_tmp(M_max);
    std::vector<std::vector<int>> M_loc_ij(dims[0], std::vector<int>(dims[1]));
    std::vector<int> N_tmp(dims[1]);

    for(int j=0; j<dims[1]; ++j){
      if(coords[0]==0 && coords[1]==j){
	N_tmp[j] = N_loc;
      }else{
	MPI_Recv((void*)&N_tmp[j], 1, MPI_UNSIGNED, get_MPI_Cart_rank(comm_grid, std::array<int,2>{0,j}.data()), j, comm_grid, 
		 &status);
      }
    }
    
    for(int i=0; i<dims[0]; ++i){
      for(int j=0; j<dims[0]; ++j){
	if(coords[0]==i && coords[1]==j){
	  M_loc_ij[i][j] = M_loc;
	}else{
	  int rank_src = get_MPI_Cart_rank(comm_grid, std::array<int,2>{i,j}.data());
	  MPI_Recv((void*)&M_loc_ij[i][j], 1, MPI_UNSIGNED, rank_src, rank_src, comm_grid, &status);
	}
      }
    }
      
    for(int j=0; j<dims[1]; ++j){      
      //Dont copy the first or last row if we are not on the boundary
      first_N = j==0 ? 0 : 1;
      N_tmp[j] = j==(dims[1]-1) ? N_tmp[j] : N_tmp[j]-1;
      for(int k=first_N; k<N_tmp[j]; ++k){
	
      	for(int i=0; i<dims[0]; ++i){   
      	  if(coords[0]==i && coords[1]==j){
      	    std::copy_n(U.begin()+k*M_loc, M_loc, U_tmp.begin());
      	  }else{
	    int rank_src = get_MPI_Cart_rank(comm_grid, std::array<int,2>{i,j}.data());
      	    MPI_Recv((void*)U_tmp.data(), M_loc_ij[i][j], MPI_DOUBLE, rank_src, nprocs*nprocs*rank_src+nprocs*k, comm_grid, 
		     &status);
      	    MPI_Recv((void*)V_tmp.data(), M_loc_ij[i][j], MPI_DOUBLE, rank_src, nprocs*nprocs*nprocs*rank_src+nprocs*k, comm_grid, 
		     &status);
      	  }	  
      	  first_M = i==0 ? 0 : 1;
      	  M_tmp = i==(dims[0]-1) ? M_loc_ij[i][j] : M_loc_ij[i][j]-1;
      	  for(int l=first_M; l<M_tmp; ++l){
      	    fs << U_tmp[l] << " " << V_tmp[l] << " 0.0" << std::endl;
      	  }	  
      	}
      }
    }

    fs << std::endl;  
    fs.close();
  }else{
    if(coords[0]==0){
      MPI_Send((void*)&N_loc, 1, MPI_UNSIGNED, 0, coords[1], comm_grid);
    }
    MPI_Send((void*)&M_loc, 1, MPI_UNSIGNED, 0, rank, comm_grid);

    unsigned first_N = coords[1]==0 ? 0 : 1;
    unsigned N_tmp = coords[1]==(dims[1]-1) ? N_loc : N_loc-1;

    for(unsigned k=first_N; k<N_tmp; ++k){
      MPI_Send((void*)(U.begin()+k*M_loc), M_loc, MPI_DOUBLE, 0, nprocs*nprocs*rank+nprocs*k, comm_grid);
      MPI_Send((void*)(V.begin()+k*M_loc), M_loc, MPI_DOUBLE, 0, nprocs*nprocs*nprocs*rank+nprocs*k, comm_grid);
    }
  }
}
