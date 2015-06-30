#ifndef PARALLEL_H
#define PARALLEL_H
#include <array>
#include <string>
#include "mpi.h"

class Matrix;

int get_MPI_Comm_size(const MPI_Comm& comm);
int get_MPI_Comm_rank(const MPI_Comm& comm);
int get_MPI_Cart_rank(const MPI_Comm& comm, const int* coords);
std::array<int,2> get_MPI_Dims_create(const MPI_Comm& comm, int ndims);
std::array<int,2> get_MPI_Cart_coords(const MPI_Comm& comm_grid, int ndims);
std::array<int,2> get_MPI_Cart_dims(const MPI_Comm& comm);


int init_MPI_Grid(MPI_Comm& comm_grid);

void Matrix_exchange(const MPI_Comm& comm_grid, Matrix& M);
void MPI_VectorFieldVTK(const MPI_Comm& comm_grid, const std::string& filename, const std::string& descr, 
			const Matrix& U, const Matrix& V, const double dx, const double dy);

#endif
