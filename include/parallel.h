#ifndef PARALLEL_H
#define PARALLEL_H
#include <vector>
#include "mpi.h"

class Matrix;

int get_MPI_Comm_size(MPI_Comm comm);
int get_MPI_Comm_rank(MPI_Comm comm);
int get_MPI_Cart_rank(MPI_Comm comm, const int* coords);
std::vector<int> get_MPI_Dims_create(MPI_Comm comm, int ndims);
std::vector<int> get_MPI_Cart_coords(MPI_Comm comm_grid, int ndims);


int init_MPI_Grid(MPI_Comm& comm_grid);

void Matrix_exchange(const MPI_Comm& comm_grid, Matrix& M);

#endif
