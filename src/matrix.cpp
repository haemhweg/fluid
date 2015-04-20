#include <stdlib.h>
#include <time.h>
#include <stdexcept>
#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <ios>

#include "real.h"
#include "matrix.h"
#include "vector.h"

void Matrix::allocateData()
{
  _data = new REAL*[M];
  for (size_t i = 0; i < M; i++)
    {
      _data[i] = new REAL[N];
    }
}

Matrix::Matrix(size_t M_, size_t N_) : M(M_), N(N_)
{
  allocateData();
}

Matrix::Matrix(const std::string & filename)
{
  std::ifstream fs(filename, std::ios::binary);

  size_t real_size = sizeof(REAL);

  fs.read((char *)&M, sizeof(size_t));
  fs.read((char *)&N, sizeof(size_t));

  std::cout << M << " " << N << std::endl;

  allocateData();

  for(size_t i=0; i<M; ++i) {
    fs.read((char *)_data[i], N*real_size);
  }
}

Matrix::~Matrix()
{
  for(size_t i=0; i<M; ++i) {
    delete [] _data[i];
  }
  delete _data;
}

void Matrix::fill(REAL value)
{
  for (size_t i = 0; i < M; i++)
    {
      for (size_t j = 0; j < N; j++)
	{
	  _data[i][j] = value;
	}
    }
}

void Matrix::fillrand()
{
  srand(time(NULL));
  for (size_t i = 0; i < M; i++)
    {
      for (size_t j = 0; j < N; j++)
	{
	  _data[i][j] = rand();
	}
    }
}

bool operator==(const Matrix & v1, const Matrix & v2)
{

  if (v1.getRows() != v2.getRows() || v1.getCols() != v2.getCols())
    {
      throw std::runtime_error("Illegal comparison of matrices");
    }

  for (size_t i = 0; i < v1.getRows(); i++)
    {
      for (size_t j = 0; j < v1.getCols(); j++)
	{
	  if (v1.get(i, j) != v2.get(i, j))
	    {
	      return false;
	    }
	}
    }
  return true;
}

REAL Matrix::get(size_t M_, size_t N_) const
{
  return _data[M_][N_];
}

void Matrix::set(size_t M_, size_t N_, REAL v)
{
  _data[M_][N_] = v;
}

size_t Matrix::getRows() const
{
  return M;
}

size_t Matrix::getCols() const
{
  return N;
}

void Matrix::apply(const std::function<REAL(REAL)> & f)
{
  for (size_t i = 0; i < M; i++)
    {
      for (size_t j = 0; j < N; j++)
	{
	  _data[i][j] = f(_data[i][j]);
	}
    }
}

void Matrix::print()
{
  std::cout << std::endl;
  std::cout << std::setprecision(2) << std::setw(7);

  for (size_t i = 0; i < M; i++)
    {
      for (size_t j = 0; j < N; j++)
	{
	  std::cout << _data[i][j] << "\t";
	}
      std::cout << std::endl;
    }
}

void Matrix::write(const std::string & filename)
{
  std::ofstream fs(filename, std::ios::binary);
  
  size_t real_size = sizeof(REAL);

  fs.write((char *)&M, sizeof(size_t));
  fs.write((char *)&N, sizeof(size_t));

  for (size_t i = 0; i < M; i++)
    {
      fs.write((char *)_data[i], N*real_size);
    }
}

void Matrix::writeVTKfile(const std::string & filename, const std::string& descr, const double dx, const double dy)
{
  std::ofstream fs(filename);

  fs << "# vtk DataFile Version 3.0\n"
     << "Scalar Field\n"
     << "ASCII\n"
     << "DATASET RECTILINEAR_GRID\n"
     << "DIMENSIONS " << M << " " << N << " 1\n"
     << "X_COORDINATES " << M << " double\n";

  for(size_t i=0; i<M; ++i) {
    fs << dx*i << " ";
  }

  fs << "\nY_COORDINATES " << N << " double\n";
  
  for(size_t j=0; j<N; ++j) {
    fs << dy*j << " ";
  }

  fs << "\nZ_COORDINATES 1 double\n"
     << "0.0\n"
     << "POINT_DATA " << M*N << "\n"
     << "SCALARS " << descr << " double 1\n"
     << "LOOKUP_TABLE default\n";

  for(size_t j=0; j<N; ++j) {
    for(size_t i=0; i<M; ++i) {
      fs << _data[i][j] << "\n";
    }
  }

  fs << std::endl;    
}

void writeVectorFieldVTK(const std::string& filename, const std::string& descr, 
			 const Matrix& U, const Matrix& V, const double dx, 
			 const double dy)
{
  std::ofstream fs(filename);

  size_t M = U.getRows(), N = U.getCols();

  fs << "# vtk DataFile Version 3.0\n"
     << "Vector Field\n"
     << "ASCII\n"
     << "DATASET RECTILINEAR_GRID\n"
     << "DIMENSIONS " << M << " " << N << " 1\n"
     << "X_COORDINATES " << M << " double\n";

  for(size_t i=0; i<M; ++i) {
    fs << dx*i << " ";
  }

  fs << "\nY_COORDINATES " << N << " double\n";
  
  for(size_t j=0; j<N; ++j) {
    fs << dy*j << " ";
  }

  fs << "\nZ_COORDINATES 1 double\n"
     << "0.0\n"
     << "POINT_DATA " << M*N << "\n"
     << "VECTORS " << descr << " double\n";

  for(size_t j=0; j<N; ++j) {
    for(size_t i=0; i<M; ++i) {
      fs << U.get(i,j) << " " << V.get(i,j) << " 0.0\n";
    }
  }

  fs << std::endl;  
}

Vector * operator*(const Matrix & A, const Vector & v)
{

  if (v.getSize() != A.getCols())
    {
      throw std::runtime_error("Illegal operation on incompatible structures.");
    }

  Vector * r = new Vector(v.getSize());

  for (size_t i = 0; i < A.getRows(); i++)
    {
      REAL p = 0;
      for (size_t j = 0; j < A.getCols(); j++)
	{
	  p += A.get(i, j) * v.get(j);
	}
      r->set(i, p);
    }

  return r;

}

Matrix * operator*(REAL a, const Matrix & A)
{

  Matrix * B = new Matrix(A.getRows(), A.getCols());

  for (size_t i = 0; i < A.getRows(); i++)
    {
      for (size_t j = 0; j < A.getCols(); j++)
	{
	  B->set(i, j, a * A.get(i, j));
	}
    }

  return B;

}

Matrix * operator+(const Matrix & A1, const Matrix & A2)
{

  if (A1.getRows() != A2.getRows() || A1.getCols() != A2.getCols())
    {
      throw std::runtime_error("Illegal operation on incompatible matrices.");
    }

  Matrix * B = new Matrix(A1.getRows(), A1.getCols());

  for (size_t i = 0; i < A1.getRows(); i++)
    {
      for (size_t j = 0; j < A1.getCols(); j++)
	{
	  B->set(i, j, A1.get(i, j) + A2.get(i, j));
	}
    }

  return B;

}
