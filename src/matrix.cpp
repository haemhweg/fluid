#include <cassert>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>

#include "real.h"
#include "matrix.h"
#include "vector.h"

void Matrix::allocateData()
{
  _data = new REAL[M*N];
}

Matrix::Matrix(const size_t M_, const size_t N_, const REAL val) 
  : _data(nullptr), M(M_), N(N_)
{
  allocateData();
  fill(val);
}

Matrix::Matrix(Matrix&& rhs) : _data(rhs._data), M(rhs.M), N(rhs.N)
{
  rhs._data = nullptr;
  rhs.M = 0;
  rhs.N = 0;
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
    for(size_t j=0; j<N; ++j) {
      fs.read((char *)&_data[j*M+i], real_size);
    }
  }
}

Matrix::~Matrix()
{
  if(N && M) delete [] _data;
}

void Matrix::fill(REAL val)
{
  std::fill_n(_data, M*N, val);
}


REAL Matrix::at(size_t i, size_t j) const
{
  return _data[j*M+i];
}

REAL& Matrix::at(size_t i, size_t j)
{
  return _data[j*M+i];
}

size_t Matrix::getRows() const
{
  return M;
}

size_t Matrix::getCols() const
{
  return N;
}

void Matrix::print(const std::string& descr) const
{
  std::cout << std::endl;
  std::cout << descr << std::fixed << std::setprecision(1) << std::setw(7) 
	    << std::left<< std::endl;

  for (size_t i = 0; i < M; i++)
    {
      for (size_t j = 0; j < N; j++)
	{
	  std::cout << _data[j*M+i] << "\t";
	}
      std::cout << std::endl;
    }

  std::cout << std::setprecision(6);
}

void Matrix::writeBinary(const std::string & filename) const
{
  std::ofstream fs(filename, std::ios::binary);
  
  size_t real_size = sizeof(REAL);

  fs.write((char *)&M, sizeof(size_t));
  fs.write((char *)&N, sizeof(size_t));

  fs.write((char *)_data, M*N*real_size);
}

void Matrix::writeVTK(const std::string & filename, const std::string& descr, const double dx, const double dy) const
{
  std::ofstream fs(filename);

  fs << std::fixed << std::setprecision(6);

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
  
  fs << std::scientific;

  for(size_t i=0; i<M*N; ++i) {
    fs << _data[i] << "\n";
  }

  fs << std::endl;    
}

void writeVectorFieldVTK(const std::string& filename, const std::string& descr, 
			 const Matrix& U, const Matrix& V, const double dx, 
			 const double dy)
{
  std::ofstream fs(filename);

  size_t M = U.getRows(), N = U.getCols();

  assert(M == V.getRows() && N == V.getCols());

  fs << std::fixed<< std::setprecision(6);

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

  fs << std::scientific;
  
  const REAL* u = U.begin();
  const REAL* v = V.begin();
  for(size_t i=0; i<M*N; ++i) {
    fs << u[i] << " " << v[i] << " 0.0\n";
  }

  fs << std::endl;  
}

bool operator==(const Matrix & A1, const Matrix & A2)
{
  assert(A1.getRows() == A2.getRows() && A1.getCols() == A2.getCols());
  
  return std::equal(std::begin(A1), std::end(A1), std::begin(A2), REAL_equal);
}

Vector operator*(const Matrix & A, const Vector & v)
{
  assert(v.getSize() == A.getCols());

  Vector r{v.getSize()};

  for (size_t i = 0; i < A.getRows(); i++)
    {
      for (size_t j = 0; j < A.getCols(); j++)
	{
	  r.at(i) += A.at(i, j) * v.at(j);
	}
    }
  
  return r;
}

Matrix operator*(REAL a, const Matrix & A)
{
  Matrix B{A.getRows(), A.getCols()};

  for (size_t i = 0; i < A.getRows(); i++)
    {
      for (size_t j = 0; j < A.getCols(); j++)
	{
	  B.at(i, j) =  a * A.at(i, j);
	}
    }

  return B;
}

Matrix operator+(const Matrix & A1, const Matrix & A2)
{
  assert(A1.getRows() == A2.getRows() && A1.getCols() == A2.getCols());

  Matrix B{A1.getRows(), A1.getCols()};

  for (size_t i = 0; i < A1.getRows(); i++)
    {
      for (size_t j = 0; j < A1.getCols(); j++)
	{
	  B.at(i, j) = A1.at(i, j) + A2.at(i, j);
	}
    }

  return B;
}
