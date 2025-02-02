#include <cassert>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>

#include "real.h"
#include "matrix.h"

void Matrix::allocateData()
{
  _data = new REAL[M*N];
}

Matrix::Matrix(const unsigned M_, const unsigned N_, const REAL val) 
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
Matrix::Matrix(const Matrix& rhs) : _data(), M(rhs.M), N(rhs.N)
{
  allocateData();
  std::copy_n(rhs._data, M*N, _data);
}

Matrix& Matrix::operator=(Matrix&& rhs)
{
  if(_data) delete [] _data;

  _data = rhs._data;
  M = rhs.M;
  N = rhs.N;

  rhs._data = nullptr;
  rhs.M = 0;
  rhs.N = 0;

  return *this;
}

Matrix::Matrix(const std::string & filename)
{
  std::ifstream fs(filename, std::ios::binary);

  unsigned real_size = sizeof(REAL);

  fs.read((char *)&M, sizeof(unsigned));
  fs.read((char *)&N, sizeof(unsigned));

  std::cout << M << " " << N << std::endl;

  allocateData();

  for(unsigned i=0; i<M; ++i) {
    for(unsigned j=0; j<N; ++j) {
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


REAL Matrix::at(unsigned i, unsigned j) const
{
  return _data[j*M+i];
}

REAL& Matrix::at(unsigned i, unsigned j)
{
  return _data[j*M+i];
}

unsigned Matrix::getRows() const
{
  return M;
}

unsigned Matrix::getCols() const
{
  return N;
}

REAL Matrix::getMax() const
{
  const auto absLess = [] (const REAL a, const REAL b) -> bool 
    { return std::fabs(a)<std::fabs(b); };
  return std::fabs(*std::max_element(_data, _data + M*N, absLess));
}


void Matrix::print(const std::string& descr) const
{
  std::cout << std::endl;
  std::cout << descr << std::fixed << std::setprecision(1) << std::setw(7) 
	    << std::left<< std::endl;

  for (unsigned i = 0; i < M; i++)
    {
      for (unsigned j = 0; j < N; j++)
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
  
  unsigned real_size = sizeof(REAL);

  fs.write((char *)&M, sizeof(unsigned));
  fs.write((char *)&N, sizeof(unsigned));

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
     << "SCALARS " << descr << " double 1\n"
     << "LOOKUP_TABLE default\n";
  
  fs << std::scientific;

  for(unsigned i=0; i<M*N; ++i) {
    fs << _data[i] << "\n";
  }

  fs << std::endl;    
}

void writeVectorFieldVTK(const std::string& filename, const std::string& descr, 
			 const Matrix& U, const Matrix& V, const double dx, 
			 const double dy)
{
  std::ofstream fs(filename);

  unsigned M = U.getRows(), N = U.getCols();

  assert(M == V.getRows() && N == V.getCols());

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
  
  const REAL* u = U.begin();
  const REAL* v = V.begin();
  for(unsigned i=0; i<M*N; ++i) {
    fs << u[i] << " " << v[i] << " 0.0\n";
  }

  fs << std::endl;  
}
