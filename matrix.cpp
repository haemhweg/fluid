#include <stdlib.h>
#include <time.h>
#include <stdexcept>
#include <functional>
<<<<<<< HEAD
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <ios>
=======
#include <string>
#include <iostream>
>>>>>>> 82a3f4f7a3a60d87a8d660b2399ecbe03886241c

#include "real.h"
#include "matrix.h"
#include "vector.h"

void Matrix::allocateData()
{
	_data = (REAL **)malloc(sizeof(REAL *) * M);
	for (size_t i = 0; i < M; i++)
	{
		_data[i] = (REAL *)malloc(sizeof(REAL) * N);
	}
}

<<<<<<< HEAD
Matrix::Matrix(size_t M, size_t N) : M(M), N(N)
{
	allocateData();
=======
Matrix::Matrix(Matrix&& rhs) : _data(rhs._data), M(rhs.M), N(rhs.N) {
  rhs.M = 0;
  rhs.N = 0;
  rhs._data = nullptr;
}

Matrix::~Matrix() {
  for(size_t i=0; i<M; ++i) {
    free(_data[i]);
  }
  free(_data);
>>>>>>> 82a3f4f7a3a60d87a8d660b2399ecbe03886241c
}

Matrix::Matrix(const std::string & filename)
{
  std::ifstream fs(filename, std::ios::binary);

  size_t real_size = sizeof(REAL);

  fs.read(&M, real_size);
  fs.read(&N, real_size);

  allocateData();

  for(size_t i=0; i<M; ++i) {
    fs.read(_data[i], N*real_size);
  }
}

Matrix::~Matrix()
{
  for(size_t i=0; i<M; ++i) {
    free(_data[i]);
  }
  free(_data);
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

REAL Matrix::get(size_t M, size_t N) const
{
	return _data[M][N];
}

void Matrix::set(size_t M, size_t N, REAL v)
{
	_data[M][N] = v;
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

<<<<<<< HEAD
Vector * operator*(const Matrix & A, const Vector & v)
{
=======
void Matrix::print(const std::string& prefix, const std::string& delim) const {
  std::cout << prefix;

  std::cout << std::setprecision(2);
  for(unsigned i=0; i<M; ++i) {
    for(unsigned j=0; j<N; ++j) {
      std::cout << "\t" << _data[i][j] << delim;
    }
    std::cout << "\n";
  }
  std::cout << std::endl;
}

void Matrix::writeBinary(const std::string& fileName) const {
  std::ofstream file(fileName);

  file << M << " " << N;

  for(unsigned i=0; i<M; ++i) {
    for(unsigned j=0; j<N; ++j) {
      file << " " << _data[i][j];
    }
  }
  file << std::endl;
}

void Matrix::writeVTK(const std::string& fileName, const std::string& descr, const REAL dx, const REAL dy) const {
  std::ofstream file(fileName);

  // Header information
  file << "# vtk DataFile Version 3.0\n"
       << "Scalar Field\n"
       << "ASCII\n"
       << "DATASET RECTILIEAR_GRID\n"
       << "DIMENSIONS " << M << " " << N << " 1\n"
       << "X_COORDINATES " << M << " double\n";
  
  for(size_t i=0; i<M; ++i) {
    file << dx*i << " ";
  }

  file << "\nY_COORDINATES " << N << " double\n";
  
  for(size_t j=0; j<N; ++j) {
    file << dy*j << " ";
  }

  file << "\nZ_COORDINATES 1 double\n"
       << "0.0\n"
       << "POINT_DATA " << M*N << "\n"
       << "SCALARS " << descr << " double\n"
       << "LOOKUP_TABLE default\n";

  // Data
  for(size_t j=0; j<N; ++j) {
    for(size_t i=0; i<M; ++i) {
      file << _data[i][j] << "\n";
    }	
  }
  
  file << std::endl;
}

Vector * operator*(const Matrix & A, const Vector & v) {
>>>>>>> 82a3f4f7a3a60d87a8d660b2399ecbe03886241c

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

<<<<<<< HEAD
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

  fs.write(&M, real_size);
  fs.write(&N, real_size);

  for (size_t i = 0; i < M; i++)
    {
      fs.write(_data[i], N*real_size);
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

void writeVectorFieldVTK(const std::string& filename, const std::string& descr, const Matrix& U, const Matrix& V, const double dx, const double dy)
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
=======
Matrix readMatrixFromBinary(const std::string& fileName) {
  std::ifstream file{fileName};
  size_t M{0}, N{0};

  file >> M >> N;

  Matrix A{M, N};

  for(size_t i=0; i<M; ++i) {
    for(size_t j=0; j<N; ++j) {
      REAL val{0};
      file >> val;
      A.set(i, j, val);
    }
  }

  return A;
}

void writeVectorFieldVTK(const std::string& fileName, const std::string& descr,
			 const Matrix& U, const Matrix& V, const REAL dx, const REAL dy) {
  assert(U.getRows() == V.getRows() && U.getCols() == V.getCols());

   std::ofstream file(fileName);

  // Header information
  file << "# vtk DataFile Version 3.0\n"
       << "Vector Field\n"
       << "ASCII\n"
       << "DATASET RECTILIEAR_GRID\n"
       << "DIMENSIONS " << M << " " << N << " 1\n"
       << "X_COORDINATES " << M << " double\n";
  
  for(size_t i=0; i<M; ++i) {
    file << dx*i << " ";
  }

  file << "\nY_COORDINATES " << N << " double\n";
  
  for(size_t j=0; j<N; ++j) {
    file << dy*j << " ";
  }

  file << "\nZ_COORDINATES 1 double\n"
       << "0.0\n"
       << "POINT_DATA " << M*N << "\n"
       << "VECTORS " << descr << " double\n";

  // Data
  for(size_t j=0; j<N; ++j) {
    for(size_t i=0; i<M; ++i) {
      file << U.get(i,j) << " " << V.get(i,j) << " 0.0\n";
    }	
  }
  
  file << std::endl;
>>>>>>> 82a3f4f7a3a60d87a8d660b2399ecbe03886241c
}
