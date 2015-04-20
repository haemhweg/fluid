#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <functional>
#include <string>

#include "real.h"
#include "vector.h"

/**
 * So I've created only matrix and no field method. Does it matter?
 */
class Matrix
{
private:
  
  REAL* _data; //!< Enträge sind Spaltenweise gespeichert
  size_t M, N;
  void allocateData();
public:
  Matrix(size_t M, size_t N, const REAL val=REAL(0.));
  Matrix(Matrix&& rhs);
  /**
   * Since we write size of the matrix in the file, we can't instaniate it on its own without reading it before.
   * So instead of read function we use appropriate constructor.
   */
  Matrix(const std::string & filename);
  ~Matrix();

  void fill(REAL value);

  /**
   * Since we want to encapsulate _data, we need getters and setters. Is design alright (don't know whether C++-compatible)?
   */
  REAL at(size_t i, size_t j) const;
  REAL& at(size_t i, size_t j);

  size_t getRows() const;
  size_t getCols() const;

  REAL* begin() { return _data; }
  REAL* end() { return _data+N*M; }
  const REAL* begin() const { return _data; }
  const REAL* end() const {return _data+N*M; }

  void print(const std::string& descr) const;
  void writeBinary(const std::string & filename) const;
  void writeVTK(const std::string & filename, const std::string& descr, 
		const double dx=0.1, const double dy=0.1) const;
};

/**
 * To make operations look more descriptive, we overload operators. The question here is, whether we should define these operations so that they return void?
 * Like we would perform those operation on the subject itself: for instance, void operator*(REAL a, Matrix & A).
 * This would spare memory in case we don't need a copy of vector / matrix. Same consideration for vectors.
 */

Vector operator*(const Matrix & A, const Vector & v);
Matrix operator*(REAL a, const Matrix & A);
Matrix operator+(const Matrix & A1, const Matrix & A2);
bool operator==(const Matrix & A1, const Matrix & A2);

/**
 *  Writes a 2D vector field to a .vtk file. The first component of the vector field is passed by the matrix
 *  U and the second component by the matrix V.
 *  The filename has to end in ".vtk"!
 */
void writeVectorFieldVTK(const std::string& filename, const std::string& descr, const Matrix& U, const Matrix& V, const double dx=0.1, const double dy=0.1);

#endif
