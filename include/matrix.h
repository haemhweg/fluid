#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <functional>
#include <string>

#include "real.h"

/**
 * So I've created only matrix and no field method. Does it matter?
 */
class Matrix
{
private:
  
  REAL* _data; //!< Enträge sind Spaltenweise gespeichert
  unsigned M, N;
  void allocateData();
public:
  Matrix(const unsigned M, const unsigned N, const REAL val=REAL(0.));
  Matrix(Matrix&& rhs);
  Matrix(const Matrix& rhs);
  Matrix() : Matrix(0, 0, 0) {}
  /**
   * Since we write size of the matrix in the file, we can't instaniate it on its own without reading it before.
   * So instead of read function we use appropriate constructor.
   */
  Matrix(const std::string & filename);
  ~Matrix();

  Matrix & operator=(Matrix&&);

  void fill(const REAL value);

  /**
   * Since we want to encapsulate _data, we need getters and setters. Is design alright (don't know whether C++-compatible)?
   */
  REAL at(const unsigned i, const unsigned j) const;
  REAL& at(const unsigned i, const unsigned j);

  unsigned getRows() const;
  unsigned getCols() const;

  REAL* begin() { return _data; }
  REAL* end() { return _data+N*M; }
  const REAL* begin() const { return _data; }
  const REAL* end() const {return _data+N*M; }

  void print(const std::string& descr) const;
  void writeBinary(const std::string & filename) const;
  void writeVTK(const std::string & filename, const std::string& descr, 
		const double dx=0.1, const double dy=0.1) const;

  REAL getMax() const;

};

/**
 *  Writes a 2D vector field to a .vtk file. The first component of the vector field is passed by the matrix
 *  U and the second component by the matrix V.
 *  The filename has to end in ".vtk"!
 */
void writeVectorFieldVTK(const std::string& filename, const std::string& descr, const Matrix& U, const Matrix& V, const double dx=0.1, const double dy=0.1);

#endif
