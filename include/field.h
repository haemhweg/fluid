#ifndef __FIELD_H__
#define __FIELD_H__

#include <functional>
#include <string>

#include "real.h"

class Field2D
{
private:
  
  REAL* _data; //!< Enträge sind Zeilenweise gespeichert
  unsigned imax_, jmax_;
  void allocateData();
public:
  Field2D(const unsigned imax, const unsigned jmax, const REAL val=REAL(0.));
  Field2D(Field2D&& rhs);
  Field2D(const Field2D& rhs);
  Field2D() : Field2D(0,0,0) { }

  ~Field2D();

  Field2D & operator=(Field2D&&);

  void fill(const REAL value);

  /**
   * Since we want to encapsulate _data, we need getters and setters. Is design alright (don't know whether C++-compatible)?
   */
  REAL at(const unsigned i, const unsigned j) const;
  REAL& at(const unsigned i, const unsigned j);

  unsigned imax() const;
  unsigned jmax() const;

  REAL* data() { return _data; }
  const REAL* data() const { return _data; }

  void print(const std::string& descr) const;
  void writeVTK(const std::string & filename, const std::string& descr, 
		const double dx=0.1, const double dy=0.1) const;

  REAL getMax() const;

};

/**
 *  Writes a 2D vector field to a .vtk file. The first component of the vector field is passed by the matrix
 *  U and the second component by the matrix V.
 *  The filename has to end in ".vtk"!
 */
void writeVectorFieldVTK(const std::string& filename, const std::string& descr, const Field2D& U, const Field2D& V, 
			 const double dx=0.1, const double dy=0.1);

#endif
