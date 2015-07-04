#include <cassert>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>

#include "real.h"
#include "field.h"

void Field2D::allocateData()
{
  _data = new REAL[(imax_+2)*(jmax_+2)];
}

Field2D::Field2D(const unsigned imax, const unsigned jmax, const REAL val) 
  : _data(nullptr), imax_(imax), jmax_(jmax)
{
  allocateData();
  fill(val);
}

Field2D::Field2D(Field2D&& rhs) : _data(rhs._data), imax_(rhs.imax_), jmax_(rhs.jmax_)
{
  rhs._data = nullptr;
  rhs.imax_ = 0;
  rhs.jmax_ = 0;
}
Field2D::Field2D(const Field2D& rhs) : _data(), imax_(rhs.imax_), jmax_(rhs.jmax_)
{
  allocateData();
  std::copy_n(rhs._data, (imax_+2)*(jmax_+2), _data);
}

Field2D& Field2D::operator=(Field2D&& rhs)
{
  if(_data) delete [] _data;

  _data = rhs._data;
  imax_ = rhs.imax_;
  jmax_ = rhs.jmax_;

  rhs._data = nullptr;
  rhs.imax_ = 0;
  rhs.jmax_ = 0;

  return *this;
}

Field2D::~Field2D()
{
  if(imax_ && jmax_) delete [] _data;
}

void Field2D::fill(REAL val)
{
  std::fill_n(_data, (imax_+2)*(jmax_+2), val);
}


REAL Field2D::at(unsigned i, unsigned j) const
{
  assert(i<imax_+2 && j<jmax_+2);
  return _data[j*(imax_+2)+i];
}

REAL& Field2D::at(unsigned i, unsigned j)
{
  assert(i<imax_+2 && j<jmax_+2);
  return _data[j*(imax_+2)+i];
}

unsigned Field2D::imax() const
{
  return imax_;
}

unsigned Field2D::jmax() const
{
  return jmax_;
}

REAL Field2D::getMax() const
{
  return *std::max_element(_data, _data + (imax_+2)*(jmax_+2));
}

void Field2D::print(const std::string& descr) const
{
  std::cout << std::endl;
  std::cout << descr << std::fixed << std::setprecision(1) << std::setw(7) 
	    << std::left<< std::endl;

  for(int j=jmax_+1; j>=0; --j){
    for(unsigned i=0; i<imax_+2; ++i){
      std::cout << at(i,j) << "\t";
    }
    std::cout << std::endl;
  }

  std::cout << std::setprecision(6);
}

void Field2D::writeBinary(const std::string & filename) const
{
  std::ofstream fs(filename, std::ios::binary);
  
  const unsigned real_size = sizeof(REAL);

  const unsigned M = imax_+2;
  const unsigned N = jmax_+2;

  fs.write((char *)&M, sizeof(unsigned));
  fs.write((char *)&N, sizeof(unsigned));

  fs.write((char *)_data, M*N*real_size);
}

void Field2D::writeVTK(const std::string & filename, const std::string& descr, const double dx, const double dy) const
{
  std::ofstream fs(filename);

  fs << std::fixed << std::setprecision(6);

  fs << "# vtk DataFile Version 3.0\n"
     << "Scalar Field2D\n"
     << "ASCII\n"
     << "DATASET RECTILINEAR_GRID\n"
     << "DIMENSIONS " << imax_+2 << " " << jmax_+2 << " 1\n"
     << "X_COORDINATES " << imax_+2 << " double\n";

  for(unsigned i=0; i<imax_+2; ++i) {
    fs << dx*i << " ";
  }

  fs << "\nY_COORDINATES " << jmax_+2 << " double\n";
  
  for(unsigned j=0; j<jmax_+2; ++j) {
    fs << dy*j << " ";
  }

  fs << "\nZ_COORDINATES 1 double\n"
     << "0.0\n"
     << "POINT_DATA " << (imax_+2)*(jmax_+2) << "\n"
     << "SCALARS " << descr << " double 1\n"
     << "LOOKUP_TABLE default\n";
  
  fs << std::scientific;

  for(unsigned i=0; i<(imax_+2)*(jmax_+2); ++i) {
    fs << _data[i] << "\n";
  }

  fs << std::endl;    
}

void writeVectorFieldVTK(const std::string& filename, const std::string& descr, 
			 const Field2D& U, const Field2D& V, const double dx, 
			 const double dy)
{
  std::ofstream fs(filename);

  const unsigned imax = U.imax(), jmax = U.jmax();

  assert(imax == V.imax() && jmax == V.jmax());

  fs << std::fixed<< std::setprecision(6);

  fs << "# vtk DataFile Version 3.0\n"
     << "Vector Field\n"
     << "ASCII\n"
     << "DATASET RECTILINEAR_GRID\n"
     << "DIMENSIONS " << imax+2 << " " << jmax+2 << " 1\n"
     << "X_COORDINATES " << imax+2 << " double\n";

  for(unsigned i=0; i<imax+2; ++i) {
    fs << dx*i << " ";
  }

  fs << "\nY_COORDINATES " << jmax+2 << " double\n";
  
  for(unsigned j=0; j<jmax+2; ++j) {
    fs << dy*j << " ";
  }

  fs << "\nZ_COORDINATES 1 double\n"
     << "0.0\n"
     << "POINT_DATA " << (imax+2)*(jmax+2) << "\n"
     << "VECTORS " << descr << " double\n";

  fs << std::scientific;
  
  const REAL* u = U.data();
  const REAL* v = V.data();
  for(unsigned i=0; i<(imax+2)*(jmax+2); ++i) {
    fs << u[i] << " " << v[i] << " 0.0\n";
  }

  fs << std::endl;  
}

