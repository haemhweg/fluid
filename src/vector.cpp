#include <string>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <cassert>

#include "real.h"
#include "vector.h"

void Vector::allocateData()
{
  _data = new REAL[size];
}

Vector::Vector(const size_t size_, const REAL val) : _data(nullptr), size(size_)
{
  allocateData();
  fill(val);
}

Vector::Vector(Vector&& rhs) : _data(rhs._data), size(rhs.size)
{
  rhs._data = nullptr;
  rhs.size = 0;
}

Vector::~Vector()
{
  if(size) delete [] _data;
}

void Vector::fill(REAL val)
{
  std::fill_n(_data, size, val);
}



REAL Vector::nrm2() const
{	
  /* Kinda works */
  return sqrt((*this) * (*this));
}

size_t Vector::getSize() const
{
  return size;
}

REAL Vector::at(size_t i) const
{
  return _data[i];
}

REAL& Vector::at(size_t i)
{
  return _data[i];
}

void Vector::print(const std::string& descr) const
{
  std::cout << std::endl;
  std::cout << descr << std::fixed << std::setprecision(1) << std::setw(7) 
	    << std::left<< std::endl;;

  for (size_t i = 0; i < size; i++)
    {
      std::cout << "\t" << _data[i] << std::endl;
    }
}

bool operator==(const Vector & v1, const Vector & v2)
{
  assert(v1.getSize() == v2.getSize());

  for (size_t i = 0; i < v1.getSize(); i++)
    {
      if (v1.at(i) != v2.at(i))
	{
	  return false;
	}
    }

  return true;

}

Vector operator*(REAL a, const Vector & v)
{

  Vector r{v.getSize()};

  for (size_t i = 0; i < v.getSize(); i++)
    {
      r.at(i) = v.at(i) * a;
    }

  return r;

}

REAL operator*(const Vector & v1, const Vector & v2)
{
  assert(v1.getSize() == v2.getSize());

  REAL sum = 0;

  for (size_t i = 0; i < v1.getSize(); i++)
    {
      sum += v1.at(i) * v2.at(i);
    }

  return sum;
}

Vector operator+(const Vector & v1, const Vector & v2)
{	
  assert(v1.getSize() == v2.getSize());

  Vector r{v1.getSize()};

  for (size_t i = 0; i < v1.getSize(); i++)
    {
      r.at(i) = v1.at(i) + v2.at(i);
    }

  return r;
}
