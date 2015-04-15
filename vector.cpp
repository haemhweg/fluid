#include <string>
#include <iostream>
#include <stdlib.h>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cmath>

#include "real.h"
#include "vector.h"

void Vector::allocateData()
{
	_data = (REAL *)malloc(sizeof(REAL) * size);
}

<<<<<<< HEAD
Vector::Vector(size_t size) : size(size)
{
	allocateData();
}

Vector::~Vector()
{
=======
Vector::Vector(Vector&& rhs) : _data(rhs._data), size(rhs.size) {
  rhs.size = 0;
  rhs._data = nullptr;
}

Vector::~Vector() {
>>>>>>> 82a3f4f7a3a60d87a8d660b2399ecbe03886241c
	free(_data);
}

void Vector::fill(REAL value)
{
	for (size_t i = 0; i < size; i++)
	{
		_data[i] = value;
	}
}

bool operator==(const Vector & v1, const Vector & v2)
{

	if (v1.getSize() != v2.getSize())
	{
		throw std::runtime_error("Illegal comparison of incompatible vectors");
	}

	for (size_t i = 0; i < v1.getSize(); i++)
	{
		if (v1.get(i) != v2.get(i))
		{
			return false;
		}
	}

	return true;

}

Vector * operator*(REAL a, const Vector & v)
{

	Vector * r = new Vector(v.getSize());

	for (size_t i = 0; i < v.getSize(); i++)
	{
		r->set(i, v.get(i) * a);
	}

	return r;

}

REAL operator*(const Vector & v1, const Vector & v2)
{

	if (v1.getSize() != v2.getSize())
	{
		throw std::runtime_error("Illegal operation on incompatible vectors");
	}

	REAL sum = 0;

	for (size_t i = 0; i < v1.getSize(); i++)
	{
		sum += v1.get(i) * v2.get(i);
	}

	return sum;

}

Vector * operator+(const Vector & v1, const Vector & v2)
{
	
	if (v1.getSize() != v2.getSize())
	{
		throw std::runtime_error("Illegal operation on incompatible vectors");
	}

	Vector * r = new Vector(v1.getSize());

	for (size_t i = 0; i < v1.getSize(); i++)
	{
		r->set(i, v1.get(i) + v2.get(i));
	}

	return r;

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

void Vector::apply(const std::function<REAL(REAL)> & f)
{
	for (size_t i = 0; i < size; i++)
	{
		_data[i] = f(_data[i]);
	}
}

REAL Vector::get(size_t N) const
{
	return _data[N];
}

void Vector::set(size_t N, REAL v)
{
	_data[N] = v;
}

void Vector::copy(const Vector & x)
{

	if (size != x.size)
	{
		throw std::runtime_error("Illegal copy of incompatible vectors.");
	}

	for (size_t i = 0; i < size; i++)
	{
		_data[i] = x._data[i];
	}

}

<<<<<<< HEAD
void Vector::print()
{

	std::cout << std::endl;
	std::cout << std::setprecision(2) << std::setw(7);

	for (size_t i = 0; i < size; i++)
	{
		std::cout << _data[i] << std::endl;
	}

}

void Vector::write(const std::string & filename)
{
  std::ofstream fs(filename, std::ios::binary);

	fs << size << std::endl;

	for (size_t i = 0; i < size; i++)
	{
		fs << _data[i] << std::endl;
	}
}

Vector::Vector(const std::string & filename)
{
  std::ifstream fs(filename, std::ios::binary);

	fs >> size;

	allocateData();

	REAL d;

	size_t i = 0;
	while (fs >> d)
	{
		_data[i++] = d;
	}
=======
void Vector::print(const std::string& prefix, const std::string& delim) const {
  std::cout << prefix << std::endl;

  std::cout << std::setprecision(2);
  for(unsigned i=0; i<size; ++i)
    std::cout << "\t" << _data[i] << delim;

  std::cout << std::endl;    
}

void Vector::writeBinary(const std::string& fileName) const {
  std::ostream file(fileName);

  file << size;

  for(unsigned i=0; i<size; ++i) {
    file << " " << _data[i];
  }
  file << std::endl;    
}

Vector readVectorFromBinary(const std::string& fileName) {
  std::ifstream file{fileName};
  size_t size;

  file >> size;

  Vector x{size};

  for(size_t i=0; i<size; ++i) {
    REAL val;
    file >> val;
    x.set(i, val);
  }

  return x;
>>>>>>> 82a3f4f7a3a60d87a8d660b2399ecbe03886241c
}
