#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <functional>

#include "real.h"

class Vector
{
private:
	/* It tends to be a bad idea to make non-const fields public because it then becomes hard to force error checking constraints and/or add side-effects to value changes in the future. */
	REAL* _data;
	size_t size;
	void allocateData();
public:
	Vector(const size_t size, const REAL val=REAL(0.));
	Vector(Vector&& rhs);
	~Vector();
	/* Fill with REAL values */
	void fill(const REAL value);

	size_t getSize() const;

	/* Since _data is private, use getters and setters */
	REAL at(const size_t N) const;
	REAL& at(const size_t N);

	REAL* begin() { return _data; }
	REAL* end() { return _data+size; }
	const REAL* begin() const { return _data; }
	const REAL* end() const { return _data+size; }

	REAL nrm2() const;

	void print(const std::string& descr) const;
	void write(const std::string & filename) const;
	void writeVTKfile(const std::string & filename) const;
};

bool operator==(const Vector & v1, const Vector & v2);
Vector operator*(const REAL a, const Vector & v);
Vector operator+(const Vector & v1, const Vector & v2);
REAL operator*(const Vector & v1, const Vector & v2);

#endif
