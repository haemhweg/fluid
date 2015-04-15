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
	Vector(size_t size);
	/**
	* Since we write size of the matrix in the file, we can't instaniate it on its own without reading it before.
	* So instead of read function we use appropriate constructor.
	*/
	Vector(const std::string & filename);
	~Vector();
	/* Fill with REAL values */
	void fill(REAL value);
	/* To make the interface more abstract (for instance, to enable usage of lambda functions), use std::function (exhausting workaround with templates was futile). */
	void apply(const std::function<REAL(REAL)> &);
	size_t getSize() const;
	/* Copy vector x into this vector */
	void copy(const Vector & x);
	/* Since _data is private, use getters and setters */
	REAL get(size_t N) const;
	void set(size_t N, REAL v);
	REAL nrm2() const;
	void print();
	void write(const std::string & filename);
	void writeVTKfile(const std::string & filename);

};

bool operator==(const Vector & v1, const Vector & v2);
Vector * operator*(REAL a, const Vector & v);
Vector * operator+(const Vector & v1, const Vector & v2);
REAL operator*(const Vector & v1, const Vector & v2);

#endif