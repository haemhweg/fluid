#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <functional>

#include "real.h"
#include "vector.h"

class Matrix {
private:
	REAL** _data;
	size_t M, N;
public:
	Matrix(size_t M, size_t N);
	~Matrix();
	void fill(REAL value);
	REAL get(size_t M, size_t N) const;
	void set(size_t M, size_t N, REAL v);
	size_t getRows() const;
	size_t getCols() const;
	/* To make the interface more abstract (for instance, to enable usage of lambda functions), use std::function (exhausting workaround with templates was futile). */
	void apply(const std::function<REAL(REAL)> &);
};

Vector * operator*(const Matrix & A, const Vector & v);
/* Maybe aX should be void operator*(REAL a, Matrix & A)? */
Matrix * operator*(REAL a, const Matrix & A);
Matrix * operator+(const Matrix & A1, const Matrix & A2);
bool operator==(const Matrix & v1, const Matrix & v2);

#endif