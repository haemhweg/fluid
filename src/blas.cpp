#include <algorithm>
#include <cassert>

#include "vector.h"
#include "matrix.h"
#include "real.h"

Vector axpy(REAL a, const Vector & x, const Vector & y)
{
	return (a * x) + y;
}

REAL dot(const Vector & x, const Vector & y)
{
	return x * y;
}

REAL nrm2(const Vector & x)
{
	return x.nrm2();
}

void copy(const Vector & source, Vector & dest)
{
  std::copy(std::begin(source), std::end(source), std::begin(dest));
}

Vector scal(REAL a, Vector & x)
{
	return a * x;
}

Vector gemv(REAL a, const Matrix & A, const Vector & x, REAL b, const Vector & y)
{

  assert(A.getCols() != x.getSize() || A.getRows() != y.getSize());
  
  return a * (A * x) + (b * y);

}

Matrix scal2Dfield(REAL a, Matrix & X)
{
	return a * X;
}

Matrix axpy2Dfield(REAL a, const Matrix & X, const Matrix & Y)
{
	return (a * X) + Y;
}
