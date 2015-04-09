/**
 * Use given function names due to standardization reasons. The functions below might use methods of the classes.
 */

#include "vector.h"
#include "matrix.h"
#include "real.h"

Vector * axpy(REAL a, const Vector & x, const Vector & y);
REAL dot(const Vector & x, const Vector & y);
REAL nrm2(const Vector & x);
void copy(const Vector & source, Vector & dest);
Vector * scal(REAL a, Vector & x);
Vector * gemv(REAL a, const Matrix & A, const Vector & x, REAL b, const Vector & y);

Matrix * scal2Dfield(REAL a, Matrix & X);
Matrix * axpy2Dfield(REAL a, const Matrix & X, const Matrix & Y);