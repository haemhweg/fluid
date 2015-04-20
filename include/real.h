#ifndef __REAL_H__
#define __REAL_H__

#include <cmath>

typedef double REAL;

const REAL PI = atan(REAL(1.))*4;

bool REAL_equal(const REAL a, const REAL b);

#endif
