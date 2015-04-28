#include <cmath>

#include "differences.h"
#include "real.h"
#include "matrix.h"

REAL d2f(const REAL h, const REAL f_l, const REAL f_m, const REAL f_r)
{
  return (f_l - 2*f_m + f_r)/h/h;
}

REAL df2(const REAL h, const REAL alpha, const REAL f_l, const REAL f_m, const REAL f_r)
{
  return ( std::pow((f_m+f_r)/2,2) - std::pow((f_m+f_l)/2,2) )/h
    +alpha/h * ( std::fabs(f_m+f_r)*(f_m-f_r)/4 - std::fabs(f_m+f_l)*(f_l-f_m)/4 );
} 

REAL dfg(const REAL h, const REAL alpha, const REAL f_l, const REAL f_m, const REAL f_r,
	 const REAL g_ll, const REAL g_lm, const REAL g_ul, const REAL g_um)
{
  return ( (g_lm+g_um)*(f_m+f_r)/4 - (g_ll+g_ul)*(f_l+f_m)/4 )/h
    + alpha/h * ( std::fabs(g_lm+g_um)*(f_m-f_r)/4 - std::fabs(g_ll-g_ul)*(f_l-f_m) );
}


REAL d2fdx(const REAL h, const Matrix& f, const unsigned i, const unsigned j)
{
  return d2f(h, f.at(i-1,j), f.at(i,j), f.at(i+1,j));
}
REAL d2fdy(const REAL h, const Matrix& f, const unsigned i, const unsigned j)
{
  return d2f(h, f.at(i,j-1), f.at(i,j), f.at(i,j+1));
}

REAL df2dx(const REAL h, const REAL alpha, const Matrix& f, const unsigned i, const unsigned j)
{
  return df2(h, alpha, f.at(i-1,j), f.at(i,j), f.at(i+1,j));
}
REAL df2dy(const REAL h, const REAL alpha, const Matrix& f, const unsigned i, const unsigned j)
{
  return df2(h, alpha, f.at(i,j-1), f.at(i,j), f.at(i,j+1));
}

REAL dfgdx(const REAL h, const REAL alpha, const Matrix& f, const Matrix& g, 
	   const unsigned i, const unsigned j)
{
  return dfg(h, alpha, g.at(i-1,j), g.at(i,j), g.at(i+1,j), f.at(i-1,j), f.at(i,j),
	     f.at(i-1,j+1), f.at(i,j+1));
}
REAL dfgdy(const REAL h, const REAL alpha, const Matrix& f, const Matrix& g, 
	   const unsigned i, const unsigned j)
{
  return dfg(h, alpha, f.at(i,j-1), f.at(i,j), f.at(i,j+1), g.at(i,j-1), g.at(i,j),
	     g.at(i+1,j-1), g.at(i+1,j));
}
  
