#include <cmath>

#include "real.h"
#include "differences.h"

REAL d2f(const REAL h, const REAL f_r, const REAL f_m, const REAL f_l)
{
  return (f_r - 2*f_m + f_l)/h/h;
}

REAL df2(const REAL h, const REAL alpha, const REAL f_r, const REAL f_m, const REAL f_l)
{
  return ( std::pow((f_m+f_r)/2,2) - std::pow((f_m+f_l)/2,2) )/h
    +alpha/h * ( std::fabs(f_m+f_r)*(f_m-f_r)/4 - std::fabs(f_m+f_l)*(f_l-f_m)/4 );
} 

REAL dfg(const REAL h, const REAL alpha, const REAL f_r, const REAL f_m, const REAL f_l,
	 const REAL g_ll, const REAL g_lm, const REAL g_ul, const REAL g_um)
{
  return ( (g_lm+g_um)*(f_m+f_r)/4 - (g_ll+g_ul)*(f_l+f_m)/4 )/h
    + alpha/h * ( std::fabs(g_lm+g_um)*(f_m-f_r)/4 - std::fabs(g_ll-g_ul)*(f_l-f_m) );
}
