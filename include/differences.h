#ifndef DIFFERENCES_H
#define DIFFERENCES_H

#include "real.h"

class Field2D;

/**
 *  @brief Computes the second derivative of a function f(x) using the discretization
 *         f''(x) = (f(x-h) - 2f(x) + f(x+h)) / h
 *  @param h Discretization size
 *  @param f_r,f_m,f_l The values f(x-h), f(x), f(x-h)
 */
REAL d2f(const REAL h, const REAL f_r, const REAL f_m, const REAL f_l);

/**
 *  @brief Computes the first derivative of f(x)^2 using the central difference and Donor-Cell scheme
 *  @param h Discretization size
 *  @param alpha If alpha=0 only the central difference scheme is choosen, if alpha=1 only the 
 *         Donor-Cell scheme is used, otherwise a weighted combination of both
 *  @param f_r,f_m,f_l The values f(x-h), f(x), f(x-h)
 */
REAL df2(const REAL h, const REAL alpha, const REAL f_r, const REAL f_m, const REAL f_l);

/**
 *  @brief Computes the first derivative of f(x,y)*g(x,y) in x-direction using the central difference and Donor-Cell scheme.
 *         For the y-direction swap f and g.
 *  @param h Discretization size
 *  @param alpha If alpha=0 only the central difference scheme is choosen, if alpha=1 only the 
 *         Donor-Cell scheme is used, otherwise a weighted combination of both
 *  @param f_r,f_m,f_l,g_ll,g_lm,g_ul,g_um The values f(x-h,y), f(x,y), f(x+h,y) and g(x-h,y),g(x,y),g(x-h,y+h),g(x,y+h)
 */
REAL dfg(const REAL h, const REAL alpha, const REAL f_r, const REAL f_m, const REAL f_l,
	 const REAL g_ll, const REAL g_lm, const REAL g_ul, const REAL g_um);



REAL d2fdx(const REAL h, const Field2D& f, const unsigned i, const unsigned j);
REAL d2fdy(const REAL h, const Field2D& f, const unsigned i, const unsigned j);

REAL df2dx(const REAL h, const REAL alpha, const Field2D& f, const unsigned i, const unsigned j);
REAL df2dy(const REAL h, const REAL alpha, const Field2D& f, const unsigned i, const unsigned j);

REAL dfgdx(const REAL h, const REAL alpha, const Field2D& f, const Field2D& g, 
	   const unsigned i, const unsigned j);
REAL dfgdy(const REAL h, const REAL alpha, const Field2D& f, const Field2D& g, 
	   const unsigned i, const unsigned j);
  


#endif
