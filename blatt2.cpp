#include <iostream>
#include <algorithm>

#include "real.h"
#include "matrix.h"
#include "config.h"

#include "differences.h"

void initUVP(Config::geo geoConfig, Config::constants constantsConfig, Matrix & U, Matrix & V, Matrix & P)
{
	// U = Matrix(geoConfig.imax + 2, geoConfig.jmax + 2, constantsConfig.UI);
	// V = Matrix(geoConfig.imax + 2, geoConfig.jmax + 2, constantsConfig.VI);
	// P = Matrix(geoConfig.imax + 2, geoConfig.jmax + 2, constantsConfig.PI);
}

REAL computeDelta(Config::time timeConfig, Config::geo geoConfig, Config::constants constantsConfig, Matrix * U, Matrix * V) {
	
	// REAL vals[3];

	// vals[0] = constantsConfig.Re / 2 / (1 / geoConfig.delx / geoConfig.delx + 1 / geoConfig.dely / geoConfig.dely);
	// vals[1] = geoConfig.delx / U->getMax();
	// vals[2] = geoConfig.dely / V->getMax();

	// return timeConfig.tau * *std::min_element(vals, vals + 3);

}

void setBoundaryConditions(Matrix * U, Matrix * V)
{
	/**
	 * Iteration laut (18)
	 */
}

void computeFG(Config::geo geoConfig, Config::time timeConfig, Config::constants constantsConfig, Matrix const & U, Matrix const & V, Matrix & F, Matrix & G)
{
	// REAL delt = timeConfig.delt;
	// REAL Re = constantsConfig.Re;

	// for (size_t i = 1; i < geoConfig.imax; ++i) {
	// 	for (size_t j = 1; j < geoConfig.jmax + 1; ++j)
	// 	{
	// 		F.at(i, j) = U.at(i, j) + delt * (1 / Re * (d2Udx2(geoConfig, constantsConfig, U, i, j) + (d2Udy2(geoConfig, constantsConfig, U, i, j))) - dU2dx(geoConfig, constantsConfig, U, i, j) - dUVdx(geoConfig, constantsConfig, U, V, i, j) + constantsConfig.GX);
	// 	}
	// }

}

int main()
{

	/**
	 * Falls man Matrizen aus diesem Scope in einer Subroutine erstellen will, so müssen hier nur Pointer definiert sein.
	 */
	Matrix U();
	Matrix V();
	Matrix P();
	Matrix RHS();
	Matrix F();
	Matrix G();

	/* Lese Konfigurationsdatei ein. */
	Config conf("config");

}
