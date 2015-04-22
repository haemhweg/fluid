#include <iostream>
#include <algorithm>

#include "real.h"
#include "matrix.h"
#include "config.h"

void initUVP(Config::geo geoConfig, Config::constants constantsConfig, Matrix * U, Matrix * V, Matrix * P)
{
	/**
	 * Da U beschreibt waagerechte Werte, braucht man eigentlich die Werte 0 und jmax + 1 nicht (analog bei V).
	 * Wir definieren aber diese Werte, um die Notation konsistent zu halten.
	 */
	U = new Matrix(geoConfig.imax + 2, geoConfig.jmax + 2, constantsConfig.UI);
	V = new Matrix(geoConfig.imax + 2, geoConfig.jmax + 2, constantsConfig.VI);
	P = new Matrix(geoConfig.imax + 2, geoConfig.jmax + 2, constantsConfig.PI);
}

REAL computeDelta(Config::time timeConfig, Config::geo geoConfig, Config::constants constantsConfig, Matrix * U, Matrix * V) {
	
	REAL vals[3];

	vals[0] = constantsConfig.Re / 2 / (1 / geoConfig.delx / geoConfig.delx + 1 / geoConfig.dely / geoConfig.dely);
	vals[1] = geoConfig.delx / U->getMax();
	vals[2] = geoConfig.dely / V->getMax();

	return timeConfig.tau * *std::min_element(vals, vals + 3);

}

void setBoundaryConditions(Matrix * U, Matrix * V)
{
	/**
	 * Iteration laut (18)
	 */
}

int main()
{

	/**
	 * Falls man Matrizen aus diesem Scope in einer Subroutine erstellen will, so m�ssen hier nur Pointer definiert sein.
	 */
	Matrix * U;
	Matrix * V;
	Matrix * P;
	Matrix * RHS;
	Matrix * F;
	Matrix * G;

	/* Lese Konfigurationsdatei ein. */
	Config conf("config.xml");

}