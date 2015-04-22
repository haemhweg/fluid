#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "real.h"

class Config
{
	/**
	 * Nur ein Vorschlag. Eigentlich muss man hier den Config irgendwie enkapsulieren.
	 */
public:
	struct geo
	{
		REAL xlength;
		REAL ylength;
		int imax;
		int jmax;
		REAL delx;
		REAL dely;
	};
	struct time
	{
		REAL t;
		REAL t_end;
		REAL delt;
		REAL tau;
		REAL del_vec;
	};
	struct pressure
	{
		int itermax;
		int it;
		REAL res;
		REAL eps;
		REAL omg;
		REAL alpha;
	};
	struct constants
	{
		REAL Re;
		REAL GX, GY;
		REAL UI, VI, PI;
	};
	
	geo _geo;
	time _time;
	pressure _pressure;
	constants _constants;

	Config(char * filename);
};

#endif