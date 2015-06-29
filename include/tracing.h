#ifndef __TRACING_H__
#define __TRACING_H__

#include <vector>
#include <utility>

#include "config.h"

class Velocity;

class Particle {
public:
	Particle(REAL x_, REAL y_) : x(x_), y(y_) {};
	REAL x, y;
};

class ParticleLine {
public:
	// constructor with initial values
	ParticleLine(Velocity * velocity_, REAL, REAL);
	void advance(REAL delt);
	void output();

private:
	Velocity * velocity;
	std::vector<Particle> particles;

};

class Tracer {
public:
	// pass configs and pointer to velocity
	Tracer(const Config::tracing tracing_, Velocity * velocity_);
	void advance(REAL delt);
	void output();

private:
	Config::tracing tracingConfig;
	Velocity * velocity;
	std::vector<ParticleLine> particleLines;

};

#endif