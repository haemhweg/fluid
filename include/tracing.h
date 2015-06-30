#ifndef __TRACING_H__
#define __TRACING_H__

#include <deque>
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
	void output();
	void inject();
	void advance(REAL delt, TracingType tr);

private:
	Velocity * velocity;
	// save the values of the first particle in case we need to create inject a new one for streaklines
	REAL x, y;
	std::deque<Particle> particles;

};

class Tracer {
public:
	// pass configs and pointer to velocity
	Tracer(const Config::tracing tracing_, Velocity * velocity_);
	void advance(REAL delt);
	void inject();
	void output();

private:
	Config::tracing tracingConfig;
	Velocity * velocity;
	std::deque<ParticleLine> particleLines;
	REAL delt_write_residual, delt_inject_residual;

};

#endif