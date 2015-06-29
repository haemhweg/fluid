#ifndef __TRACING_CPP__
#define __TRACING_CPP__

#include <vector>
#include <utility>
#include <iostream>

#include "tracing.h"
#include "velocity.h"

Tracer::Tracer(const Config::tracing tracing_, Velocity * velocity_) : velocity(velocity_) , tracingConfig(tracing_) , particleLines() {

	REAL x1 = tracingConfig.x1;
	REAL x2 = tracingConfig.x2;
	REAL y1 = tracingConfig.y1;
	REAL y2 = tracingConfig.y2;
	unsigned int N = tracingConfig.N;

	for(int i = 0; i < N; ++i) {
		ParticleLine pl(velocity, i * (x2 - x1) / (N - 1) + x1 , i * (y2 - y1) / (N - 1) + y1  );
		particleLines.push_back(pl);
	}

}

void Tracer::advance(REAL delt) {

	for(int i = 0; i < particleLines.size(); ++i) {
		particleLines[i].advance(delt);
	}

}

void Tracer::output() {

	ParticleLine first = particleLines[0];
	first.output();

}

ParticleLine::ParticleLine(Velocity * velocity_, REAL x, REAL y) : velocity(velocity_) , particles() {

	Particle p(x, y);
	particles.push_back(p);

}

void ParticleLine::advance(REAL delt) {

	unsigned int i, j;

	Particle currentParticle = particles[particles.size() - 1];

	particles.push_back(velocity->advanceParticle(currentParticle, delt));

}

void ParticleLine::output() {

	for(int i = 0; i < particles.size(); ++i) {
		std::cout << particles[i].x << " " << particles[i].y << std::endl;
	}

}

#endif