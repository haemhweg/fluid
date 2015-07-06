#ifndef __TRACING_CPP__
#define __TRACING_CPP__

#include <deque>
#include <utility>
#include <sstream>
#include <iostream>
#include <fstream>

#include "tracing.h"
#include "velocity.h"
#include "config.h"

Tracer::Tracer(const Config::tracing tracing_, const Config::geo geo_, Velocity * velocity_) : velocity(velocity_) , tracingConfig(tracing_) , geoConfig(geo_) , particleLines() {

	delt_write_residual = tracingConfig.delt_write;
	delt_inject_residual = tracingConfig.delt_inject;

	REAL x1 = tracingConfig.x1;
	REAL x2 = tracingConfig.x2;
	REAL y1 = tracingConfig.y1;
	REAL y2 = tracingConfig.y2;
	unsigned int N = tracingConfig.N;

	for(int i = 0; i < N; ++i) {
		ParticleLine pl(velocity, i * (x2 - x1) / (N - 1) + x1 , i * (y2 - y1) / (N - 1) + y1  );
		particleLines.push_front(pl);
	}

}

void Tracer::advance(REAL delt) {

	// inject only if streaklines, otherwise ignore parameter
	delt_inject_residual -= delt;	
	if(delt_inject_residual < 0 && tracingConfig.tr == STREAKLINES) {
		delt_inject_residual = tracingConfig.delt_inject;
		inject();
	}

	delt_write_residual -= delt;
	if(delt_write_residual < 0) {
		delt_write_residual = tracingConfig.delt_write;
		output(writeStep++);
	}

	for(size_t i = 0; i < particleLines.size(); ++i) {
		particleLines[i].advance(delt, tracingConfig.tr);
	}

}

void Tracer::output(size_t writeStep) {

	std::ostringstream filename;

	filename << "Streaklines" << writeStep << ".vtk";

	std::ofstream fs(filename.str());

	fs << std::scientific;

	fs << "# vtk DataFile Version 3.0\n"
     << "Streaklines\n"
     << "ASCII\n"
     << "DATASET POLYDATA\n";

    size_t nodesSize = 0;

    for(size_t i = 0; i < particleLines.size(); i++) {
    	nodesSize += particleLines[i].getSize();
    }

    fs << "POINTS " << nodesSize << " double\n";

    for(size_t i = 0; i < particleLines.size(); i++) {
    	particleLines[i].output(fs);
    }

    fs << "LINES " << particleLines.size() << " " << (nodesSize + particleLines.size()) << std::endl;

    size_t k = 0;
    for(size_t i = 0; i < particleLines.size(); i++) {
    	fs << particleLines[i].getSize() << " ";
    	for(size_t j = 0; j < particleLines[i].getSize(); j++) {
    		fs << (k + j) << " ";
    	}
    	k += particleLines[i].getSize();
    	fs << std::endl;
    }

}

void Tracer::inject() {

	for(int i = 0; i < particleLines.size(); ++i) {
		particleLines[i].inject();
	}

}

ParticleLine::ParticleLine(Velocity * velocity_, REAL x_, REAL y_) : velocity(velocity_) , particles(), x(x_), y(y_) {

	inject();

}

void ParticleLine::inject() {

	Particle p(x, y);
	particles.push_front(p);	

}

void ParticleLine::advance(REAL delt, TracingType tr) {

	unsigned int i, j;

	// so basically if we have pathlines we advance only last particle and keep the rest, and if we have streaklines we advance all of them
	if(tr == PATHLINES) {
		Particle * currentParticle =& particles[particles.size() - 1];
		std::pair<REAL, REAL> vals = velocity->advanceParticle(currentParticle->x, currentParticle->y, delt);
		particles.push_back(Particle(vals.first, vals.second));
	} else if (tr == STREAKLINES) {
		for(int i = 0; i < particles.size(); ++i) {
			Particle * currentParticle =& particles[i];
			try {
				std::pair<REAL, REAL> vals = velocity->advanceParticle(currentParticle->x, currentParticle->y, delt);
				currentParticle->x = vals.first;
				currentParticle->y = vals.second;
			} catch(int n) {
				particles.resize(i);
				return;
			}
		}
	}

}

void ParticleLine::output(std::ofstream & fs) {

	for(size_t i = 0; i < particles.size(); ++i) {
		fs << particles[i].x << " " << particles[i].y << " 0" << std::endl;
	}

}

size_t ParticleLine::getSize() {

	return particles.size();

}

#endif