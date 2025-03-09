#include "particle.h"
#include <iostream>
#include "config.h"
particle::particle()
{
	PARTICLETYPE type;
	x = 0;
	y = 0;
	u = 0;
	v = 0;
	T = 0;
	heatFlux_x = 0;
	heatFlux_y = 0;
	heatFlux_x_left = 0;
	heatFlux_x_right = 0;
	heatFlux_y_left = 0;
	heatFlux_y_right = 0;
	mass = initialDistance * initialDistance;
	density = 1;
	row = 0;
	line = 0;
}

void particle::changeCoordinate(double newX, double newY) {
	x = newX;
	y = newY;
}

void particle::changeType(PARTICLETYPE newType) {
	type = newType;
}

void particle::changeGridNumber(int newRow, int newLine) {
	row = newRow;
	line = newLine;
}

void particle::changeHeatFlux(double newFlux_x, double newFlux_y) {
	heatFlux_x = newFlux_x;
	heatFlux_y = newFlux_y;
}

void particle::changeT(double newT) {
	T = newT;
}

void particle::changeVelocity(double newU, double newV) {
	u = newU;
	v = newV;
}