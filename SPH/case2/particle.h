#pragma once
/*-------.h…˘√˜£¨.cpp µœ÷----------*/
enum PARTICLETYPE
{
	itinlet = 0,
	itoutlet = 1,
	itfluid = 2,
	itwall = 3,
};

class particle
{
public:
	PARTICLETYPE type;
	double x;
	double y;
	double u;
	double v;
	double T;
	double heatFlux_x;
	double heatFlux_y;
	double heatFlux_x_left;
	double heatFlux_x_right;
	double heatFlux_y_left;
	double heatFlux_y_right;
	double mass;
	double density;
	int row;
	int line;
public:
	particle();
	void changeCoordinate(double newX, double newY);
	void changeType(PARTICLETYPE newType);
	void changeGridNumber(int newRow, int newLine);
	void changeHeatFlux(double newFlux_x, double newFlux_y);
	void changeT(double newT);
	void changeVelocity(double newU, double newV);
};

typedef struct LNode {
	particle data;
	LNode* next;
}*p;

