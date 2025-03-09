#include "particle.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include <iostream>
#include "config.h"
using namespace std;

double* kernelFunction(double near_x, double near_y, particle centerPar, int order) {
	static double valueW[2];
	double dis = sqrt(pow(near_x - centerPar.x, 2) + pow(near_y - centerPar.y, 2));
	double x = near_x - centerPar.x;
	double y = near_y - centerPar.y;
	double q = dis / h;
	double factor;
	if (order == 0) {
		factor = 1.0f / h;
	}
	else if (order == 1) {
		factor = 15.0f / (7.0f * PI * h * h);
	}
	else if (order == 2) {
		factor = 3.0f / (2.0f * PI * h * h * h);
	}
	if (order == 0)
	{
		if (q >= 0 && q < 1)
		{
			valueW[0] = factor * (2.0 / 3.0 - q * q + 0.5 * q * q * q);
			valueW[1] = factor * (2.0 / 3.0 - q * q + 0.5 * q * q * q);
		}
		else if (q >= 1 && q < 2)
		{
			valueW[0] = factor * ((2.0 - q) * (2.0 - q) * (2.0 - q) / 6.0);
			valueW[1] = factor * ((2.0 - q) * (2.0 - q) * (2.0 - q) / 6.0);
		}
		else if (q >= 2)
		{
			valueW[0] = valueW[1] = 0.0;
		}
	}
	/*******一阶导数*******/
	if (order == 1)
	{
		if (q > 0 && q < 1)
		{
			valueW[0] = -factor * (1.5 * q * q - 2.0 * q) * (x / dis / h);
			valueW[1] = -factor * (1.5 * q * q - 2.0 * q) * (y / dis / h);
		}
		else if (q >= 1 && q < 2)
		{
			valueW[0] = -factor * (-0.5 * (2.0 - q) * (2.0 - q)) * (x / dis / h);
			valueW[1] = -factor * (-0.5 * (2.0 - q) * (2.0 - q)) * (y / dis / h);
		}
		else if (q >= 2 || q == 0)
		{
			valueW[0] = valueW[1] = 0.0;
		}
	}
	/*******二阶导数*******/
	if (order == 2)
	{
		if (q > 0 && q < 1)
		{
			valueW[0] = 1.0 / h * factor * (3.0 * q - 2.0) * 1.0 / h;//简化形式
			valueW[1] = 1.0 / h * factor * (3.0 * q - 2.0) * 1.0 / h;
		}
		else if (q >= 1 && q < 2)
		{
			valueW[0] = 1.0 / h * factor * (2.0 - q) * 1.0 / h;
			valueW[1] = 1.0 / h * factor * (2.0 - q) * 1.0 / h;
		}
		else if (q >= 2 || q == 0)
		{
			valueW[0] = valueW[1] = 0.0;
		}
	}
	return valueW;
}