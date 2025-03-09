#include "particle.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include <iostream>
#include "fstream"
#include "config.h"
#include <vector>
#include "kernelFunction.h"
using namespace std;

vector<double> calculate_heatFlux_by_Model2(struct LNode* nearPar, particle centerPar) {
	double* kernel;
	double heatFlux_x = 0;
	double heatFlux_y = 0;
	double heatFlux_x_1 = 0;
	double heatFlux_y_1 = 0;
	vector<double> res(2, 0);
	while (nearPar != NULL) {
		kernel = kernelFunction((*nearPar).data.x, (*nearPar).data.y, centerPar, 1);
		heatFlux_x += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.T - centerPar.T) * kernel[0];
		heatFlux_x_1 += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.x - centerPar.x) * kernel[0];
		heatFlux_y += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.T - centerPar.T) * kernel[1];
		heatFlux_y_1 += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.y - centerPar.y) * kernel[1];
		nearPar = (*nearPar).next;
	}
	res[0] = heatFlux_x / heatFlux_x_1;
	res[1] = heatFlux_y / heatFlux_y_1;
	return res;
}

double calculate_T_by_Model2(struct LNode* nearPar, particle centerPar) {
	double* kernel;
	double heatFlux_x = 0;
	double heatFlux_y = 0;
	double heatFlux_x_1 = 0;
	double heatFlux_y_1 = 0;
	double res;
	while (nearPar != NULL) {
		kernel = kernelFunction((*nearPar).data.x, (*nearPar).data.y, centerPar, 1);
		heatFlux_x += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.heatFlux_x - centerPar.heatFlux_x) * kernel[0];
		heatFlux_y += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.heatFlux_y - centerPar.heatFlux_y) * kernel[1];
		heatFlux_x_1 += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.x - centerPar.x) * kernel[0];
		heatFlux_y_1 += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.y - centerPar.y) * kernel[1];
		nearPar = (*nearPar).next;
	}
	res = centerPar.T + heatConductivity * (heatFlux_x / heatFlux_x_1 + heatFlux_y / heatFlux_y_1) * timeStep;
	return res;
}

vector<double> calculate_heatFlux_by_Model4(struct LNode* nearPar, particle centerPar) {
	double* kernel;
	double heatFlux_x_left = 0;
	double heatFlux_y_left = 0;
	double heatFlux_x_right = 0;
	double heatFlux_y_right = 0;
	double heatFlux_x_left_1 = powf(10, -10);
	double heatFlux_y_left_1 = powf(10, -10);
	double heatFlux_x_right_1 = powf(10, -10);
	double heatFlux_y_right_1 = powf(10, -10);
	vector<double> res(4, 0);
	while (nearPar != NULL) {
		double a = (*nearPar).data.x - centerPar.x;
		double b = (*nearPar).data.y - centerPar.y;
		kernel = kernelFunction((*nearPar).data.x, (*nearPar).data.y, centerPar, 1);
		if ((a + b) >= 0 && (a - b) <= 0) {
			heatFlux_x_left += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.T - centerPar.T) * kernel[1];
			heatFlux_x_left_1 += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.y - centerPar.y) * kernel[1];
		}
		if ((a + b) <= 0 && (a - b) >= 0) {
			heatFlux_x_right += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.T - centerPar.T) * kernel[1];
			heatFlux_x_right_1 += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.y - centerPar.y) * kernel[1];
		}
		if ((a + b) >= 0 && (a - b) >= 0) {
			heatFlux_y_left += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.T - centerPar.T) * kernel[0];
			heatFlux_y_left_1 += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.x - centerPar.x) * kernel[0];
		}
		if ((a + b) <= 0 && (a - b) <= 0) {
			heatFlux_y_right += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.T - centerPar.T) * kernel[0];
			heatFlux_y_right_1 += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.x - centerPar.x) * kernel[0];
		}
		nearPar = (*nearPar).next;
	}
	res[0] = heatFlux_x_left / heatFlux_x_left_1;
	res[1] = heatFlux_x_right / heatFlux_x_right_1;
	res[2] = heatFlux_y_left / heatFlux_y_left_1;
	res[3] = heatFlux_y_right / heatFlux_y_right_1;
	return res;
}

double calculate_T_by_Model4(struct LNode* nearPar, particle centerPar) {
	double* kernel;
	double T_x_left = 0;
	double T_x_right = 0;
	double T_y_left = 0;
	double T_y_right = 0;
	double T_x_left_1 = powf(10, -10);
	double T_x_right_1 = powf(10, -10);
	double T_y_left_1 = powf(10, -10);
	double T_y_right_1 = powf(10, -10);
	double res;
	while (nearPar != NULL) {
		double a = (*nearPar).data.x - centerPar.x;
		double b = (*nearPar).data.y - centerPar.y;
		kernel = kernelFunction((*nearPar).data.x, (*nearPar).data.y, centerPar, 1);
		if ((a + b) >= 0 && (a - b) <= 0) {
			T_x_left += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.heatFlux_x_right - centerPar.heatFlux_x_right) * kernel[1];
			T_x_left_1 += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.y - centerPar.y) * kernel[1];
		}
		if ((a + b) <= 0 && (a - b) >= 0) {
			T_x_right += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.heatFlux_x_left - centerPar.heatFlux_x_left) * kernel[1];
			T_x_right_1 += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.y - centerPar.y) * kernel[1];
		}
		if ((a + b) >= 0 && (a - b) >= 0) {
			T_y_left += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.heatFlux_y_right - centerPar.heatFlux_y_right) * kernel[0];
			T_y_left_1 += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.x - centerPar.x) * kernel[0];
		}
		if ((a + b) <= 0 && (a - b) <= 0) {
			T_y_right += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.heatFlux_y_left - centerPar.heatFlux_y_left) * kernel[0];
			T_y_right_1 += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.x - centerPar.x) * kernel[0];
		}
		nearPar = (*nearPar).next;
	}
	res = centerPar.T + heatConductivity * (T_x_left / T_x_left_1 + T_x_right / T_x_right_1 + T_y_left / T_y_left_1 + T_y_right / T_y_right_1) * timeStep / 2.0f;
	return res;
}

