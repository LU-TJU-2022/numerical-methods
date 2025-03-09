#include "particle.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include <iostream>
#include "fstream"
#include "config.h"
#include <vector>
#include "kernelFunction.h"
#include <ctime> 
using namespace std;

double calculate_T_by_Model3(struct LNode* nearPar, particle centerPar) {
	double* kernel;
	double seconder_X = 0;
	double seconder_X_1 = 0;
	double seconder_Y_1 = 0;
	double seconder_Y = 0;
	double dis;
	double res;
	while (nearPar != NULL) {
		dis = sqrt(pow((*nearPar).data.x - centerPar.x, 2) + pow((*nearPar).data.y - centerPar.y, 2));
		if (dis <= 2.0 * h) {
			if (dis != 0) {
				kernel = kernelFunction((*nearPar).data.x, (*nearPar).data.y, centerPar, 1);
				if ((centerPar.x - (*nearPar).data.x) != 0) {
					seconder_X += 2*(*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.T - centerPar.T) * ((*nearPar).data.x - centerPar.x) / pow(dis, 2) * kernel[0];
					seconder_X_1 += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.x - centerPar.x)  * kernel[0];

				}
				if ((centerPar.y - (*nearPar).data.y) != 0) {
					seconder_Y += 2*(*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.T - centerPar.T) * ((*nearPar).data.y - centerPar.y) / pow(dis, 2) * kernel[1];
					seconder_Y_1 += (*nearPar).data.mass / (*nearPar).data.density * ((*nearPar).data.y - centerPar.y) * kernel[1];
				}
			}
		}
		nearPar = (*nearPar).next;
	}
	res = centerPar.T + heatConductivity * (seconder_X/ seconder_X_1 + seconder_Y/ seconder_Y_1) * timeStep;
	return res;
}