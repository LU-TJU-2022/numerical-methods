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

double get_boundary_T(struct LNode* nearPar, particle centerPar) {
	double* kernel;
	double t_x = 0;
	double T_x_1 = 0;
	double t_y = 0;
	double T_y_1 = 0;
	double T = 0;
	double seconder_Y_1 = 0;
	double seconder_Y = 0;
	double dis;
	double res;
	while (nearPar != NULL) {
		dis = sqrt(pow((*nearPar).data.x - centerPar.x, 2) + pow((*nearPar).data.y - centerPar.y, 2));
		if (dis <= 2.0 * h) {
			if (dis != 0) {
				kernel = kernelFunction((*nearPar).data.x, (*nearPar).data.y, centerPar, 0);
				if ((centerPar.x - (*nearPar).data.x) != 0) {
					t_x += (*nearPar).data.T * kernel[0];
					T_x_1 += kernel[0];
				}
				/*else if ((centerPar.y - (*nearPar).data.y) != 0)
				{
					t_y += (*nearPar).data.T * kernel[1];
					T_y_1 += kernel[1];
				}*/
			}
			/*if (abs(nearPar.data.y - 0) < 1e-8)
			{
				T = nearPar.data.T;
			}
			else if(abs(nearPar.data.y - 1) < 1e-8)
			{
				T = nearPar.data.T;
			}*/
		}
		nearPar = (*nearPar).next;
	}
	if (T_x_1 != 0 && T_y_1 != 0)
	{
		res = t_x / T_x_1 + t_y / T_y_1;
	}
	else if (T_x_1 != 0 && T_y_1 == 0)
	{
		res = t_x / T_x_1;
	}
	else if (T_x_1 == 0 && T_y_1 != 0)
	{
		res = t_y / T_y_1;
	}
	else
	{
		res = 0;
	}
	return res;
}
