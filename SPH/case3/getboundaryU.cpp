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

double get_boundary_u(struct LNode* nearPar, particle centerPar) {
	double* kernel;
	double u = 0;
	double u_1 = 0;
	double u_y = 0;
	double u_1_y = 0;
	double U = 0;
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
					u += (*nearPar).data.u * kernel[0];
					u_1 += kernel[0];
				}
			}
			/*if (abs(nearPar.data.y - 0) < 1e-8)
			{
				U = nearPar.data.u;
			}
			else if(abs(nearPar.data.y - 1) < 1e-8)
			{
				U = nearPar.data.u;
			}*/
		}
		nearPar = (*nearPar).next;
	}
	if (u_1 != 0)
	{
		if (centerPar.y <0)
		{
			res = - u / u_1;
		}
		else 
		{
			res = 2 - u / u_1;
		}
		/*else
		{
			res = 2 * 0.05 - u / u_1;
		}*/
	}
	else if(u_1 == 0)
	{
		if (centerPar.y < 0)
		{
			res = 0.0f;
		}
		else
		{
			res = 2.0f;
		}
	}
	return res;
}