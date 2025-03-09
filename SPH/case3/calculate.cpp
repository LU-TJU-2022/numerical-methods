#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include <iostream>
#include "fstream"
#include "calculate.h"
#include "config.h"
#include "particle.h"
#include <ctime>
#include "particleGeneration.h"
#include "particleSearch.h"
#include "updatePosition.h"
#include "gridReset.h"
#include "model2.h"
#include "model3.h"
#include "model4.h"
#include "getboundaryT.h"
#include "getboundaryU.h"

using namespace std;

void saveData(struct LNode** grid, int timeStep) {/*--保存计算结果--*/
	cout.setf(ios::fixed);
	cout.precision(10);
	char szName[100] = {'\0'};
	sprintf_s(szName, "model%d时间步%d结果10000粒子case3.txt", model,timeStep);
	ofstream dataFile;
	dataFile.open(szName, ios::out);
	if (!dataFile) {
		cout << "error" << endl;
		exit(0);
	}
	int rowNum = grid[0][0].data.row;
	int lineNum = grid[0][0].data.line;
	struct LNode* particleNow;
	for (int i = 0; i < rowNum; i++) {
		for (int j = 0; j < lineNum; j++) {
			particleNow = grid[i][j].next;
			while (particleNow != NULL) {

					dataFile << (*particleNow).data.x << '\t' << (*particleNow).data.y << '\t' << (*particleNow).data.T << '\t' << (*particleNow).data.row << '\t' << (*particleNow).data.line << '\t' << (*particleNow).data.u << '\t' << (*particleNow).data.type << '\t' << (*particleNow).data.heatFlux_x_left << '\t' << (*particleNow).data.heatFlux_x_right << '\t' << (*particleNow).data.heatFlux_y_left << '\t' << (*particleNow).data.heatFlux_y_right << endl;
					particleNow = (*particleNow).next;

			}
		}
	}
}

void calculate() {
	
	time_t start, stop;
	start = time(NULL);
	struct LNode** grid;
	struct LNode* particleNow;
	struct LNode* particlePrev;
	struct LNode* particleNear;
	struct LNode* nearPar;
	struct LNode* nearPar_next;
	grid = Generation(coordination, initialDistance, initialDistancey, h, hy);
	vector<double> newHeatFlux;
	double newT;
	double newU;
	int rowNum = grid[0][0].data.row;
	int lineNum = grid[0][0].data.line;

	struct LNode** gridCopy;
	gridCopy = new LNode * [rowNum];
	for (int n = 0; n < rowNum; n++) {
		gridCopy[n] = new LNode[lineNum];
		for (int m = 0; m < lineNum; m++) {
			gridCopy[n][m].next = NULL;
		}
	}
	saveData(grid, 0);
	for (int i = 0; i <= 1200; i++) {
		cout <<"第"<< i<<"时间步" << endl;
		/*for (int j = 0; j < rowNum; j++) {
			for (int k = 0; k < lineNum; k++) {
				particleNow = grid[j][k].next;
				while (particleNow != NULL) {
					if ((*particleNow).data.type == itwall) {
						particleNear = particleSearch(grid, (*particleNow).data, 2 * h);

						if (model == 2) {
							newHeatFlux = calculate_heatFlux_by_Model2(particleNear, (*particleNow).data);
							(*particleNow).data.heatFlux_x = newHeatFlux[0];
							(*particleNow).data.heatFlux_y = newHeatFlux[1];
						}
						if ((*particleNow).data.type == itwall) {
							if (abs((*particleNow).data.y - 0) < 1e-8)
							{
								newT = 0.0f;
								newU = 0.0f;
								(*particleNow).data.changeVelocity(newU, 0.0f);

							}
							else if (abs((*particleNow).data.y - 1) < 1e-8)
							{
								newT = 1.0f;
								newU = 1.0f;
								(*particleNow).data.changeVelocity(newU, 0.0f);
							}
							else
							{
								if ((*particleNow).data.y < 0)
								{
									newT = 0.0f;
								}
								else
								{
									newT = 1.0f;
								}
								//newT = get_boundary_T(particleNear, (*particleNow).data);
								newU = get_boundary_u(particleNear, (*particleNow).data);
								(*particleNow).data.changeVelocity(newU, 0.0f);
							}
							(*particleNow).data.heatFlux_x = newT;
							(*particleNow).data.changeT((*particleNow).data.heatFlux_x);
							
						}
						else {
							newHeatFlux = calculate_heatFlux_by_Model4_1(particleNear, (*particleNow).data);
							(*particleNow).data.heatFlux_x_left = newHeatFlux[0];
							(*particleNow).data.heatFlux_x_right = newHeatFlux[1];
							(*particleNow).data.heatFlux_y_left = newHeatFlux[2];
							(*particleNow).data.heatFlux_y_right = newHeatFlux[3];
						}
						nearPar = particleNear;
						while (nearPar != NULL) {
							nearPar_next = (*nearPar).next;
							delete(nearPar);
							nearPar = nearPar_next;
						}

					}
					particleNow = (*particleNow).next;
				}
			}
		}*/
		for (int j = 0; j < rowNum; j++) {
			for (int k = 0; k < lineNum; k++) {
				particleNow = grid[j][k].next;
				while (particleNow != NULL) {
					if ((*particleNow).data.type == itfluid || (*particleNow).data.type == itboundary) {
						particleNear = particleSearch(grid, (*particleNow).data, 2 * h);
						
						if (model == 2) {
							newHeatFlux = calculate_heatFlux_by_Model2(particleNear, (*particleNow).data);
							(*particleNow).data.heatFlux_x = newHeatFlux[0];
							(*particleNow).data.heatFlux_y = newHeatFlux[1];
						}
						if (model == 3) {

						}
						else if (model == 4)
						{
							newHeatFlux = calculate_heatFlux_by_Model4_1(particleNear, (*particleNow).data);
							(*particleNow).data.heatFlux_x_left = newHeatFlux[0];
							(*particleNow).data.heatFlux_x_right = newHeatFlux[1];
							(*particleNow).data.heatFlux_y_left = newHeatFlux[2];
							(*particleNow).data.heatFlux_y_right = newHeatFlux[3];
						}
						nearPar = particleNear;
						while (nearPar != NULL) {
							nearPar_next = (*nearPar).next;
							delete(nearPar);
							nearPar = nearPar_next;
						}

					}
					particleNow = (*particleNow).next;
				}
			}
		}
		for (int j = 0; j < rowNum; j++) {
			for (int k = 0; k < lineNum; k++) {
				particleNow = grid[j][k].next;
				while (particleNow != NULL) {
					if ((*particleNow).data.type == itfluid) {
						particleNear = particleSearch(grid, (*particleNow).data, 2 * h);
						if (model == 2) {
							newT = calculate_T_by_Model2(particleNear, (*particleNow).data);
							(*particleNow).data.T = newT;
						}
						else if (model == 3) {
							newT = calculate_T_by_Model3(particleNear, (*particleNow).data);
							(*particleNow).data.changeT(newT);
						}
						else if (model == 4)
						{
							newT = calculate_T_by_Model4_1(particleNear, (*particleNow).data);
							(*particleNow).data.T = newT;
						}						
						nearPar = particleNear;
						while (nearPar != NULL) {
							nearPar_next = (*nearPar).next;
							delete(nearPar);
							nearPar = nearPar_next;
						}

					}
					particleNow = (*particleNow).next;
				}
			}
		}
		updatePosition(grid, gridCopy);
		gridReset(grid, gridCopy);
		
			
		
	}
	
	for (int j = 0; j < rowNum; j++) {
		for (int k = 0; k < lineNum; k++) {
			particleNow = grid[j][k].next;
			particlePrev = &(grid[j][k]);
			while (particleNow != NULL) {
				if ((*particleNow).data.x <= -6400 || (*particleNow).data.x >= 6400 || (*particleNow).data.y >= 6400 || (*particleNow).data.y <= -6400)
				{
					(*particlePrev).next = (*particleNow).next;
					delete(particleNow);
					particleNow = (*particlePrev).next;
					continue;
				}
				particleNow = (*particleNow).next;
			}
		}
	}

	saveData(grid, 1200);
	return;
};
