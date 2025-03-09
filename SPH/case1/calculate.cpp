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

using namespace std;

void saveData(struct LNode** grid, int timeStep) {/*--保存计算结果--*/
	cout.setf(ios::fixed);
	cout.precision(10);
	char szName[100] = {'\0'};
	sprintf_s(szName, "SPH时间步%d结果case1.txt", timeStep);
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
				dataFile << (*particleNow).data.x << '\t' << (*particleNow).data.y << '\t' << (*particleNow).data.T << '\t' << (*particleNow).data.row << '\t' << (*particleNow).data.line << endl;
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
	grid = Generation(coordination, initialDistance, h);
	vector<double> newHeatFlux;
	double newT;
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
	for (int i = 1; i <= 560; i++) {
		cout <<"第"<< i<<"时间步" << endl;
		for (int j = 0; j < rowNum; j++) {
			for (int k = 0; k < lineNum; k++) {
				particleNow = grid[j][k].next;
				while (particleNow != NULL) {
					if ((*particleNow).data.type == itfluid) {
						particleNear = particleSearch(grid, (*particleNow).data, 2 * h);
						
						/*if (model == 2) {
							newHeatFlux = calculate_heatFlux_by_Model2(particleNear, (*particleNow).data);
							(*particleNow).data.heatFlux_x = newHeatFlux[0];
							(*particleNow).data.heatFlux_y = newHeatFlux[1];
						}*/
						if (model == 3) {
							newT = calculate_T_by_Model3(particleNear, (*particleNow).data);
							(*particleNow).data.heatFlux_x = newT;
							(*particleNow).data.changeT((*particleNow).data.heatFlux_x);
						}
						/*else {
							newHeatFlux = calculate_heatFlux_by_Model4_1(particleNear, (*particleNow).data);
							(*particleNow).data.heatFlux_x_left = newHeatFlux[0];
							(*particleNow).data.heatFlux_x_right = newHeatFlux[1];
							(*particleNow).data.heatFlux_y_left = newHeatFlux[2];
							(*particleNow).data.heatFlux_y_right = newHeatFlux[3];
						}*/
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
		/*for (int j = 0; j < rowNum; j++) {
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
							(*particleNow).data.changeT((*particleNow).data.heatFlux_x);
						}
						else {
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
		}*/
		updatePosition(grid, gridCopy);
		gridReset(grid, gridCopy);
		
			
		
	}
	saveData(grid, 560);
	return;
};
