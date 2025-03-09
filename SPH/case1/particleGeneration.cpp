#include "particle.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include <iostream>

using namespace std;
struct LNode** Generation(double wall[4][2], double initialDistance, double h) {
	double x1 = wall[0][0];
	double y1 = wall[0][1];
	double x2 = wall[1][0];
	double y2 = wall[1][1];
	double x3 = wall[2][0];
	double y3 = wall[2][1];
	double x4 = wall[3][0];
	double y4 = wall[3][1];
	double dis = initialDistance;   /*---粒子间距---*/
	double len = 2 * h;   /*---流体域边界外设两层粒子---*/
	double gNode_x1 = x1 - len;   /*---边界四个顶点的坐标---*/
	double gNode_y1 = y1 - len;
	double gNode_x2 = x2 + len;
	double gNode_y2 = gNode_y1;
	double gNode_x3 = gNode_x2;
	double gNode_y3 = y3 +len;
	double gNode_x4 = gNode_x1;
	double gNode_y4 = gNode_y3;

	int rowNum = (int)(fabs(gNode_y4 - gNode_y1) / len) + 1;   
	int lineNum = (int)(fabs(gNode_x2 - gNode_x1) / len) + 1;

	struct LNode** grid;  
	grid = new LNode * [rowNum];
	for (int n = 0; n < rowNum; n++) {
		grid[n] = new LNode[lineNum];
		for (int m = 0; m < lineNum; m++) {
			grid[n][m].next = NULL;
		}
	}
	int row;
	int line;
	//int count=0;
	double x, y;
	struct LNode* p;
	for (x = x1 - 2 * dis; x <= gNode_x2; x += dis) {
		for (y = y1 - 2 * dis; y <= gNode_y4; y += dis){
			row = (int)((y - gNode_y1) / len);
			line = (int)((x - gNode_x1) / len);
			p = new LNode;
			(*p).data.changeCoordinate(x, y); 
			(*p).data.changeGridNumber(row, line);
			if (x >= x1 && x <= x2 && y >= y1 && y <= y4) {

				//(*p).data.T = exp(-powf((x - 2), 2) / 2 / 0.1 / 0.1 - powf((y - 2), 2) / 2 / 0.4 / 0.4);
				if (x >= 3 * x1 / 4 && x <= x1 / 4 && y >= 3 * y1 / 4 && y <= y1 / 4)
				{
					(*p).data.T = 10.0f;
				}

				(*p).data.changeType(itfluid);
				(*p).data.changeVelocity(1.0f, 1.0f);
			}
			else if(!(y >= y1 && x >= x1)){
				if (x >= gNode_x1 && y >= gNode_y1 &&  y - y4 <= x - x4 && y - y2 >= x - x2) {
					(*p).data.T = 0.0f;
					(*p).data.changeType(itinlet);
					(*p).data.changeVelocity(1.0f, 1.0f);
				}
				else {
					(*p).data.T = 0.0f;
					(*p).data.changeType(itinlet);
					(*p).data.changeVelocity(0.0f, 0.0f);
				}
			}
			else if(!(y <= y3 && x <= x3)) {
				if (x <= gNode_x3 && y <= gNode_y3 &&  y-y4<=x - x4  && y - y2>= x-x2) {
					(*p).data.T = 0.0f;
					/*(*p).data.changeType(itinlet);*/
					(*p).data.changeType(itoutlet);
					(*p).data.changeVelocity(1.0f, 1.0f);
				}
				else {
					(*p).data.T = 0.0f;
					/*(*p).data.changeType(itinlet);*/
					(*p).data.changeType(itoutlet);
					(*p).data.changeVelocity(0.0f, 0.0f);
				}
			}
			else {
				(*p).data.T = 0.0f;
				(*p).data.changeType(itwall);
				(*p).data.changeVelocity(0.0f, 0.0f);
			}
			(*p).next = grid[row][line].next;
			grid[row][line].next = p;
		}
	}
	grid[0][0].data.row = rowNum;
	grid[0][0].data.line = lineNum;
	return grid;
};