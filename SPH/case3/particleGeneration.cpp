#include "particle.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "config.h"
#include <iostream>

using namespace std;
struct LNode** Generation(double wall[4][2], double initialDistance, double initialDistancey, double h, double hy) {
	double x1 = wall[0][0];
	double y1 = wall[0][1];
	double x2 = wall[1][0];
	double y2 = wall[1][1];
	double x3 = wall[2][0];
	double y3 = wall[2][1];
	double x4 = wall[3][0];
	double y4 = wall[3][1];
	double disx = initialDistance;   /*---粒子间距---*/
	double disy = initialDistancey;
	double r = 0;
	double lenx = 2 * h;   /*---流体域边界外设两层粒子---*/
	double leny = 2 * hy;
	double gNode_x1 = x1 - lenx;   /*---边界四个顶点的坐标---*/
	double gNode_y1 = y1 - leny;
	double gNode_x2 = x2 + lenx;
	double gNode_y2 = gNode_y1;
	double gNode_x3 = gNode_x2;
	double gNode_y3 = y3 + leny;
	double gNode_x4 = gNode_x1;
	double gNode_y4 = gNode_y3;

	int rowNum = (int)(fabs(gNode_y4 - gNode_y1) / lenx) + 1;   
	int lineNum = (int)(fabs(gNode_x2 - gNode_x1) / lenx) + 1;

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
	for (x = x1; x <= x2; x += disx) {
		for (y = y1; y <= y4; y += disy){
			row = (int)((y - gNode_y1) / lenx);
			line = (int)((x - gNode_x1) / lenx);
			p = new LNode;
			(*p).data.changeCoordinate(x, y); 
			(*p).data.changeGridNumber(row, line);
			double u = - PI / 6000 * ((*p).data.y - 0);
			double v = PI / 6000 * ((*p).data.x - 0);
			
			if (x > x1 - 1e-8 && x < x2 + 1e-8 && y > y1 - 1e-8 && y < y3 + 1e-8) {

				(*p).data.T = exp(- powf((x - 0), 2) / 2 / 707.15 / 707.15 - powf((y - 3200), 2) / 2 / 500 / 500);
				(*p).data.changeType(itfluid);
				r = sqrt(powf((x - 0), 2) + powf((y - 0), 2));
				(*p).data.changeVelocity(u, v);
			}
			else if (((x < coordination[0][0] && y > 0.5 * coordination[2][1]) || (x > 0.5 * coordination[2][0] && y > coordination[2][1]) || (x > coordination[2][0] && y < 0.5 * coordination[2][1]) || (x < 0.5 * coordination[2][0] && y < coordination[0][0])))
			{
				(*p).data.T = 0.0f;
				(*p).data.changeType(itoutlet);
				r = sqrt(powf((x - 0), 2) + powf((y - 0), 2));

			}
			else if (((x < coordination[0][0] && y <= 0.5 * coordination[2][1]) || (x <= 0.5 * coordination[2][0] && y > coordination[2][1]) || (x > coordination[2][0] && y >= 0.5 * coordination[2][1]) || (x >= 0.5 * coordination[2][0] && y < coordination[0][0])) && abs(x - x1) > 1e-8 && abs(x - x2) > 1e-8)
			{
				(*p).data.T = 0.0f;
				(*p).data.changeType(itinlet);
				r = sqrt(powf((x - 0), 2) + powf((y - 0), 2));

			}
			/*else {
				if (abs(y - y1) < 1e-8)
				{
					(*p).data.T = 0.0f;
					(*p).data.changeVelocity(0.0f, 0.0f);

				}
				else if (abs(y - y3) < 1e-8)
				{
					(*p).data.T = 0.0f;
					(*p).data.changeVelocity(0.0f, 0.0f);
				}
				else if (abs(x - x1) < 1e-8)
				{
					(*p).data.T = 0.0f;
					(*p).data.changeVelocity(0.0f, 0.0f);
				}
				else if (abs(x - x2) < 1e-8)
				{
					(*p).data.T = 0.0f;
					(*p).data.changeVelocity(0.0f, 0.0f);
				}
				/*else
				{
					(*p).data.T = 0.0f;
					(*p).data.changeVelocity(0.0f, 0.0f);
				}
				(*p).data.changeType(itboundary);
			}*/
			(*p).next = grid[row][line].next;
			grid[row][line].next = p;
		}
	}
	grid[0][0].data.row = rowNum;
	grid[0][0].data.line = lineNum;
	return grid;
};