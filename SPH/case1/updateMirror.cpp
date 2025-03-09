#include "stdio.h"
#include "stdlib.h"
#include "particle.h"
#include "config.h"

void updateMirror(struct LNode** grid, struct LNode** gridCopy) {
	double x1 = coordination[0][0];
	double y1 = coordination[0][1];
	double x2 = coordination[1][0];
	double y2 = coordination[1][1];
	double x3 = coordination[2][0];
	double y3 = coordination[2][1];
	double x4 = coordination[3][0];
	double y4 = coordination[3][1];
	double dis = initialDistance;   /*---粒子间距---*/
	double len = 2 * h;   /*---流体域边界外设两层粒子---*/
	double gNode_x1 = x1 - len;   /*---边界四个顶点的坐标---*/
	double gNode_y1 = y1 - len;
	double gNode_x2 = x2 + len;
	double gNode_y2 = gNode_y1;
	double gNode_x3 = gNode_x2;
	double gNode_y3 = y3 + len;
	double gNode_x4 = gNode_x1;
	double gNode_y4 = gNode_y3;
	int rowNum = (int)grid[0][0].data.row;
	int lineNum = (int)grid[0][0].data.line;
	int k, j;
	struct LNode* particleNow;
	struct LNode* mirrorPar;
	double x, y;
	for (k = 0; k < rowNum; k++) {
		for (j = 0; j < lineNum; j++)
		{
			particleNow = grid[k][j].next;
			while (particleNow != NULL)
			{
				if ((*particleNow).data.type == itfluid || (*particleNow).data.type == itinlet || (*particleNow).data.type == itoutlet) {
					mirrorPar = new LNode;
					

					
				}
			}
		}
	}
}