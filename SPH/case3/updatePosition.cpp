#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "iostream"
#include "fstream"
#include "particle.h"
#include "config.h"


void updatePosition(struct LNode** grid, struct LNode** gridCopy) {/*每次计算完后再更新位置*/
	int rowNum = (int)grid[0][0].data.row;
	int lineNum = (int)grid[0][0].data.line;
	struct LNode* particleNow;
	struct LNode* particlePrev;
	struct LNode* newPar;
	double newX;
	double newY;
	int newRow;
	int newLine;
	for (int j = 0; j < rowNum; j++) {
		for (int k = 0; k < lineNum; k++) {
			particleNow = grid[j][k].next;
			particlePrev = &(grid[j][k]);
			while (particleNow != NULL) {
				if ((*particleNow).data.type == itfluid) {
					newX = (*particleNow).data.x + (*particleNow).data.u * timeStep;
					newY = (*particleNow).data.y + (*particleNow).data.v * timeStep;
					(*particleNow).data.changeCoordinate(newX, newY);
					newRow = (int)((newY - calculate_bottom_boundary) / (2 * h));
					newLine = (int)((newX - calculate_left_boundary) / (2 * h));
					(*particleNow).data.changeGridNumber(newRow, newLine);
					double u = - PI / 6000 * ((*particleNow).data.y - 0);
					double v = PI / 6000 * ((*particleNow).data.x - 0);
					double r = sqrt(powf(((*particleNow).data.x - 0), 2) + powf(((*particleNow).data.y - 0), 2));

					(*particleNow).data.changeVelocity(u, v);

					if (r >= 6400 * sqrt(2))
					{
						(*particlePrev).next = (*particleNow).next;
						delete(particleNow);
						particleNow = (*particlePrev).next;
						continue;
					}
					else if (newX < coordination[0][0] || newX > coordination[2][0] || newY < coordination[0][0] || newY > coordination[2][1])
					{
						(*particleNow).data.changeType(itoutlet);
						(*particleNow).data.changeT(0.0f);
					}
					else if (((*particleNow).data.x < coordination[0][0] && (*particleNow).data.y <= 0.5 * coordination[2][1]) || ((*particleNow).data.x <= 0.5 * coordination[2][0] && (*particleNow).data.y > coordination[2][1]) || ((*particleNow).data.x > coordination[2][0] && (*particleNow).data.y >= 0.5 * coordination[2][1]) || ((*particleNow).data.x >= 0.5 * coordination[2][0] && (*particleNow).data.y < coordination[0][0]))
					{
						(*particleNow).data.changeType(itinlet);
						(*particleNow).data.changeT(0.0f);
					}
				}
				else if ((*particleNow).data.type == itoutlet) {
					newX = (*particleNow).data.x + (*particleNow).data.u * timeStep;
					newY = (*particleNow).data.y + (*particleNow).data.v * timeStep;
					double u = - PI / 6000 * ((*particleNow).data.y - 0);
					double v = PI / 6000 * ((*particleNow).data.x - 0);
					double r = sqrt(powf(((*particleNow).data.x - 0), 2) + powf(((*particleNow).data.y - 0), 2));
					(*particleNow).data.changeVelocity(u, v);
					//先判断！如果粒子流出流体域，则删除
					if (r >= 6400 * sqrt(2))
					{
						(*particlePrev).next = (*particleNow).next;
						delete(particleNow);
						particleNow = (*particlePrev).next;
						continue;
					}
					else if (((*particleNow).data.x < coordination[0][0] && (*particleNow).data.y <= 0.5 * coordination[2][1]) || ((*particleNow).data.x <= 0.5 * coordination[2][0] && (*particleNow).data.y > coordination[2][1]) || ((*particleNow).data.x > coordination[2][0] && (*particleNow).data.y >= 0.5 * coordination[2][1]) || ((*particleNow).data.x >= 0.5 * coordination[2][0] && (*particleNow).data.y < coordination[0][0]))
					{
						(*particleNow).data.changeType(itinlet);
						(*particleNow).data.changeCoordinate(newX, newY);

					}
					else if (newX >= coordination[0][0] && newX <= coordination[2][0] && newY >= coordination[0][0] && newY <= coordination[2][1])
					{
						(*particleNow).data.changeType(itfluid);
					}
					else
					{
						(*particleNow).data.changeCoordinate(newX, newY);
					}
					newRow = (int)((newY - calculate_bottom_boundary) / (2 * h));
					newLine = (int)((newX - calculate_left_boundary) / (2 * h));
					(*particleNow).data.changeGridNumber(newRow, newLine);

				}
				else if ((*particleNow).data.type == itinlet) {
					newX = (*particleNow).data.x + (*particleNow).data.u * timeStep;
					newY = (*particleNow).data.y + (*particleNow).data.v * timeStep;
					double u = - PI / 6000 * ((*particleNow).data.y - 0);
					double v = PI / 6000 * ((*particleNow).data.x - 0);
					double r = sqrt(powf(((*particleNow).data.x - 0), 2) + powf(((*particleNow).data.y - 0), 2));
					(*particleNow).data.changeVelocity(u, v);
					if (r >= 6400 * sqrt(2))
					{
						(*particlePrev).next = (*particleNow).next;
						delete(particleNow);
						particleNow = (*particlePrev).next;
						continue;
					}
					else if (newX >= coordination[0][0] && newX <= coordination[2][0] && newY >= coordination[0][0] && newY <= coordination[2][1])
					{
						(*particleNow).data.changeType(itfluid);
					}
					else if (((*particleNow).data.x < coordination[0][0] && (*particleNow).data.y <= 0.5 * coordination[2][1]) || ((*particleNow).data.x <= 0.5 * coordination[2][0] && (*particleNow).data.y > coordination[2][1]) || ((*particleNow).data.x > coordination[2][0] && (*particleNow).data.y >= 0.5 * coordination[2][1]) || ((*particleNow).data.x >= 0.5 * coordination[2][0] && (*particleNow).data.y < coordination[0][0]))
					{
						(*particleNow).data.changeType(itoutlet);
					}
					(*particleNow).data.changeCoordinate(newX, newY);
					newRow = (int)((newY - calculate_bottom_boundary) / (2 * h));
					newLine = (int)((newX - calculate_left_boundary) / (2 * h));
					(*particleNow).data.changeGridNumber(newRow, newLine);
				}
				/*else if((*particleNow).data.type == itboundary)
				{
					newX = (*particleNow).data.x + (*particleNow).data.u * timeStep;
					newY = (*particleNow).data.y + (*particleNow).data.v * timeStep;
					(*particleNow).data.changeCoordinate(newX, newY);
					newRow = (int)((newY - calculate_bottom_boundary) / (2 * h));
					newLine = (int)((newX - calculate_left_boundary) / (2 * h));
					(*particleNow).data.changeGridNumber(newRow, newLine);
				}*/
				else {

					newRow = (*particleNow).data.row;
					newLine = (*particleNow).data.line;
				}
				newPar = new LNode;
				(*newPar).data = (*particleNow).data;
				(*newPar).next = gridCopy[newRow][newLine].next;
				gridCopy[newRow][newLine].next = newPar;
				particleNow = (*particleNow).next;
				particlePrev = (*particlePrev).next;
			}
		}
	}
}
