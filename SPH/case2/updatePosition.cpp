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
					if (newX > coordination[2][0] || newY > coordination[2][1]) {
						(*particleNow).data.changeType(itoutlet);
					}
				}
				else if ((*particleNow).data.type == itoutlet) {
					newX = (*particleNow).data.x + (*particleNow).data.u * timeStep;
					newY = (*particleNow).data.y + (*particleNow).data.v * timeStep;
					if (newX > coordination[2][0] + 2 * h || newY > coordination[2][1] + 2 * h) {
						(*particlePrev).next = (*particleNow).next;
						delete(particleNow);
						particleNow = (*particlePrev).next;
						continue;
					}
					/*(*particleNow).data.changeCoordinate(newX, newY);
					(*particleNow).data.changeVelocity(1.0f, 1.0f);
					newRow = (int)((newY - calculate_bottom_boundary) / (2 * h));
					newLine = (int)((newX - calculate_left_boundary) / (2 * h));
					(*particleNow).data.changeGridNumber(newRow, newLine);*/
					
				}
				else if ((*particleNow).data.type == itinlet) {
					newX = (*particleNow).data.x + (*particleNow).data.u * timeStep;
					newY = (*particleNow).data.y + (*particleNow).data.v * timeStep;
					(*particleNow).data.changeCoordinate(newX, newY);
					newRow = (int)((newY - calculate_bottom_boundary) / (2 * h));
					newLine = (int)((newX - calculate_left_boundary) / (2 * h));
					(*particleNow).data.changeGridNumber(newRow, newLine);
					if (newX >= coordination[0][0] && newY >= coordination[0][1]) {
						if (newX <= coordination[2][0] && newY <= coordination[2][1]) {
							(*particleNow).data.changeType(itfluid);
						}
						else {
							(*particleNow).data.changeType(itoutlet);
						}
						newPar = new LNode;
						(*newPar).data.changeType(itinlet);
						(*newPar).data.changeCoordinate(newX - 2 * initialDistance, newY - 2 * initialDistance);
						(*newPar).data.changeVelocity(0.2f, 0.0f);
						newRow = (int)(((*newPar).data.y - calculate_bottom_boundary) / (2 * h));
						newLine = (int)(((*newPar).data.x - calculate_left_boundary) / (2 * h));
						(*newPar).data.changeGridNumber(newRow, newLine);
						(*newPar).next = gridCopy[newRow][newLine].next;
						gridCopy[newRow][newLine].next = newPar;
					}
				}
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
