#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "iostream"
#include "fstream"
#include "particle.h"

using namespace std;

struct LNode* particleSearch(struct LNode** grid, particle center, double distance) {
	struct LNode* neighborParticle = NULL;
	struct LNode particleCur;
	int rowNum = (int)grid[0][0].data.row;
	int lineNum = (int)grid[0][0].data.line;
	LNode* p;
	double distanceTrue;
	for (int j = -1; j <= 1; j++) {
		if (center.row + j < 0 || center.row + j >= rowNum) { 
			continue;
		}
		for (int k = -1; k <= 1; k++) {
			if (center.line + k < 0 || center.line + k >= lineNum) {
				continue;
			}
			particleCur = grid[center.row + j][center.line + k]; 
			while (particleCur.next != NULL) {
				particleCur = *(particleCur.next);
				distanceTrue = sqrt(pow(center.x - particleCur.data.x, 2) + pow(center.y - particleCur.data.y, 2));
				if (distanceTrue <= distance) {
					p = new LNode;
					(*p).data = particleCur.data;
					(*p).next = neighborParticle;
					neighborParticle = p;
				}
			}
		}
	}
	return neighborParticle;
}