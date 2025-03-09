#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "iostream"
#include "fstream"
#include "particle.h"
#include "config.h"

void gridReset(struct LNode** grid, struct LNode** gridCopy) {
    int rowNum = (int)grid[0][0].data.row;
    int lineNum = (int)grid[0][0].data.line;
    int k, j;
    struct LNode* particleNow;
    struct LNode* particlePrev;
    struct LNode* newPar;
    for (k = 0; k < rowNum; k++) {
        for (j = 0; j < lineNum; j++)
        {
            particlePrev = &(grid[k][j]);
            particleNow = grid[k][j].next;
            while (particleNow != NULL)
            {
                (*particlePrev).next = (*particleNow).next;
                delete(particleNow);
                particleNow = (*particlePrev).next;
            }
        }
    }
    for (k = 0; k < rowNum; k++) {
        for (j = 0; j < lineNum; j++)
        {
            particlePrev = &(gridCopy[k][j]);
            particleNow = gridCopy[k][j].next;
            while (particleNow != NULL)
            {
                newPar = new LNode;  
                (*newPar).data = (*particleNow).data;
                (*newPar).next = grid[(*particleNow).data.row][(*particleNow).data.line].next;
                grid[(*particleNow).data.row][(*particleNow).data.line].next = newPar;

                (*particlePrev).next = (*particleNow).next;
                delete(particleNow);
                particleNow = (*particlePrev).next;
            }
        }
    }
}