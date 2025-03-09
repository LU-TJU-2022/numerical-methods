#pragma once
#include "particle.h"

vector<double> calculate_heatFlux_by_Model4_1(struct LNode* nearPar, particle centerPar);
double calculate_T_by_Model4_1(struct LNode* nearPar, particle centerPar);