#pragma once
#include "particle.h"
#include <vector>
using namespace std;

vector<double> calculate_heatFlux_by_Model2(struct LNode* nearPar, particle centerPar);
double calculate_T_by_Model2(struct LNode* nearPar, particle centerPar);

vector<double> calculate_heatFlux_by_Model4(struct LNode* nearPar, particle centerPar);
double calculate_T_by_Model4(struct LNode* nearPar, particle centerPar);

vector<double> calculate_heatFlux_by_Model4_1(struct LNode* nearPar, particle centerPar);
double calculate_T_by_Model4_1(struct LNode* nearPar, particle centerPar);