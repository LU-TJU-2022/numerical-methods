#pragma once
using namespace std;

static double initialDistance = 0.09375;   /*---初始化粒子间距---*/
static double coordination[4][2] = { {-3, -3},{3, -3},{3, 3},{-3, 3} };   /*---流体域坐标区域---*/
static double smooth_r = 1.33;
static double h = smooth_r * initialDistance;
static double heatConductivity = 0.0000;
static int problemType = 1;
static double calculate_left_boundary = coordination[0][0] - 2 * h;
static double calculate_bottom_boundary = coordination[0][1] - 2 * h;
static double timeStep = 0.005f;
static double timeEnd = 2.8f;
static int model = 3;
static int timeSum = (int)(timeEnd / timeStep);
static double PI = 3.1415926535898;