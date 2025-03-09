#pragma once
using namespace std;

static double initialDistance = 128;   /*---初始化粒子间距---*/
static double initialDistancey = 128;
static double coordination[4][2] = { {-9088, -9088},{9088, -9088},{9088, 9088},{-9088, 9088}};   /*---流体域坐标区域---*/
static double smooth_r = 1.33;
static double h = smooth_r * initialDistance;
static double hy = smooth_r * initialDistancey;
static double heatConductivity = 1.0;
static double Pe = 10;
static double Pr = 0.5;
static double Re = 10;
static double Ec = 0.1;
static int problemType = 1;
static double calculate_left_boundary = coordination[0][0] - 2 * h;
static double calculate_bottom_boundary = coordination[0][1] - 2 * hy;
static double timeStep = 1.0f;
static double timeEnd = 12000.0f;
static int model = 3;
static int timeSum = (int)(timeEnd / timeStep);
static double PI = 3.1415926535898;