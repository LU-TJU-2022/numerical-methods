#pragma once
using namespace std;

static double initialDistance = 5;   /*---��ʼ�����Ӽ��---*/
static double coordination[4][2] = { {0, 0},{500, 0},{500, 200},{0, 200} };   /*---��������������---*/
static double smooth_r = 1.33;
static double h = smooth_r * initialDistance;
static double heatConductivity = 0.1;
static int problemType = 1;
static double calculate_left_boundary = coordination[0][0] - 2 * h;
static double calculate_bottom_boundary = coordination[0][1] - 2 * h;
static double timeStep = 1.0f;
static double timeEnd = 3500.0f;
static int model = 4;
static int timeSum = (int)(timeEnd / timeStep);
static double PI = 3.1415926535898;