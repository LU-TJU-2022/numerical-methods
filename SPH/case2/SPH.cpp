// SPH.cpp 
//

#include <iostream>
#include "calculate.h"
#include<ctime>
using namespace std;
int main()
{
    clock_t startTime, endTime;
    startTime = clock();
    cout << "Hello World!\n";
    calculate();
    endTime = clock();
    cout << "The run time is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
}

