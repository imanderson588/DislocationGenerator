#ifndef EDGEDISPLACEMENT_H
#define EDGEDISPLACEMENT_H

#include <string>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

double edgeDisplacement_x(double x, double y, double loc1, double loc2, double burgers, double nu);

double edgeDisplacement_y(double x, double y, double loc1, double loc2, double burgers, double nu);

#endif