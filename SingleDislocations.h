#ifndef SINGLEDISLOCATIONS_H
#define SINGLEDISLOCATIONS_H

#include <string>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

double singleEdgeDisplacement_x(double x, double y, double burgers, double nu);
double singleEdgeDisplacement_y(double x, double y, double burgers, double nu);

#endif