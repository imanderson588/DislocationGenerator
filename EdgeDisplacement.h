#ifndef EDGEDISPLACEMENT_H
#define EDGEDISPLACEMENT_H

#include <string>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

const double pi = M_PI;

double edgeDisplacement_x(double x, double y, double loc1, double loc2, double burgers, double nu);

double edgeImage_x(double x, double y, double a, double b, double loc1, double loc2, double burgers, double nu, int N);

double edgeCorrection_x(double x, double y, double a, double b, double loc1, double loc2, double burgers, double nu, int N);

double edgeDisplacement_y(double x, double y, double loc1, double loc2, double burgers, double nu);


#endif

