#ifndef SINGLEDISLOCATIONS_H
#define SINGLEDISLOCATIONS_H

#include <string>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

double singleEdgeDisplacement_x(double x, double y, double burgers, double nu);
double singleEdgeImage_x(double x, double y, double a, double b, double burgers, double nu, int N);
double singleEdgeCorrection_x(double x, double y, double a, double b, double burgers, double nu, int N);
double totSingleEdge_x(double x, double y, double a, double b, double burgers, double nu, int N);

double singleEdgeDisplacement_y(double x, double y, double burgers, double nu);
double singleEdgeImage_y(double x, double y, double a, double b, double burgers, double nu, int N);
double singleEdgeCorrection_y(double x, double y, double a, double b, double burgers, double nu, int N);
double totSingleEdge_y(double x, double y, double a, double b, double burgers, double nu, int N);

#endif