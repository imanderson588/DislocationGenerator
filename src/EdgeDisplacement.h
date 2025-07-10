#ifndef EDGEDISPLACEMENT_H
#define EDGEDISPLACEMENT_H

#include <string>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

double edgeDisplacement_x(double x, double y, double x1, double y1, double x2, double y2, double burgers, double nu);

double edgeImage_x(double x, double y, double a, double b, double x1, double y1, double x2, double y2, double burgers, double nu, int N);

double edgeCorrection_x(double x, double y, double a, double b, double x1, double y1, double x2, double y2, double burgers, double nu, int N);

double totEdge_x(double x, double y, double a, double b, double x1, double y1, double x2, double y2, double burgers, double nu, int N);

double edgeDisplacement_y(double x, double y, double x1, double y1, double x2, double y2, double burgers, double nu);

double edgeImage_y(double x, double y, double a, double b, double x1, double y1, double x2, double y2, double burgers, double nu, int N);

double edgeDipoleTilt(double x, double y, double a, double b, double burgers, double loc1, double loc2, int N);

double edgeCorrection_y(double x, double y, double a, double b, double x1, double y1, double x2, double y2, double burgers, double nu, int N);

double totEdge_y(double x, double y, double a, double b, double x1, double y1, double x2, double y2, double burgers, double nu, int N);

#endif