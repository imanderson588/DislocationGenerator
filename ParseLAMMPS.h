#ifndef PARSELAMMPS_H
#define PARSELAMMPS_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

void displaceAtoms(string &inputFile, string &outputFilePath, double a, double b, double burgers, double x1, double y1, double x2, double y2, double nu, int N);

vector<string> splitBySpaces(const string &line);

string recombine(const vector<string> &words);

#endif