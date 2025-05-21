#ifndef PARSEVASP_H
#define PARSEVASP_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

void displaceAtomsVasp(string &inputFile, string &outputFilePath, double a, double b, double burgers, double loc1, double loc2, double nu, int N);

#endif