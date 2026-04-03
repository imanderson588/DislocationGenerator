#ifndef PARSEVASP_H
#define PARSEVASP_H

#include <string>

void displaceAtomsVasp(const std::string& inputFile, const std::string& outputFilePath, double a, double b, double burgers, double x1, double y1, double x2, double y2, double nu, int N);

#endif
