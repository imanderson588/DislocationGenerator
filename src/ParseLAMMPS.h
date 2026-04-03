#ifndef PARSELAMMPS_H
#define PARSELAMMPS_H

#include <string>
#include <vector>

void displaceAtoms(int dislocationType, const std::string& inputFile, const std::string& outputFilePath, double a, double b, double burgers, double x1, double y1, double x2, double y2, double nu, int N, double bwidth = 0.0);

std::vector<std::string> splitBySpaces(const std::string& line);

std::string recombine(const std::vector<std::string>& words);

#endif
