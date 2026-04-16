#include "ParseLAMMPS.h"
#include "EdgeDisplacement.h"
#include "ScrewDisplacement.h"
#include "SingleDislocations.h"
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

vector<string> splitBySpaces(const string &line)
{
    istringstream iss(line);
    vector<string> tokens;
    string word;

    while (iss >> word)
    {
        tokens.push_back(word);
    }

    return tokens;
}

string recombine(const vector<string> &words)
{
    ostringstream oss;
    for (size_t i = 0; i < words.size(); ++i)
    {
        oss << words[i];
        if (i != words.size() - 1)
        {
            oss << " ";
        }
    }
    return oss.str();
}

static bool atomOverlapsDislocation(double ax, double ay, double x1, double y1, double x2, double y2, int dislocationType, double tol = 1e-6)
{
    double dx1 = ax - x1, dy1 = ay - y1;
    if (dx1 * dx1 + dy1 * dy1 < tol * tol)
        return true;

    if (dislocationType == 0 || dislocationType == 1)
    {
        double dx2 = ax - x2, dy2 = ay - y2;
        if (dx2 * dx2 + dy2 * dy2 < tol * tol)
            return true;
    }

    return false;
}

void displaceAtoms(int dislocationType, const string &inputFile, const string &outputFilePath, double a, double b, double burgers, double x1, double y1, double x2, double y2, double nu, int N, double bwidth)
{
    string lineContents;
    ifstream inputData(inputFile);

    // First pass: check that no atomic position coincides with a dislocation core
    {
        int checkFlag = 0;
        while (getline(inputData, lineContents))
        {
            if (lineContents == "Atoms # atomic")
            {
                checkFlag = 1;
                continue;
            }
            if (lineContents == "Velocities")
                break;
            if (lineContents.empty() || !checkFlag)
                continue;

            vector<string> words = splitBySpaces(lineContents);
            if (words.size() < 4)
                continue;
            double ax = stod(words[2]);
            double ay = stod(words[3]);

            if (atomOverlapsDislocation(ax, ay, x1, y1, x2, y2, dislocationType))
            {
                cerr << "Error: dislocation core overlaps with an atomic position (atom at x=" << ax << ", y=" << ay << "). Aborting.\n";
                inputData.close();
                return;
            }
        }
        // Rewind for second pass
        inputData.clear();
        inputData.seekg(0);
    }

    ofstream outputFile(outputFilePath);
    int modifyFlag = 0;

    while (getline(inputData, lineContents))
    {
        if (lineContents == "Atoms # atomic")
        {
            outputFile << lineContents << "\n";
            modifyFlag = 1;
            continue;
        }

        if (lineContents == "Velocities")
        {
            outputFile << lineContents << "\n";
            modifyFlag = 0;
            continue;
        }

        if (lineContents == "")
        {
            outputFile << lineContents << "\n";
            continue;
        }

        if (modifyFlag)
        {
            vector<string> words = splitBySpaces(lineContents);

            double x_value = stod(words[2]);
            double y_value = stod(words[3]);
            double z_value = stod(words[4]);

            if (dislocationType == 0)
            {
                // Edge dipole displacements
                double u_x = totEdge_x(x_value, y_value, a, b, x1, y1, x2, y2, burgers, nu, N);
                double u_y = totEdge_y(x_value, y_value, a, b, x1, y1, x2, y2, burgers, nu, N);
                words[2] = to_string(x_value + u_x);
                words[3] = to_string(y_value + u_y);
            }
            else if (dislocationType == 1)
            {
                // Screw dipole displacements
                double u_z = screwDipole(x_value, y_value, a, b, burgers, x1, x2, N);
                words[4] = to_string(z_value + u_z);
            }
            else if (dislocationType == 2)
            {
                // Single edge dislocation at (x1, y1)
                // x1 is position along Burgers vector (file y-col), y1 along glide plane normal (file z-col)
                double dx = x_value - x1;
                double dy = y_value - y1;
                double u_x = singleEdgeDisplacement_x(dx, dy, burgers, nu);
                double u_y = singleEdgeDisplacement_y(dx, dy, burgers, nu);
                words[2] = to_string(x_value + u_x);
                words[3] = to_string(y_value + u_y);
            }
            else if (dislocationType == 3)
            {
                // Single screw dislocation at (x1, y1)
                double u_z = singleScrew(x_value - x1, y_value - y1, burgers);
                words[4] = to_string(z_value + u_z);
            }

            string newLine = recombine(words);
            outputFile << newLine << "\n";
        }
        else
        {
            outputFile << lineContents << "\n";
        }
    }

    inputData.close();
    outputFile.close();
}
