#include "ParseLAMMPS.h"
#include "EdgeDisplacement.h"
#include "ScrewDisplacement.h"
#include "SingleDislocations.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

static constexpr double PI = M_PI;

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

static bool positionOverlapsAnyAtom(double px, double py, const vector<double> &atomX, const vector<double> &atomY, double tol)
{
    for (size_t i = 0; i < atomX.size(); ++i)
    {
        double dx = atomX[i] - px, dy = atomY[i] - py;
        if (dx * dx + dy * dy < tol * tol)
            return true;
    }
    return false;
}

static void adjustPosition(double &x, double &y, const vector<double> &atomX, const vector<double> &atomY, double tol, double step)
{
    for (int ring = 1; ring <= 10000; ++ring)
    {
        double r = ring * step;
        int nAngles = max(8, ring * 4);
        for (int k = 0; k < nAngles; ++k)
        {
            double angle = 2.0 * PI * k / nAngles;
            double tx = x + r * cos(angle);
            double ty = y + r * sin(angle);
            if (!positionOverlapsAnyAtom(tx, ty, atomX, atomY, tol))
            {
                x = tx;
                y = ty;
                return;
            }
        }
    }
    cerr << "Warning: could not find a non-overlapping position within search range. Proceeding with original coordinates.\n";
}

void displaceAtoms(int dislocationType, const string &inputFile, const string &outputFilePath, double a, double b, double burgers, double x1, double y1, double x2, double y2, double nu, int N, double bwidth)
{
    string lineContents;
    ifstream inputData(inputFile);

    // First pass: collect all atom positions and adjust dislocation cores if needed
    {
        vector<double> atomX, atomY;
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
            atomX.push_back(stod(words[2]));
            atomY.push_back(stod(words[3]));
        }

        // Tolerance for overlap and step size for adjustment
        const double tol = 1e-6;
        const double step = 1e-4;

        if (positionOverlapsAnyAtom(x1, y1, atomX, atomY, tol))
        {
            cerr << "Warning: dislocation 1 position (" << x1 << ", " << y1
                 << ") overlaps with an atomic position.\n";
            adjustPosition(x1, y1, atomX, atomY, tol, step);
            cerr << "         Adjusted dislocation 1 position to (" << x1 << ", " << y1 << ").\n";
        }

        if (dislocationType == 0 || dislocationType == 1)
        {
            if (positionOverlapsAnyAtom(x2, y2, atomX, atomY, tol))
            {
                cerr << "Warning: dislocation 2 position (" << x2 << ", " << y2
                     << ") overlaps with an atomic position.\n";
                adjustPosition(x2, y2, atomX, atomY, tol, step);
                cerr << "         Adjusted dislocation 2 position to (" << x2 << ", " << y2 << ").\n";
            }
        }

        // Rewind for displacement pass
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
