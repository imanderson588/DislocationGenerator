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

static void adjustPosition(double &x, double &y, const vector<double> &atomX, const vector<double> &atomY, double tol, double step,
                           double xlo, double xhi, double ylo, double yhi)
{
    // Search only within a local radius (20 steps = 2*a) to stay near the given position.
    // The interstitial region between neighboring atoms is always within a few fractions of
    // the lattice parameter, so there is no need to search further than this.
    const int maxRings = 20;
    for (int ring = 1; ring <= maxRings; ++ring)
    {
        double r = ring * step;
        int nAngles = max(8, ring * 4);
        for (int k = 0; k < nAngles; ++k)
        {
            double angle = 2.0 * PI * k / nAngles;
            double tx = x + r * cos(angle);
            double ty = y + r * sin(angle);
            if (tx < xlo || tx > xhi || ty < ylo || ty > yhi)
                continue;
            if (!positionOverlapsAnyAtom(tx, ty, atomX, atomY, tol))
            {
                x = tx;
                y = ty;
                return;
            }
        }
    }
    cerr << "Warning: could not find a non-overlapping interstitial position within " << maxRings * step
         << " of the given coordinates. Proceeding with original coordinates.\n";
}

void displaceAtoms(int dislocationType, const string &inputFile, const string &outputFilePath, double a, double b, double burgers, double x1, double y1, double x2, double y2, double nu, int N, double bwidth)
{
    string lineContents;
    ifstream inputData(inputFile);

    // First pass: collect cell bounds and atom positions, then adjust dislocation cores if needed
    {
        vector<double> atomX, atomY;
        double xlo = -1e300, xhi = 1e300, ylo = -1e300, yhi = 1e300;
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

            if (!checkFlag)
            {
                // Parse box bounds from header lines like "lo hi xlo xhi"
                vector<string> words = splitBySpaces(lineContents);
                if (words.size() == 4 && words[2] == "xlo" && words[3] == "xhi")
                {
                    xlo = stod(words[0]);
                    xhi = stod(words[1]);
                }
                else if (words.size() == 4 && words[2] == "ylo" && words[3] == "yhi")
                {
                    ylo = stod(words[0]);
                    yhi = stod(words[1]);
                }
            }

            if (lineContents.empty() || !checkFlag)
                continue;

            vector<string> words = splitBySpaces(lineContents);
            if (words.size() < 4)
                continue;
            atomX.push_back(stod(words[2]));
            atomY.push_back(stod(words[3]));
        }

        // Compute the minimum 2-D nearest-neighbour distance between atoms so that
        // the overlap tolerance and search step are always correctly scaled to the
        // actual inter-atomic spacing rather than to the (much larger) supercell
        // dimension `a`.  Using `a` directly caused the tolerance to be orders of
        // magnitude too large for dense structures (e.g. rocksalt), making it
        // impossible to find any non-overlapping interstitial placement.
        constexpr double LARGE_DISTANCE = 1e300;
        double minDist2 = LARGE_DISTANCE;
        for (size_t i = 0; i < atomX.size(); ++i)
        {
            for (size_t j = i + 1; j < atomX.size(); ++j)
            {
                double dx = atomX[i] - atomX[j];
                double dy = atomY[i] - atomY[j];
                double d2 = dx * dx + dy * dy;
                if (d2 < minDist2)
                    minDist2 = d2;
            }
        }
        double minDist = (minDist2 < LARGE_DISTANCE && atomX.size() > 1) ? sqrt(minDist2) : a;

        // Tolerance for overlap and step size for adjustment (relative to nearest-neighbour distance)
        const double tol = 0.3 * minDist;
        const double step = 0.1 * minDist;

        if (positionOverlapsAnyAtom(x1, y1, atomX, atomY, tol))
        {
            cerr << "Warning: dislocation 1 position (" << x1 << ", " << y1
                 << ") overlaps with an atomic position.\n";
            adjustPosition(x1, y1, atomX, atomY, tol, step, xlo, xhi, ylo, yhi);
            cerr << "         Adjusted dislocation 1 position to (" << x1 << ", " << y1 << ").\n";
        }

        if (dislocationType == 0 || dislocationType == 1)
        {
            if (positionOverlapsAnyAtom(x2, y2, atomX, atomY, tol))
            {
                cerr << "Warning: dislocation 2 position (" << x2 << ", " << y2
                     << ") overlaps with an atomic position.\n";
                adjustPosition(x2, y2, atomX, atomY, tol, step, xlo, xhi, ylo, yhi);
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
