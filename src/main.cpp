//
//  main.cpp
//  dislocationStructures
//
//  Created by Ian Anderson on 4/29/25.
//

#include <iostream>
#include <string>
#include "ParseLAMMPS.h"

using namespace std;

int main()
{
    double a;
    double b;
    double burgers;
    double x1 = 0, y1 = 0;
    double x2 = 0, y2 = 0;
    double nu = 0;
    double bwidth = 0;
    int N;
    int dislocationType;
    string inputFile;
    string outputFile;

    cout << "Enter type of dislocation:\n";
    cout << "  0: Edge dipole\n";
    cout << "  1: Screw dipole\n";
    cout << "  2: Single edge\n";
    cout << "  3: Single screw\n";
    cin >> dislocationType;

    cout << "Enter length of a vector: \n";
    cin >> a;

    cout << "Enter length of b vector: \n";
    cin >> b;

    cout << "Enter length of burgers vector: \n";
    cin >> burgers;

    if (dislocationType == 0 || dislocationType == 2)
    {
        cout << "Enter Poisson ratio: \n";
        cin >> nu;
    }

    if (dislocationType == 2)
    {
        cout << "Enter boundary width at free surfaces (same units as coordinates): \n";
        cin >> bwidth;
    }

    if (dislocationType == 0 || dislocationType == 1)
    {
        cout << "Enter x location of dislocation 1: \n";
        cin >> x1;

        cout << "Enter y location of dislocation 1: \n";
        cin >> y1;

        cout << "Enter x location of dislocation 2: \n";
        cin >> x2;

        cout << "Enter y location of dislocation 2: \n";
        cin >> y2;
    }
    else
    {
        cout << "Enter x location of dislocation: \n";
        cin >> x1;

        cout << "Enter y location of dislocation: \n";
        cin >> y1;
    }

    cout << "Enter number of images: \n";
    cin >> N;

    cout << "Enter path to input file: \n";
    cin >> inputFile;

    cout << "Enter path to output file: \n";
    cin >> outputFile;

    displaceAtoms(dislocationType, inputFile, outputFile, a, b, burgers, x1, y1, x2, y2, nu, N, bwidth);

    return 0;
}
