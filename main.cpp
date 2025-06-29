//
//  main.cpp
//  dislocationStructures
//
//  Created by Ian Anderson on 4/29/25.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include "EdgeDisplacement.h"
#include "ScrewDisplacement.h"
#include "ParseLAMMPS.h"
#include "ParseVASP.h"
#include "SingleDislocations.h"

using namespace std;

int main()
{
    double a;
    double b;
    double burgers;
    double x1;
    double y1;
    double x2;
    double y2;
    double nu;
    int N;
    int dislocationType;
    string inputFile;
    string outputFile;

    cout << "Enter type of dislocation. 0 for Edge 1 for screw \n";
    cin >> dislocationType;

    cout << "Enter length of a vector: \n";
    cin >> a;

    cout << "Enter length of b vector: \n";
    cin >> b;

    cout << "Enter length of burgers vector: \n";
    cin >> burgers;

    cout << "Enter x location of dislocation 1: \n";
    cin >> x1;

    cout << "Enter y location of dislocation 1: \n";
    cin >> y1;

    cout << "Enter x location of dislocation 2: \n";
    cin >> x2;

    cout << "Enter y location of dislocation 2: \n";
    cin >> y2;

    cout << "Enter poisson ratio: \n";
    cin >> nu;

    cout << "Enter number of images: \n";
    cin >> N;

    cout << "Enter path to input file: \n";
    cin >> inputFile;

    cout << "Enter path to output file: \n";
    cin >> outputFile;

    displaceAtoms(dislocationType, inputFile, outputFile, a, b, burgers, x1, y1, x2, y2, nu, N);

    return 0;
}
