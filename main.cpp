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

using namespace std;

int main()
{
    double a;
    double b;
    double burgers;
    double loc1;
    double loc2;
    double nu;
    int N;
    string inputFile;
    string outputFile;

    cout << "Enter length of a vector: \n";
    cin >> a;

    cout << "Enter length of b vector: \n";
    cin >> b;

    cout << "Enter length of burgers vector: \n";
    cin >> burgers;

    cout << "Enter location of dislocation 1: \n";
    cin >> loc1;

    cout << "Enter location of dislocation 2: \n";
    cin >> loc2;

    cout << "Enter poisson ratio: \n";
    cin >> nu;

    cout << "Enter number of images: \n";
    cin >> N;

    cout << "Enter path to input file: \n";
    cin >> inputFile;

    cout << "Enter path to output file: \n";
    cin >> outputFile;

    displaceAtomsVasp(inputFile, outputFile, a, b, burgers, loc1, loc2, nu, N);

    return 0;
}
