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


using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




vector<string> splitBySpaces(const string& line) {
    
    istringstream iss(line);
    vector<string> tokens;
    string word;

    while (iss >> word) {
        tokens.push_back(word);
    }

    return tokens;
}

string recombine(const vector<string>& words) {
    ostringstream oss;
    for (size_t i = 0; i < words.size(); ++i) {
        oss << words[i];
        if (i != words.size() - 1) {
            oss << " ";  // Add space between words
        }
    }
    return oss.str();
}

void displaceAtoms(string &inputFile, string &outputFilePath, double a, double b, double burgers, double loc1, double loc2, int N){
    
    
    string lineContents;
    ofstream outputFile(outputFilePath);
    ifstream inputData(inputFile);
    int modifyFlag = 0;
    
    while (getline(inputData, lineContents)){
        
        if (lineContents=="Atoms # atomic"){
            outputFile << lineContents << "\n";
            modifyFlag = 1;
            continue;
        }
        
        if (lineContents=="Velocities"){
            outputFile << lineContents << "\n";
            modifyFlag = 0;
            continue;
        }
        
        if (lineContents==""){
            outputFile << lineContents << "\n";
            continue;
        }
        
        if (modifyFlag){

            vector<string> words = splitBySpaces(lineContents);

            double x_value = stod(words[2]);
            double y_value = stod(words[3]);
           
            double u_z = screwDipole(x_value, y_value, a, b, burgers, loc1, loc2, N);
            
            words[4] = to_string(stof(words[4])+u_z);

            string newLine = recombine(words);
            outputFile << newLine << "\n";
        }
        else {
            outputFile << lineContents << "\n";
        }
        
    }
    
    inputData.close();
    outputFile.close();
}



int main() {
    double a;
    double b;
    double burgers;
    double loc1;
    double loc2;
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
   cin >>  loc1;
   
   cout << "Enter location of dislocation 2: \n";
   cin >>  loc2;

   cout << "Enter number of images: \n";
   cin >>  N;

   cout << "Enter path to input file: \n";
   cin >>  inputFile;
 
   cout << "Enter path to output file: \n";
   cin >>  outputFile;

   displaceAtoms(inputFile, outputFile, a, b, burgers, loc1, loc2, N);

   cout << edgeDisplacement_y(0.5, 0.5, 1, -1, 1,0.25) << "\n";
    
    return 0;
}
