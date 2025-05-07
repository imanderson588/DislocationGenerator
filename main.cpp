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


using namespace std;

const double pi = M_PI;


//Functions for screw dislocation dipole //////////////////////////////////////////////////////////////////

double screwDipoleDisplacement(double x, double y, double burgers, double loc1, double loc2){
    
    double theta1 = atan2(y, (x-loc1));
    double theta2 = atan2(y, (x-loc2));
    double frac = (theta1 - theta2)/(2*pi);

    double displacement = burgers*frac;

    return displacement;

}


double screwDipoleImage(double x, double y, double a, double b, double burgers, double loc1, double loc2, int N){

    double image_sum = 0;
    double x_update = 0;
    double y_update = 0;


    for (double i = -N; i <= N; i++){
        for(double j = -N; j <= N; j++){
            vector<double> R = {i*a, j*b};
            x_update = x - R[0];
            y_update = y - R[1];
            image_sum += screwDipoleDisplacement(x_update, y_update, burgers, loc1, loc2);



        }
    }
    return image_sum;

}

double screwDipoleCorrection(double x, double y, double a, double b, double burgers, double loc1, double loc2, int N){

    double s_x = screwDipoleImage(a/2, -b/2, a, b, burgers, loc1, loc2, N) - screwDipoleImage(-a/2, -b/2, a, b, burgers, loc1, loc2, N);
    double s_y = screwDipoleImage(-a/2, b/2, a, b, burgers, loc1, loc2, N) - screwDipoleImage(-a/2, -b/2, a, b, burgers, loc1, loc2, N);
    double c_x = x;
    double c_y = y;
    double correction = (s_x*c_x) + (s_y*c_y);

    return correction;


}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



double screwDisplacement(double x, double y, double b){

    double theta = atan2(y,x);
    double u_z = b*(theta/(2*pi));
    
    return u_z;
}


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

void displaceAtoms(const string& inputFile){
    
    
    string lineContents;
    ofstream outputFile("");
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
            double u_z = screwDisplacement(x_value, y_value, 5.718);

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
    double x = 2;
    double y = 1;
    double a = 2;
    double b = 1;
    double burgers = 1;
    double loc1 = 0.5;
    double loc2 = -0.5;
    int N = 10;
    
    double u_z  = screwDisplacement(x, y, b);
    
    cout << "Displacement in z-direction: " << u_z << "\n";
    
    double u_z_inf = screwDipoleCorrection(x, y, a, b, burgers, loc1, loc2, N);
    cout << u_z_inf << "\n";

    displaceAtoms("");
    
    return 0;
}
