#include "ParseLAMMPS.h"
#include "EdgeDisplacement.h"
#include "ScrewDisplacement.h"
#include "SingleDislocations.h"

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
            oss << " "; // Add space between words
        }
    }
    return oss.str();
}

void displaceAtoms(int dislocationType, string &inputFile, string &outputFilePath, double a, double b, double burgers, double x1, double y1, double x2, double y2, double nu, int N)
{

    string lineContents;
    ofstream outputFile(outputFilePath);
    ifstream inputData(inputFile);
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
                // Edge displacements
                double u_x = totEdge_x(x_value, y_value, a, b, x1, y1, x2, y2, burgers, nu, N);
                double u_y = totEdge_y(x_value, y_value, a, b, x1, y1, x2, y2, burgers, nu, N);
                // double u_y2 = edgeDipoleTilt(x_value, y_value, a, b, burgers, x1, x2, N);
                words[2] = to_string(stof(words[2]) + u_x); //+ 95.3025);
                words[3] = to_string(stof(words[3]) + u_y);
                // words[3] = to_string(stof(words[3]) + u_y2);
            }
            if (dislocationType == 1)
            {
                // Screw Displacements
                // words[2] = to_string(stof(words[2]) + u_z);

                double u_z = screwDipole(x_value, y_value, a, b, burgers, x1, x2, N);

                words[4] = to_string(stof(words[4]) + u_z);
            }

            if (dislocationType == 3)

            {
                // Single edge displacement
                // double u_x = totSingleEdge_x(x_value, y_value, a, b, burgers, nu, N);
                // double u_y = totSingleEdge_y(x_value, y_value, a, b, burgers, nu, N);
                if (y_value > 0)
                {
                    double u_x = 0;
                    words[2] = to_string(stof(words[2]) + u_x);
                }
                else
                {
                    double u_x = -1.575;
                    words[2] = to_string(stof(words[2]) + u_x);
                    words[2] = to_string(stof(words[2]) + (stof(words[2]) * 0.0123));
                }
                // words[3] = to_string(stof(words[3]) + u_y);
            }

            if (dislocationType == 4)

            {
                // Single edge displacement
                double u_z = singleScrew(x_value, y_value, burgers);

                words[4] = to_string(stof(words[4]) + u_z);
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
