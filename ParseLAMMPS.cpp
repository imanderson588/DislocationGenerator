#include "ParseLAMMPS.h"
#include "EdgeDisplacement.h"
#include "ScrewDisplacement.h"

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

void displaceAtoms(string &inputFile, string &outputFilePath, double a, double b, double burgers, double x1, double y1, double x2, double y2, double nu, int N)
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

            double u_x = totEdge_x(x_value, y_value, a, b, x1, y1, x2, y2, burgers, nu, N);
            double u_y = totEdge_y(x_value, y_value, a, b, x1, y1, x2, y2, burgers, nu, N);

            words[2] = to_string(stof(words[2]) + u_x);
            words[3] = to_string(stof(words[3]) + u_y);

            // double u_z = screwDipole(x_value, y_value, a, b, burgers, x1, x2, N);

            // words[4] = to_string(stof(words[4]) + u_z);

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
