#include "ParseVASP.h"
#include "ParseLAMMPS.h"
#include "EdgeDisplacement.h"
#include "ScrewDisplacement.h"

void displaceAtomsVasp(string &inputFile, string &outputFilePath, double a, double b, double burgers, double x1, double y1, double x2, double y2, double nu, int N)
{

    string lineContents;
    ofstream outputFile(outputFilePath);
    ifstream inputData(inputFile);
    int lineNumber = 0;
    double scaleFactor = 0;
    double xLength = 0;
    double yLength = 0;
    double zLength = 0;
    int numAtoms = 0;

    while (getline(inputData, lineContents))
    {

        if (lineNumber == 0)
        {
            outputFile << lineContents << "\n";
            // modifyFlag = 1;
            lineNumber++;
            continue;
        }

        if (lineNumber == 1)
        {
            vector<string> words = splitBySpaces(lineContents);
            scaleFactor = stod(words[0]);
            outputFile << lineContents << "\n";
            lineNumber++;
            continue;
        }

        if (lineNumber == 2 || lineNumber == 3 || lineNumber == 4)
        {
            vector<string> words = splitBySpaces(lineContents);

            if (lineNumber == 2)
            {
                xLength = scaleFactor * stod(words[0]);
            }

            if (lineNumber == 3)
            {
                yLength = scaleFactor * stod(words[1]);
            }

            if (lineNumber == 4)
            {
                zLength = scaleFactor * stod(words[2]);
            }

            outputFile << lineContents << "\n";
            lineNumber++;
            continue;
        }

        if (lineNumber == 6)
        {
            vector<string> words = splitBySpaces(lineContents);
            numAtoms = stoi(words[0]); //+ stoi(words[1]);
            outputFile
                << lineContents << "\n";
            lineNumber++;
            continue;
        }

        if (lineNumber > 7 && lineNumber <= (7 + numAtoms))
        {

            vector<string> words = splitBySpaces(lineContents);

            double x_value = stod(words[0]) * xLength;
            double y_value = stod(words[1]) * yLength;

            double u_x = totEdge_x(x_value, y_value, a, b, x1, y1, x2, y2, burgers, nu, N);
            double u_y = totEdge_y(x_value, y_value, a, b, x1, y1, x2, y2, burgers, nu, N);

            words[0] = to_string((x_value + u_x) / xLength);
            words[1] = to_string((y_value + u_y) / yLength);

            string newLine = recombine(words);
            outputFile << newLine << "\n";
            lineNumber++;
            continue;
        }
        else
        {
            outputFile << lineContents << "\n";
            lineNumber++;
            continue;
        }
    }

    inputData.close();
    outputFile.close();
}
