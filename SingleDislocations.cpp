#include "SingleDislocations.h"
const double pi = M_PI;

double singleEdgeDisplacement_x(double x, double y, double burgers, double nu)
{

    double factor = burgers / (2 * pi);
    double theta = atan2(y, x);
    double term = (x * y) / (2 * (1 - nu) * (pow(x, 2) + pow(y, 2)));

    return (factor * (theta + term));
}

double singleEdgeDisplacement_y(double x, double y, double burgers, double nu)
{

    double factor1 = -burgers / (2 * pi);
    double term11 = ((1 - 2 * nu) / (4 * (1 - nu))) * log(pow(x, 2) + pow(y, 2));
    double term12 = (pow(x, 2) - pow(y, 2)) / (4 * (1 - nu) * (pow(x, 2) + pow(y, 2)));

    return (factor1 * (term11 + term12));
}