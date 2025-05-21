#include "EdgeDisplacement.h"
#include "ScrewDisplacement.h"

double edgeDisplacement_x(double x, double y, double loc1, double loc2, double burgers, double nu)
{

    double factor1 = burgers / (2 * pi);
    double theta1 = atan2(y, x - loc1);
    double term1 = ((x - loc1) * y) / (2 * (1 - nu) * (pow(x - loc1, 2) + pow(y, 2)));

    double factor2 = -burgers / (2 * pi);
    double theta2 = atan2(y, x - loc2);
    double term2 = ((x - loc2) * y) / (2 * (1 - nu) * (pow(x - loc2, 2) + pow(y, 2)));

    return (factor1 * (theta1 + term1)) + (factor2 * (theta2 + term2));
}

double edgeDisplacement_y(double x, double y, double loc1, double loc2, double burgers, double nu)
{

    double factor1 = -burgers / (2 * pi);
    double term11 = ((1 - 2 * nu) / (4 * (1 - nu))) * log(pow(x - loc1, 2) + pow(y, 2));
    double term12 = (pow(x - loc1, 2) - pow(y, 2)) / (4 * (1 - nu) * (pow(x - loc1, 2) + pow(y, 2)));
    double factor2 = burgers / (2 * pi);
    double term21 = ((1 - 2 * nu) / (4 * (1 - nu))) * log(pow(x - loc2, 2) + pow(y, 2));
    double term22 = (pow(x - loc2, 2) - pow(y, 2)) / (4 * (1 - nu) * (pow(x - loc2, 2) + pow(y, 2)));
    return (factor1 * (term11 + term12)) + (factor2 * (term21 + term22));
}
