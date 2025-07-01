#include "SingleDislocations.h"
const double pi = M_PI;

double singleScrew(double x, double y, double burgers)
{

    double u_z = burgers * (atan2(y, x) / (2 * pi));
    return u_z;
}

double singleEdgeDisplacement_x(double x, double y, double burgers, double nu)
{

    double factor = burgers / (2 * pi);
    double theta = atan2(y, x);
    double term = (x * y) / (2 * (1 - nu) * (pow(x, 2) + pow(y, 2)));

    return (factor * (theta + term));
}

double singleEdgeImage_x(double x, double y, double a, double b, double burgers, double nu, int N)
{
    double image_sum = 0;
    double x_update = 0;
    double y_update = 0;

    for (double i = -N; i <= N; i++)
    {
        for (double j = -N; j <= N; j++)
        {
            vector<double> R = {i * a, j * b};
            x_update = x - R[0];
            y_update = y - R[1];
            image_sum += singleEdgeDisplacement_x(x_update, y_update, burgers, nu);
        }
    }
    return image_sum;
}

double singleEdgeCorrection_x(double x, double y, double a, double b, double burgers, double nu, int N)
{
    double s_xx = (singleEdgeImage_x(a / 2, y, a, b, burgers, nu, N) - singleEdgeImage_x(-a / 2, y, a, b, burgers, nu, N)) / a;
    double s_xy = (singleEdgeImage_x(a/2, -b / 2, a, b, burgers, nu, N) - singleEdgeImage_x(a/2, b / 2, a, b, burgers, nu, N));
    return  s_xy * y;
}

double totSingleEdge_x(double x, double y, double a, double b, double burgers, double nu, int N)
{

    return singleEdgeDisplacement_x(x, y, burgers, nu);// + singleEdgeCorrection_x(x, y, a, b, burgers, nu, N);
}

double singleEdgeDisplacement_y(double x, double y, double burgers, double nu)
{

    double factor1 = -burgers / (2 * pi);
    double term11 = ((1 - 2 * nu) / (4 * (1 - nu))) * log(pow(x, 2) + pow(y, 2));
    double term12 = (pow(x, 2) - pow(y, 2)) / (4 * (1 - nu) * (pow(x, 2) + pow(y, 2)));

    return (factor1 * (term11 + term12));
}

double singleEdgeImage_y(double x, double y, double a, double b, double burgers, double nu, int N)
{
    double image_sum = 0;
    double x_update = 0;
    double y_update = 0;

    for (double i = -N; i <= N; i++)
    {
        for (double j = -N; j <= N; j++)
        {
            vector<double> R = {i * a, j * b};
            x_update = x - R[0];
            y_update = y - R[1];
            image_sum += singleEdgeDisplacement_y(x_update, y_update, burgers, nu);
        }
    }
    return image_sum;
}

double singleEdgeCorrection_y(double x, double y, double a, double b, double burgers, double nu, int N)
{
    double s_yx = (singleEdgeImage_y(a / 2, y, a, b, burgers, nu, N) - singleEdgeImage_y(-a / 2, y, a, b, burgers, nu, N)) / a;
    double s_yy = (singleEdgeImage_y(x, b / 2, a, b, burgers, nu, N) - singleEdgeImage_y(x, -b / 2, a, b, burgers, nu, N)) / b;
    return s_yx * x + s_yy * y;
}

double totSingleEdge_y(double x, double y, double a, double b, double burgers, double nu, int N)
{
    double output = singleEdgeDisplacement_y(x, y, burgers, nu) ;//- singleEdgeImage_y(a/2, b/2, a, b, burgers, nu, N);
    return output;
}