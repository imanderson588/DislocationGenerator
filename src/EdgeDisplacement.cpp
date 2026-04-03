#include "EdgeDisplacement.h"
#include <cmath>

static constexpr double pi = M_PI;

double edgeDisplacement_x(double x, double y, double x1, double y1, double x2, double y2, double burgers, double nu)
{
    double factor1 = burgers / (2 * pi);
    double theta1 = atan2(y - y1, -(x - x1));
    double term1 = ((x - x1) * (y - y1)) / (2 * (1 - nu) * (pow(x - x1, 2) + pow(y - y1, 2)));

    double factor2 = -burgers / (2 * pi);
    double theta2 = atan2(y - y2, -(x - x2));
    double term2 = ((x - x2) * (y - y2)) / (2 * (1 - nu) * (pow(x - x2, 2) + pow(y - y2, 2)));

    return (factor1 * (theta1 + term1)) + (factor2 * (theta2 + term2));
}

double edgeImage_x(double x, double y, double a, double b, double x1, double y1, double x2, double y2, double burgers, double nu, int N)
{
    double image_sum = 0;

    for (int i = -N; i <= N; i++)
    {
        for (int j = -N; j <= N; j++)
        {
            image_sum += edgeDisplacement_x(x - i * a, y - j * b, x1, y1, x2, y2, burgers, nu);
        }
    }
    return image_sum;
}

double edgeCorrection_x(double x, double y, double a, double b, double x1, double y1, double x2, double y2, double burgers, double nu, int N)
{
    double s_xx = edgeImage_x(a / 2, -b / 2, a, b, x1, y1, x2, y2, burgers, nu, N) - edgeImage_x(-a / 2, -b / 2, a, b, x1, y1, x2, y2, burgers, nu, N);
    double s_xy = edgeImage_x(-a / 2, b / 2, a, b, x1, y1, x2, y2, burgers, nu, N) - edgeImage_x(-a / 2, -b / 2, a, b, x1, y1, x2, y2, burgers, nu, N);
    return (s_xx / a) * x + (s_xy / b) * y;
}

double totEdge_x(double x, double y, double a, double b, double x1, double y1, double x2, double y2, double burgers, double nu, int N)
{
    return edgeImage_x(x, y, a, b, x1, y1, x2, y2, burgers, nu, N) - edgeCorrection_x(x, y, a, b, x1, y1, x2, y2, burgers, nu, N);
}

double edgeDisplacement_y(double x, double y, double x1, double y1, double x2, double y2, double burgers, double nu)
{
    double factor1 = -burgers / (2 * pi);
    double term11 = ((1 - 2 * nu) / (4 * (1 - nu))) * log(pow(x - x1, 2) + pow(y - y1, 2));
    double term12 = (pow(x - x1, 2) - pow(y - y1, 2)) / (4 * (1 - nu) * (pow(x - x1, 2) + pow(y - y1, 2)));

    double factor2 = burgers / (2 * pi);
    double term21 = ((1 - 2 * nu) / (4 * (1 - nu))) * log(pow(x - x2, 2) + pow(y - y2, 2));
    double term22 = (pow(x - x2, 2) - pow(y - y2, 2)) / (4 * (1 - nu) * (pow(x - x2, 2) + pow(y - y2, 2)));
    return (factor1 * (term11 + term12)) + (factor2 * (term21 + term22));
}

double edgeImage_y(double x, double y, double a, double b, double x1, double y1, double x2, double y2, double burgers, double nu, int N)
{
    double image_sum = 0;

    for (int i = -N; i <= N; i++)
    {
        for (int j = -N; j <= N; j++)
        {
            image_sum += edgeDisplacement_y(x - i * a, y - j * b, x1, y1, x2, y2, burgers, nu);
        }
    }
    return image_sum;
}

double edgeCorrection_y(double x, double y, double a, double b, double x1, double y1, double x2, double y2, double burgers, double nu, int N)
{
    double s_yx = edgeImage_y(a / 2, -b / 2, a, b, x1, y1, x2, y2, burgers, nu, N) - edgeImage_y(-a / 2, -b / 2, a, b, x1, y1, x2, y2, burgers, nu, N);
    double s_yy = edgeImage_y(-a / 2, b / 2, a, b, x1, y1, x2, y2, burgers, nu, N) - edgeImage_y(-a / 2, -b / 2, a, b, x1, y1, x2, y2, burgers, nu, N);
    return (s_yx / a) * x + (s_yy / b) * y;
}

double edgeDipoleTilt(double x, double y, double a, double b, double burgers, double loc1, double loc2, int N)
{
    double u_y = b * (y / b);
    return u_y;
}

double totEdge_y(double x, double y, double a, double b, double x1, double y1, double x2, double y2, double burgers, double nu, int N)
{
    return edgeImage_y(x, y, a, b, x1, y1, x2, y2, burgers, nu, N) - edgeCorrection_y(x, y, a, b, x1, y1, x2, y2, burgers, nu, N);
}
