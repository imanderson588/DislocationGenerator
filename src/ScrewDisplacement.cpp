#include "ScrewDisplacement.h"
#include <cmath>

static constexpr double pi = M_PI;

double screwDisplacement(double x, double y, double b)
{
    double theta = atan2(y, x);
    double u_z = b * (theta / (2 * pi));

    return u_z;
}

double screwDipoleDisplacement(double x, double y, double burgers, double loc1, double loc2)
{
    double theta1 = atan2(y, (x - loc1));
    double theta2 = atan2(y, (x - loc2));
    double frac = (theta1 - theta2) / (2 * pi);

    double displacement = burgers * frac;

    return displacement;
}

double screwDipoleImage(double x, double y, double a, double b, double burgers, double loc1, double loc2, int N)
{
    double image_sum = 0;

    for (int i = -N; i <= N; i++)
    {
        for (int j = -N; j <= N; j++)
        {
            image_sum += screwDipoleDisplacement(x - i * a, y - j * b, burgers, loc1, loc2);
        }
    }
    return image_sum;
}

double screwDipoleCorrection(double x, double y, double a, double b, double burgers, double loc1, double loc2, int N)
{
    double s_x = screwDipoleImage(a / 2, -b / 2, a, b, burgers, loc1, loc2, N) - screwDipoleImage(-a / 2, -b / 2, a, b, burgers, loc1, loc2, N);
    double s_y = screwDipoleImage(-a / 2, b / 2, a, b, burgers, loc1, loc2, N) - screwDipoleImage(-a / 2, -b / 2, a, b, burgers, loc1, loc2, N);
    double correction = (s_x / a) * x + (s_y / b) * y;

    return correction;
}

double screwDipoleTilt(double x, double y, double a, double b, double burgers, double loc1, double loc2, int N)
{
    double u_z = 0.5 * burgers * y / b;
    return u_z;
}

double screwDipole(double x, double y, double a, double b, double burgers, double loc1, double loc2, int N)
{
    return screwDipoleImage(x, y, a, b, burgers, loc1, loc2, N) - screwDipoleCorrection(x, y, a, b, burgers, loc1, loc2, N) - screwDipoleTilt(x, y, a, b, burgers, loc1, loc2, N);
}
