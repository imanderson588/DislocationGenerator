#include "EdgeDisplacement.h"

double edgeDisplacement_x(double x, double y, double loc1, double loc2, double burgers, double nu){

    double factor1 = burgers/(2*pi);
    double theta1 = atan2(y, x-loc1);
    double term1 = ((x-loc1)*y)/(2*(1-nu)*(pow(x-loc1, 2)+pow(y,2)));

    double factor2 = -burgers/(2*pi);
    double theta2 = atan2(y, x-loc2);
    double term2 = ((x-loc2)*y)/(2*(1-nu)*(pow(x-loc2, 2)+pow(y,2)));

    return (factor1*(theta1 + term1))+(factor2*(theta2 + term2));

}

double edgeImage_x(double x, double y, double a, double b, double loc1, double loc2, double burgers, double nu, int N){


    double image_sum = 0;
    double x_update = 0;
    double y_update = 0;


    for (double i = -N; i <= N; i++){
        for(double j = -N; j <= N; j++){
            vector<double> R = {i*a, j*b};
            x_update = x - R[0];
            y_update = y - R[1];
            image_sum += edgeDisplacement_x(x_update, y_update, loc1, loc2, burgers, nu);



        }
    }
    return image_sum;

}

double edgeCorrection_x(double x, double y, double a, double b, double loc1, double loc2, double burgers, double nu, int N){

    double s_x = edgeImage_x(a/2, -b/2, a, b, loc1, loc2, burgers, nu, N) - edgeImage_x(-a/2, -b/2, a, b, loc1, loc2, burgers, nu, N);
    double s_y = edgeImage_x(-a/2, b/2, a, b, loc1, loc2, burgers, nu, N) - edgeImage_x(-a/2, -b/2, a, b, loc1, loc2, burgers, nu, N);
    double c_x = x;
    double c_y = y;
    double correction = (s_x*c_x) + (s_y*c_y);

    return correction;



}

double edgeDisplacement_y(double x, double y, double loc1, double loc2, double burgers, double nu){

    double factor1 = -burgers/(2*pi);
    double term11 = ((1-2*nu)/(4*(1-nu)))*log(pow(x-loc1, 2)+pow(y, 2));
    double term12 = (pow(x-loc1, 2)- pow(y, 2))/(4*(1-nu)*(pow(x-loc1, 2)+ pow(y, 2)));
    double factor2 = burgers/(2*pi);
    double term21 = ((1-2*nu)/(4*(1-nu)))*log(pow(x-loc2, 2)+pow(y, 2));
    double term22 = (pow(x-loc2, 2)-pow(y, 2))/(4*(1-nu)*(pow(x-loc2, 2)+pow(y, 2)));
    return (factor1*(term11+term12))+(factor2*(term21+term22));


}

double edgeImage_y(double x, double y, double a, double b, double loc1, double loc2, double burgers, double nu, int N){
    double image_sum = 0;
    double x_update = 0;
    double y_update = 0;


    for (double i = -N; i <= N; i++){
        for(double j = -N; j <= N; j++){
            vector<double> R = {i*a, j*b};
            x_update = x - R[0];
            y_update = y - R[1];
            image_sum += edgeDisplacement_y(x_update, y_update, loc1, loc2, burgers, nu);



        }
    }
    return image_sum;

}


double edgeCorrection_y(double x, double y, double a, double b, double loc1, double loc2, double burgers, double nu, int N){

    double s_x = edgeImage_y(a/2, -b/2, a, b, loc1, loc2, burgers, nu, N) - edgeImage_y(-a/2, -b/2, a, b, loc1, loc2, burgers, nu, N);
    double s_y = edgeImage_y(-a/2, b/2, a, b, loc1, loc2, burgers, nu, N) - edgeImage_y(-a/2, -b/2, a, b, loc1, loc2, burgers, nu, N);
    double c_x = x;
    double c_y = y;
    double correction = (s_x*c_x) + (s_y*c_y);

    return correction;

}
