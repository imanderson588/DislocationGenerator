#ifndef SCREWDISPLACEMENT_H
#define SCREWDISPLACEMENT_H

double screwDisplacement(double x, double y, double b);

double screwDipoleDisplacement(double x, double y, double burgers, double loc1, double loc2);

double screwDipoleImage(double x, double y, double a, double b, double burgers, double loc1, double loc2, int N);

double screwDipoleCorrection(double x, double y, double a, double b, double burgers, double loc1, double loc2, int N);

double screwDipoleTilt(double x, double y, double a, double b, double burgers, double loc1, double loc2, int N);

double screwDipole(double x, double y, double a, double b, double burgers, double loc1, double loc2, int N);

#endif
