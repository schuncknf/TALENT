#ifndef SPHERICALHOFUNC_H
#define SPHERICALHOFUNC_H

#include <cmath>
#include <string>
#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_errno.h>

using namespace std;

class SphericalHOFunc
{
public:
    SphericalHOFunc();
    ~SphericalHOFunc();

    double eval(int n, int l, double b, double r);

private:
    double hoRadial (int n, int l, double b, double r);
    double norm (int n, int l, double b);
//    double hoEigenValue (int n, int l, double b, double m);

};

#endif // SPHERICALHOFUNC_H
