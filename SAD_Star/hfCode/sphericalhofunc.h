#ifndef SPHERICALHOFUNC_H
#define SPHERICALHOFUNC_H

#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>

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

    double eval(int n, int l, double r);
    void setB(double j);
    double getB();

private:
    double hoRadial (int n, int l, double r);
    double norm (int n, int l);
    double b_;
    double logFac(int n);
    double logDoubleFac(int n);
};

#endif // SPHERICALHOFUNC_H
