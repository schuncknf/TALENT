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

#include <boost/math/special_functions/laguerre.hpp>
#include <boost/math/special_functions/hermite.hpp>


#include "constants.h"

using namespace std;

class SphericalHOFunc
{
public:
    SphericalHOFunc();
//    SphericalHOFunc( SphericalHOFunc& orig);
    ~SphericalHOFunc();

    double eval(int n, int l, double r);
    double eval(int n, int l, double r, double& norm);

    void setB(double j);
    double getB();
    double norm (int n, int l);

private:
    double hoRadial (int n, int l, double r);

    double b_;
    double logb_;
    double logFac(int n);
    int fac(int n);
    double laguerrel0(int n, double x);



//    double logDoubleFac(int n);
};

#endif // SPHERICALHOFUNC_H
