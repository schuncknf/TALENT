#include "sphericalhofunc.h"


//------------------------------------------------------------------------------
SphericalHOFunc::SphericalHOFunc()
{
}


//------------------------------------------------------------------------------
SphericalHOFunc::~SphericalHOFunc()
{
}


//------------------------------------------------------------------------------
double SphericalHOFunc::eval(int n, int l, double b, double r){
  return this->hoRadial( n,  l,  b,  r);
}


//------------------------------------------------------------------------------
double SphericalHOFunc::hoRadial (int n, int l, double b, double r){
    double q = r / b;
    double qsq = q * q;
    double a = (double)( l + 1. / 2.);

    return
      norm (n, l, b) * gsl_pow_int (q, l) * exp (-qsq / 2.)
      * gsl_sf_laguerre_n (n, a, qsq);
}


//------------------------------------------------------------------------------
double SphericalHOFunc::norm(int n, int l, double b){
    int arg = 2*n + 2*l + 1;

    return sqrt (pow(2., n + l + 2.) * gsl_sf_fact((unsigned)n) /
             (b * sqrt(M_PI) * gsl_sf_doublefact((unsigned)arg)));
}
