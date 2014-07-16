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
double  SphericalHOFunc::hoRadial (int n, int l, double b, double r){
    double q = r / b;
    double qsq = q * q;
    double a = (double) l + 1. / 2.;

    return
      norm (n, l, b) * gsl_pow_int (q, (l + 1)) * exp (-qsq / 2.)
      * gsl_sf_laguerre_n ((n - 1), a, qsq);
}


//------------------------------------------------------------------------------
double SphericalHOFunc::norm(int n, int l, double b){
    double arg = (double) n + (double) l + 1. / 2.;

    return sqrt (2. * gsl_sf_fact ((unsigned) (n - 1)) /
             (b * gsl_sf_gamma (arg)));
}
