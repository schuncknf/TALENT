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
double SphericalHOFunc::eval(int n, int l, double r) {
  return this->hoRadial(n,  l,  r);
}


//------------------------------------------------------------------------------
void SphericalHOFunc::setB(double j) {
    b_ = j;
}


//------------------------------------------------------------------------------
double SphericalHOFunc::getB(){
    return b_;
}


//------------------------------------------------------------------------------
double SphericalHOFunc::hoRadial (int n, int l, double r){
    double q = r / b_;
    double qsq = q * q;
    double a = (double)( l + 1. / 2.);

    return
      norm (n, l) * gsl_pow_int (q, l) * exp (-qsq / 2.)
      * gsl_sf_laguerre_n (n, a, qsq);
}


//------------------------------------------------------------------------------
double SphericalHOFunc::norm(int n, int l){
    int arg = 2*n + 2*l + 1;

    return sqrt (pow(2., n + l + 2.) * gsl_sf_fact((unsigned)n) /
             (pow(b_,3.) * sqrt(M_PI) * gsl_sf_doublefact((unsigned)arg)));
}
