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
  return hoRadial(n,  l,  r);
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
// Before modif
double SphericalHOFunc::hoRadial (int n, int l, double r){
    double q = r / b_;
    double qsq = q * q;
    return norm (n, l) * pow(q, double(l)) * exp(-qsq / 2.) * gsl_sf_laguerre_n(n, l+0.5, qsq);
}


//------------------------------------------------------------------------------
double SphericalHOFunc::norm(int n, int l){
    double logFacVal= logFac(n);

    double logOddProd=0;
    for(int i=1; i< 2*n+2*l+2; i+=2){
        logOddProd+=log(i);
    }

    double log2Pow= log(2.)*(n + l + 2.);

    double logPiSqr= 0.5* log(M_PI);

    double res= exp( logFacVal +  log2Pow - logPiSqr - logOddProd -3.*log(b_));

    return sqrt (res);
}
//------------------------------------------------------------------------------
double SphericalHOFunc::logFac(int n){
  double res(0);
  for(int i=1; i<n+1; i++){
      res+= log(i);
  }
    return res;
}


//------------------------------------------------------------------------------
double SphericalHOFunc::logDoubleFac(int n){
    double res(0);
    int k= n/2;
    int i0=1;
    if(2*k == n){
        i0=2;
    }
      for(int i=i0; i<n+1; i+=2){
        res+= log(i);
    }
      return res;
}
