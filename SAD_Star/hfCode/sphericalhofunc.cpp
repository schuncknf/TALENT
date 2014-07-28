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
  return hoRadial(n,  l,  r) * norm(n,l);
}

//------------------------------------------------------------------------------
double SphericalHOFunc::eval(int n, int l, double r, double& norm){
    return hoRadial(n, l, r) * norm;
}


//------------------------------------------------------------------------------
void SphericalHOFunc::setB(double j) {
    b_ = j;
    logb_= log(b_);
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
    if(l == 0 ){
        return exp(-qsq / 2.) * laguerrel0(n, qsq);
    }
    else{
        return pow(q, double(l)) * exp(-qsq / 2.) * gsl_sf_laguerre_n(n, l+0.5, qsq);
    }
}


//------------------------------------------------------------------------------
double SphericalHOFunc::norm(int n, int l){
    double logFacVal= logFac(n);

    int oddProd=1;
    for(int i=3; i< 2*n+2*l+2; i+=2){
        oddProd*=i;
    }
    double logOddProd= log(oddProd);

    double res= exp( logFacVal + LOG2*(2. + n + l) - LOGPI/2. - logOddProd -3.*logb_);

    return sqrt (res);
}


//------------------------------------------------------------------------------
double SphericalHOFunc::logFac(int n){
  int prod(1);
  for(int i=2; i<n+1; i++){
      prod*= i;
  }
    return log(prod);
}


//------------------------------------------------------------------------------
int SphericalHOFunc::fac(int n){
  int prod(1);
  for(int i=2; i<n+1; i++){
      prod*= i;
  }
  return prod;
}


//------------------------------------------------------------------------------
double SphericalHOFunc::laguerrel0(int n, double x){
    int pow2=2;
    for(int i=1; i<2*n+1; i++){
        pow2*=2;
    }
    double sqrtx= sqrt(x);
    return (1-n%2 *2)/(fac(n)* pow2* sqrtx)*  boost::math::hermite(2*n+1, sqrtx);
}


