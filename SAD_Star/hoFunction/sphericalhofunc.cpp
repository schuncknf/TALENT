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
//double SphericalHOFunc::hoRadial (int n, int l, double r){
//    double q = r / b_;
//    double qsq = q * q;
//    double a = (double)( l + 1./ 2.);

//    return
//      norm (n, l) * gsl_pow_int(q, l) * exp (-qsq / 2.)
//      * gsl_sf_laguerre_n (n, a, qsq);
//}


////------------------------------------------------------------------------------
//double SphericalHOFunc::norm(int n, int l){
//    int arg = 2*n + 2*l + 1;
//    //double logFac= log(gsl_sf_fact((unsigned)n));
//    double logFacVal= this->myLogFac(n);

//    //double logFacFac= log(gsl_sf_doublefact((unsigned)arg));
//    double logFacFac= logDoubleFac(arg);

//    double log2Pow= log(pow(2., n + l + 2.));

//    double prod= logFacVal+ log2Pow - logFacFac;
//    double res= exp(prod)/ pow(b_,3.) / sqrt(M_PI);

//    return sqrt (res);
//}

//------------------------------------------------------------------------------
// After modif
double SphericalHOFunc::hoRadial (int n, int l, double r){
    double q = r * b_;
    double qsq = q * q;
    double a = (double)( l + 1./2.);

//    cout<<r<<"  norm " <<norm (n, l)<<"   q^l "<<pow(q,double(l))<< "  exp(-qsi/2) "<<exp(-qsq / 2.)<<"  "<<gsl_sf_laguerre_n (n, a, qsq)<<endl;

    return norm (n, l) * pow(q,double(l)) * exp(-qsq / 2.) * gsl_sf_laguerre_n (n, a, qsq);
}

//------------------------------------------------------------------------------
double SphericalHOFunc::norm(int n, int l){

    double logOddProd=0;
    for(int i=1; i< 2*n+2*l+2; i+=2){
        logOddProd+=log(i);
    }

    double logC= myLogFac(n) + (n+l+1)*log(2) - logOddProd;

//    cout<<n<<" "<<l<<" "<<b_<<" "<<sqrt(2.*b_*b_*b_)/pow(M_PI, 0.25)* exp(logC/2.)<<endl;
//    cout<<"logOddProd "<<logOddProd<<endl;
    return sqrt(2.*b_*b_*b_)/pow(M_PI, 0.25)* exp(logC/2.);
}



//------------------------------------------------------------------------------
double SphericalHOFunc::myLogFac(int n){
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
