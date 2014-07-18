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

////------------------------------------------------------------------------------
//// After modif
//double SphericalHOFunc::hoRadial (int n, int l, double r){
//    double q = r * b_;
//    double qsq = q * q;
//    double a = (double)( l + 1./2.);

////    cout<<r<<"  norm " <<norm (n, l)<<"   q^l "<<pow(q,double(l))<< "  exp(-qsi/2) "<<exp(-qsq / 2.)<<"  "<<gsl_sf_laguerre_n (n, a, qsq)<<endl;

//    return norm (n, l) * pow(q,double(l)) * exp(-qsq / 2.) * gsl_sf_laguerre_n (n, a, qsq);
//}

////------------------------------------------------------------------------------
//double SphericalHOFunc::norm(int n, int l){

//    double logOddProd=0;
//    for(int i=1; i< 2*n+2*l+2; i+=2){
//        logOddProd+=log(i);
//    }

//    double logC= myLogFac(n) + (n+l+1)*log(2) - logOddProd;

////    cout<<n<<" "<<l<<" "<<b_<<" "<<sqrt(2.*b_*b_*b_)/pow(M_PI, 0.25)* exp(logC/2.)<<endl;
////    cout<<"logOddProd "<<logOddProd<<endl;
//    return sqrt(2.*b_*b_*b_)/pow(M_PI, 0.25)* exp(logC/2.);
//}



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


////------------------------------------------------------------------------------
//double radiallogs(int n, int l, double r, double nu)
//{
//  // v = nu = m*omega/(2*h_bar) this is different from practioners!
//  double logcoeff, R,coeff;
//  //v = 4.0;  // this is the nu from wikipedia.
//  // v = 0.5 will diagonalize the HO potential

//  //logcoeff = 0.5*(1.5*log(v) + (l-n+2)*log(2) + logoddfact(2*n+2*l+1))
//   // -0.5*(0.5*log(M_PI) + logfact(n) + 2*logoddfact(2*l+1));


//  logcoeff = 0.5*((1.5+l)*log(nu) + (n+2*l+3.5)*log(2) + logfact(n)
//		  - 0.5*log(M_PI) - logoddfact(2*n+2*l+1));



//  R = exp(logcoeff)*pow(r,l)*recursive(n,l+0.5,2.0*nu*r*r)*exp(-nu*r*r);

//  return R;
//}
