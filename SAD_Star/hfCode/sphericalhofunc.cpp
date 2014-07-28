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

//    double logOddProd=0;
//    for(int i=3; i< 2*n+2*l+2; i+=2){
//        logOddProd+= log(i);
//    }

    double log2Pow= LOG2 * (2. + n + l);

    double logPiSqr= LOGPI/2.;

    double res= exp( logFacVal + log2Pow - logPiSqr - logOddProd -3.*logb_);

    return sqrt (res);
}
//------------------------------------------------------------------------------
double SphericalHOFunc::logFac(int n){
  int prod(1);
  for(int i=2; i<n+1; i++){
      prod*= i;
  }
    return log(prod);
//    double sum(0);
//    for(int i=2; i<n+1; i++){
//        sum+= log(i);
//    }
//    return sum;
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
//    return gsl_sf_laguerre_n(n, l+0.5, qsq);
    //cout<<gsl_sf_laguerre_n(n, 0.5, x)<<"  "<<double(1-2*n%2)/(fac(n)*pow(2,2*n+1)*sqrt(x))* boost::math::hermite(2*n+1, sqrt(x))<<endl;
    //cout<<"n"<<n<<" fac"<<fac(n)<<" s"<<double(1-2*n%2)<<"po "<<pow(2,n+1)<<endl;

    double res= double(1-n%2 *2)/(fac(n)*pow(2,2*n+1)*sqrt(x))* boost::math::hermite(2*n+1, sqrt(x));
//    if(res != gsl_sf_laguerre_n(n, 0.5, x)){
//        cout<<n<<"  "<<res<<"  "<<gsl_sf_laguerre_n(n, 0.5, x)<<endl;
//    }
    return res;
}


////------------------------------------------------------------------------------
//double SphericalHOFunc::logDoubleFac(int n){
//    double res(0);
//    int k= n/2;
//    int i0=1;
//    if(2*k == n){
//        i0=2;
//    }
//      for(int i=i0; i<n+1; i+=2){
//        res+= log(i);
//    }
//      return res;
//}
