#ifndef SPHERICALHOFUNCTEST_H
#define SPHERICALHOFUNCTEST_H


#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include<iostream>
#include<iomanip>

#include <functional>

#include "sphericalhofunc.h"
#include "integratorGaussLegendre.h"
#include "integratorGaussLaguerre.h"
#include "constants.h"

using namespace std;

//------------------------------------------------------------------------------
struct SphericalHOFuncFixture{

    SphericalHOFuncFixture()
    {
    }

    ~SphericalHOFuncFixture(){
    }
};


BOOST_FIXTURE_TEST_SUITE( sphericalHOFunc, SphericalHOFuncFixture )


//------------------------------------------------------------------------------
double hoProdTest(double r, SphericalHOFunc& func, int n1, int n2){
       return func.eval(n1, 0, r) * func.eval(n2, 0, r) * r * r;
}

//------------------------------------------------------------------------------
//BOOST_AUTO_TEST_CASE( orthonormalTest )
//{
//    IntegratorGaussLegendre integrator;
//    integrator.readTables("lgvalues-weights.php", "lgvalues-abscissa.php");
//    integrator.setOrder(60);

//    SphericalHOFunc funcR;
//    funcR.setB(0.5);


//    auto integrand= bind(hoProdTest, placeholders::_1, funcR, 19, 19);
//    cout<<hoProdTest(0.5, funcR, 19, 19)<<"  "<<integrand(0.5)<<endl;
//    cout<<hoProdTest(2., funcR, 19, 19)<<"  "<<integrand(2.)<<endl;


//    double cumulError= 0;
//    int nMax= 50;
//    for(int n1=nMax-nMax; n1<nMax; n1++){
//        for(int n2=0; n2<nMax; n2++){
//            auto integrand= bind(hoProdTest, placeholders::_1, funcR, n1, n2);
//            double res= integrator.integrate0ToInf(integrand);
//            double error= abs(res - int(n1==n2));
//            if(n1==n2){
//                cout<<n1<<"  " <<n2<<"  "<<res<<endl;
//            }
//            cumulError+= error;
//        }
//    }

//    BOOST_CHECK_CLOSE(cumulError +1., 1., 1e-6);
//}


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( orthonormal2Test )
{
    IntegratorGaussLaguerre integrator;
    integrator.setTableDir("../gen_laguerre");
    //integrator.readTables( );
    integrator.setOrder(500);

    SphericalHOFunc funcR;
    funcR.setB(0.5);


    double cumulError= 0;

    int nMax= 4;
    for(int n1=0; n1<nMax; n1++){
        for(int n2= 0; n2< nMax; n2++){
            auto integrand= bind(hoProdTest, placeholders::_1, funcR, n1, n2);
            double res= integrator.integrate0ToInf(integrand);
            double error= abs(res - double(n1==n2));
//            if(n1==n2){
//                cout<<n1<<"  " <<n2<<"  "<<res<<endl;
//            }
            cumulError+= error;
        }
    }

    BOOST_CHECK_CLOSE(cumulError +1., 1., 1e-12);
}


////------------------------------------------------------------------------------
//BOOST_AUTO_TEST_CASE( plotTest ){
//    SphericalHOFunc funcR;
//    funcR.setB(0.5);

//    ofstream output("r19.dat");
//    for(double r=0; r<40; r+=0.01){
//        output<<r <<"   "<<setprecision(10)<<funcR.eval(19,0, r)<<endl;
//    }
//    output.close();
//}


//------------------------------------------------------------------------------
// Anton function:
double sho_wf(double r, double mw, int n, int l)
{
  int i, npl;
  double x, mwr2, cur, last, last2, factor;
  if ((r < 0) || (mw < 0) || (n < 0) || (l < 0))
    return 0.;
  if (r == 0) {
    if (l != 0)
      return 0.;
    factor = mw;
    for (i = 1; i <= n; i++)
      factor *= (i + 0.5) / i;
    return 2 * sqrt(sqrt(mw / M_PI) * factor);
  }
  x = sqrt(mw) * r;
  mwr2 = x * x;
  npl = n + l;
  factor = mw;
  for (i = 3; i <= 2 * npl + 1; i += 2)
    factor *= 2 * mwr2 / i;
  cur = 0.;
  last2 = sqrt(factor);
  if (l == npl)
    cur = last2;
  last = ((npl + 0.5) / x - x) * last2;
  if (l == npl - 1)
    cur = last;
  for (i = npl - 2; i >= l; i--) {
    cur = (((i + 1.5) / x - x) * last - sqrt(npl - i - 1) * last2) / sqrt(npl - i);
    last2 = last;
    last = cur;
  }
  return 2 * pow(mw / M_PI, 0.25) * exp(-0.5*mwr2) * cur;
}

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( antonBenchmark ){
    double hbarOmega= 10.;
    double b= HBARC/sqrt(MNC2*hbarOmega);
    double mw= MNC2* hbarOmega/HBARC/HBARC;

    SphericalHOFunc funcR;
    funcR.setB(b);

    double sum=0;
    for(int n=0; n<10; n++){
        int l=0;
        for(double r=0.1 ; r<100*b; r+=0.1){
            double RnBench= sho_wf(r, mw, n, l);
            double Rn= funcR.eval(n, l, r);
            double diff= abs(Rn- RnBench);
            cout<<Rn<<"  "<<RnBench<<endl;
            sum+= diff;
        }
    }

    BOOST_CHECK_LT(sum , 1e-13);
}


BOOST_AUTO_TEST_SUITE_END()

#endif // SPHERICALHOFUNCTEST_H
