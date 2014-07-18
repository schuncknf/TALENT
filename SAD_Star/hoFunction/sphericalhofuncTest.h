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
    return func.eval(n1, 0, r) * func.eval(n2, 0, r) * r*r;
}

//------------------------------------------------------------------------------
double hoProd2Test(double r, SphericalHOFunc& func, int n1, int n2){
    return func.eval(n1, 0, r) * func.eval(n2, 0, r) *r *r;
}

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( orthonormalTest )
{
    IntegratorGaussLegendre integrator;
    integrator.readTables("lgvalues-weights.php", "lgvalues-abscissa.php");
    integrator.setOrder(60);
    double xMin=1e-6;
    double xMax= 1e3;

    SphericalHOFunc funcR;
    funcR.setB(0.05);


    double cumulError= 0;

    int nMax= 20;
    for(int n1=0; n1<nMax; n1++){
        for(int n2= 0; n2< nMax; n2++){
            auto integrand= bind(hoProdTest, placeholders::_1, funcR, n1, n2);
            double res= integrator.integrate0ToInf(integrand);
            double error= abs(res - int(n1==n2));
            if(n1==n2){
                cout<<n1<<"  " <<n2<<"  "<<error<<endl;
            }
            cumulError+= error;
        }
    }

    BOOST_CHECK_CLOSE(cumulError +1., 1., 1e-6);
}


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( orthonormal2Test )
{
    IntegratorGaussLaguerre integrator;
    integrator.readTables( "../gen_laguerre/gauss-laguerre_n100_w.txt", "../gen_laguerre/gauss-laguerre_n100_x.txt", 100);
    integrator.setOrder(100);

    SphericalHOFunc funcR;
    funcR.setB(0.2);


    double cumulError= 0;

    int nMax= 20;
    for(int n1=0; n1<nMax; n1++){
        for(int n2= 0; n2< nMax; n2++){
            auto integrand= bind(hoProdTest, placeholders::_1, funcR, n1, n2);
            double res= integrator.integrate0ToInf(integrand);
            double error= abs(res - double(n1==n2));
            if(n1==n2){
                cout<<n1<<"  " <<n2<<"  "<<error<<endl;
            }

//            cout<<n1<<"  " <<n2<<"  "<<res<<endl;
            cumulError+= error;
        }
    }

    BOOST_CHECK_CLOSE(cumulError +1., 1., 1e-6);
}


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( plotTest ){
    SphericalHOFunc funcR;
    funcR.setB(1.);

    ofstream output("r19.dat");
    for(double r=0; r<10; r+=0.1){
        output<<funcR.eval(19,0, r)<<endl;
    }
}


BOOST_AUTO_TEST_SUITE_END()

#endif // SPHERICALHOFUNCTEST_H
