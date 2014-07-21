#ifndef IntegratorGaussLaguerreTEST_H
#define IntegratorGaussLaguerreTEST_H

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include<iostream>
#include<iomanip>

#include "integratorGaussLaguerre.h"


using namespace std;

//------------------------------------------------------------------------------
struct IntegratorGaussLaguerreFixture{
    IntegratorGaussLaguerre integrator_;
  
    IntegratorGaussLaguerreFixture()
    {
        integrator_.readTables("../gen_laguerre/gauss-laguerre_n20_w.txt", "../gen_laguerre/gauss-laguerre_n20_x.txt", 20);
        integrator_.readTables("../gen_laguerre/gauss-laguerre_n100_w.txt", "../gen_laguerre/gauss-laguerre_n100_x.txt", 100);
        integrator_.readTables("../gen_laguerre/gauss-laguerre_n150_w.txt", "../gen_laguerre/gauss-laguerre_n150_x.txt", 150);
        integrator_.readTables("../gen_laguerre/gauss-laguerre_n200_w.txt", "../gen_laguerre/gauss-laguerre_n200_x.txt", 200);
        integrator_.readTables("../gen_laguerre/gauss-laguerre_n500_w.txt", "../gen_laguerre/gauss-laguerre_n500_x.txt", 500);
    }
 
    ~IntegratorGaussLaguerreFixture(){
    }
};


BOOST_FIXTURE_TEST_SUITE( integratorGaussLaguerre, IntegratorGaussLaguerreFixture )


//------------------------------------------------------------------------------
double funcTest2(double x){
    return x*exp(-x);
}
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( func2Test )
{
    integrator_.setOrder(20);
    double res= integrator_.integrate0ToInf(funcTest2);
    BOOST_CHECK_CLOSE(res, 1., 1e-6);
}

//------------------------------------------------------------------------------
double funcTest3(double x){
    return sin(x)/x;
}
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( func3Test )
{
    integrator_.setOrder(500);
    double res= integrator_.integrate0ToInf(funcTest3);
    BOOST_CHECK_CLOSE(res, M_PI/2., 0.1);
}




BOOST_AUTO_TEST_SUITE_END()



#endif
