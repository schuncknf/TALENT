#ifndef INTEGRATORGAUSSLEGENDRETEST_H
#define INTEGRATORGAUSSLEGENDRETEST_H

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include<iostream>
#include<iomanip>
#include "integratorGaussLegendre.h"


using namespace std;

//------------------------------------------------------------------------------
struct IntegratorGaussLegendreFixture{
    IntegratorGaussLegendre integrator_;
  
    IntegratorGaussLegendreFixture()
    {
        integrator_.readTables("lgvalues-weights.php", "lgvalues-abscissa.php");
    }
 
    ~IntegratorGaussLegendreFixture(){
    }
};

//------------------------------------------------------------------------------
//class Xexpx{
//public:
//    eval(double x, int dummy){
//        return x*exp(x);
//    }
//}


BOOST_FIXTURE_TEST_SUITE( integratorGaussLegendre, IntegratorGaussLegendreFixture )


//------------------------------------------------------------------------------
double funcTest0(double x){
    return x;
}
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( func0Test )
{
    integrator_.setOrder(2);
    double res0= integrator_.integrate(funcTest0, 0., 10);
    BOOST_CHECK_CLOSE(res0, 50., 1e-6);
}



//------------------------------------------------------------------------------
double funcTest1(double x){
    return x*x;
}
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( func1Test )
{
    integrator_.setOrder(2);
    double res= integrator_.integrate(funcTest1, -10, 10);
    BOOST_CHECK_CLOSE(res, 2.e3/3., 1e-6);
}

//------------------------------------------------------------------------------
double funcTest2(double x){
    return x*exp(-x);
}
//------------------------------------------------------------------------------
double funcTest2Prim(double x, int dummy){
    return x*exp(-x);
}
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( func2Test )
{
    integrator_.setOrder(60);

    double res= integrator_.integrate(funcTest2, 0., 1e3);
    BOOST_CHECK_CLOSE(res, 1., 1e-2);


    res= integrator_.integrate0ToInf(funcTest2);
    BOOST_CHECK_CLOSE(res, 1., 1e-6);

    auto integrand= bind(funcTest2Prim, placeholders::_1, 0);
    res= integrator_.integrate0ToInf(integrand);
}


//------------------------------------------------------------------------------
double funcTest3(double x){
    return 1./x;
}
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( func3Test )
{
    integrator_.setOrder(60);

    double res= integrator_.integrate(funcTest3, 1., 1e3);
    BOOST_CHECK_CLOSE(res,  6.90775527898, 1e-1);
}




BOOST_AUTO_TEST_SUITE_END()



#endif
