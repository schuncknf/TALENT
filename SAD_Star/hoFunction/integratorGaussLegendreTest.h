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


BOOST_AUTO_TEST_SUITE_END()



#endif
