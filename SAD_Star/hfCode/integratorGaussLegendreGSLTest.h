#ifndef integratorGaussLegendreGSLTEST_H
#define integratorGaussLegendreGSLTEST_H

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include<iostream>
#include<iomanip>

#include "integratorGaussLegendreGSL.h"


using namespace std;

//------------------------------------------------------------------------------
struct integratorGaussLegendreGSLFixture{
    IntegratorGaussLegendreGSL integrator_;
  
    integratorGaussLegendreGSLFixture()
    {
    }
 
    ~integratorGaussLegendreGSLFixture(){
    }
};


BOOST_FIXTURE_TEST_SUITE( integratorGaussLegendreGSL, integratorGaussLegendreGSLFixture )


//------------------------------------------------------------------------------
double funcTest0(double x, void* param){
    return x;
}
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( func0Test )
{
    gsl_function F;
    F.function = funcTest0;
    F.params = NULL;

    double res0= integrator_.integrate(F, 0., 10);
    BOOST_CHECK_CLOSE(res0, 50., 1e-6);
}


//------------------------------------------------------------------------------
double funcTest1(double x, void*){
    return x*x;
}
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( func1Test )
{
    gsl_function F;
    F.function = funcTest1;
    F.params = NULL;
    double res= integrator_.integrate(F, -10, 10);
    BOOST_CHECK_CLOSE(res, 2.e3/3., 1e-6);
}

////------------------------------------------------------------------------------
//double funcTest2(double x, void*){
//    return x*exp(-x);
//}
////------------------------------------------------------------------------------
//BOOST_AUTO_TEST_CASE( func2Test )
//{
//    gsl_function F;
//    F.function = funcTest0;
//    F.params = NULL;
//    double res= integrator_.integrate0ToInf(F);
//    BOOST_CHECK_CLOSE(res, 1., 1e-6);
//}



BOOST_AUTO_TEST_SUITE_END()



#endif
