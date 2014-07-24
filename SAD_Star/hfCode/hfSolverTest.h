#ifndef HfSolverTEST_H
#define HfSolverTEST_H

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include<iostream>
#include<iomanip>
#include "hfSolver.h"


using namespace std;

//------------------------------------------------------------------------------
struct HfSolverFixture{
  
    HfSolverFixture()
    {
    }
 
    ~HfSolverFixture(){
    }
};



BOOST_FIXTURE_TEST_SUITE( hfSolver, HfSolverFixture )


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( benchmarkTest )
{
    int nMax= 3;
    int nPart= 2;
    double hbarOmega= 10.;
    double hbarc= 197.32891;
    double mc2=938.9059;
    double b= hbarc/sqrt(mc2 * hbarOmega);
    cout<<"b= "<<b<<endl;

    double hfEnergy=0.;

    HfSolver solver;
    solver.setParam(b, nMax, nPart);
    solver.run(hfEnergy);

    BOOST_CHECK_CLOSE(hfEnergy, 25.1393, 1e-3);
}




BOOST_AUTO_TEST_SUITE_END()



#endif
