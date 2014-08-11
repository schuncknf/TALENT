#ifndef HfSolverTEST_H
#define HfSolverTEST_H

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include<iostream>
#include<iomanip>
#include "constants.h"

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
    int nMax= 4;
    int nPart= 2;
    double hbarOmega= 10.;
    double b= HBARC/sqrt(MNC2 * hbarOmega);

    double hfEnergy=0.;

    HfSolver solver;
    solver.setParam(b, nMax, nPart);
    solver.run(hfEnergy);

    BOOST_CHECK_CLOSE(hfEnergy, 25.139312, 4e-6);
}




BOOST_AUTO_TEST_SUITE_END()



#endif
