#ifndef VMinnesotaMatrixGeneratorTEST_H
#define VMinnesotaMatrixGeneratorTEST_H

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include<iostream>
#include<iomanip>

#include<armadillo>

#include "constants.h"
#include "VMinnesotaMatrixGenerator.h"


using namespace std;
using namespace VMinnesotaMatrixGenerator;

//------------------------------------------------------------------------------
struct VMinnesotaMatrixGeneratorFixture{

    VMinnesotaMatrixGeneratorFixture(){
    }
 
    ~VMinnesotaMatrixGeneratorFixture(){
    }
};


BOOST_FIXTURE_TEST_SUITE( vMinnesotaMatrixGenerator, VMinnesotaMatrixGeneratorFixture )


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( benchmarkTest )
{
    int nMax=4;
    double hbarOmega= 10.;
    double b= HBARC/sqrt(MNC2*hbarOmega);

    // Vabsd
    TwoBodyMat Vabcd(nMax+1, vector<vector<vector<double> > >(nMax+1, vector<vector<double> >(nMax+1, vector<double>(nMax+1, 0.))));
    calc2BodyMat(Vabcd, b);

    // Benchmark
    double error=0.;

    ifstream input("./2BodyMatBenchmark.dat");
    int i1, i2, i3, i4;
    double elemBench;

    if(!input){
        throw logic_error("File 2BodyMatBenchmark.dat does not exist");
    }
    while(!input.eof()){
        input>>i1>>i2>>i3>>i4>>elemBench;
        double eVal= abs(Vabcd[i1][i2][i3][i4] - elemBench);
        cout<<Vabcd[i1][i2][i3][i4]<<"  "<<elemBench<<endl;
        error+= eVal;
        if(eVal > 1e-5){
//            cout<<i1/2<<" "<<i2/2<<" "<<i3/2<<" "<<i4/2<<" e="<<eVal<<"  eB= "<<elemBench<<""<<endl;
        }
        input>>elemBench;
    }


    BOOST_CHECK_CLOSE(error+1., 1., 1e-6);
}





BOOST_AUTO_TEST_SUITE_END()



#endif
