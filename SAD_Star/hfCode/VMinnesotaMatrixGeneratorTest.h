#ifndef VMinnesotaMatrixGeneratorTEST_H
#define VMinnesotaMatrixGeneratorTEST_H

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include<iostream>
#include<iomanip>

#include<armadillo>

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
    int nMax=2;
    int nPart=2;
    double hbarOmega= 10.;
    double hbarc= 197.32891; //MeV.fm
    double mc2=938.9059;
    double c=1e23; // fm/s
    double b=sqrt(mc2/hbarOmega);

    // density
    mat density= zeros(nMax, nMax);
    for(int i=0; i<nPart; i++){
        density(i,i)=1.;
    }
    cout<<"density size= "<<density.n_cols<<endl;

    // Vabsd
    TwoBodyMat Vabcd(nMax, vector<vector<vector<double> > >(nMax, vector<vector<double> >(nMax, vector<double>(nMax, 0.))));
    calc2BodyMat(Vabcd, b);
    cout<<"Vabcd"<<endl;

    // Hamiltonian
    mat h= zeros(nMax,nMax);
    fillHMatrix(h, density, Vabcd, b);

    cout<<"hMatrix:"<<endl;
    cout<<h<<endl;
    BOOST_CHECK(true);
}





BOOST_AUTO_TEST_SUITE_END()



#endif
