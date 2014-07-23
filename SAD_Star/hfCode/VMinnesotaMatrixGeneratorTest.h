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
BOOST_AUTO_TEST_CASE( notEmptyMatrixTest )
{
    int dim=2;
    int nPart=2;
    double b=9.6897156;

    // density
    mat density= zeros(dim, dim);
    for(int i=0; i<nPart; i++){
        density(i,i)=1.;
    }
    cout<<"density size= "<<density.n_cols<<endl;

    // Gamma
    mat gamma= zeros(dim,dim);
    calcGammaMat(gamma, b);
    cout<<"gamma"<<endl;
    cout<<gamma<<endl;

    // Hamiltonian
    mat h= zeros(dim,dim);
    fillHMatrix(h, density, gamma, b);

    cout<<"hMatrix:"<<endl;
    cout<<h<<endl;
    BOOST_CHECK(true);
}





BOOST_AUTO_TEST_SUITE_END()



#endif
