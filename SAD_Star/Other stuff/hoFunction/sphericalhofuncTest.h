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
       return func.eval(n1, 0, r) * func.eval(n2, 0, r) * r * r;
}

//------------------------------------------------------------------------------
//BOOST_AUTO_TEST_CASE( orthonormalTest )
//{
//    IntegratorGaussLegendre integrator;
//    integrator.readTables("lgvalues-weights.php", "lgvalues-abscissa.php");
//    integrator.setOrder(60);

//    SphericalHOFunc funcR;
//    funcR.setB(0.5);


//    auto integrand= bind(hoProdTest, placeholders::_1, funcR, 19, 19);
//    cout<<hoProdTest(0.5, funcR, 19, 19)<<"  "<<integrand(0.5)<<endl;
//    cout<<hoProdTest(2., funcR, 19, 19)<<"  "<<integrand(2.)<<endl;


//    double cumulError= 0;
//    int nMax= 50;
//    for(int n1=nMax-nMax; n1<nMax; n1++){
//        for(int n2=0; n2<nMax; n2++){
//            auto integrand= bind(hoProdTest, placeholders::_1, funcR, n1, n2);
//            double res= integrator.integrate0ToInf(integrand);
//            double error= abs(res - int(n1==n2));
//            if(n1==n2){
//                cout<<n1<<"  " <<n2<<"  "<<res<<endl;
//            }
//            cumulError+= error;
//        }
//    }

//    BOOST_CHECK_CLOSE(cumulError +1., 1., 1e-6);
//}


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( orthonormal2Test )
{
    IntegratorGaussLaguerre integrator;
    integrator.readTables( "../gen_laguerre");
    integrator.setOrder(1000);

    SphericalHOFunc funcR;
    funcR.setB(0.5);


    double cumulError= 0;

    int nMax= 20;
    for(int n1=0; n1<nMax; n1++){
        for(int n2= 0; n2< nMax; n2++){
            auto integrand= bind(hoProdTest, placeholders::_1, funcR, n1, n2);
            double res= integrator.integrate0ToInf(integrand);
            double error= abs(res - double(n1==n2));
//            if(n1==n2){
//                cout<<n1<<"  " <<n2<<"  "<<res<<endl;
//            }
            cumulError+= error;
        }
    }

    BOOST_CHECK_CLOSE(cumulError +1., 1., 0.1);
}


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( plotTest ){
    SphericalHOFunc funcR;
    funcR.setB(0.5);

    ofstream output("r19.dat");
    for(double r=0; r<40; r+=0.01){
        output<<r <<"   "<<setprecision(10)<<funcR.eval(19,0, r)<<endl;
    }
    output.close();
}


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( r19Test ){
    SphericalHOFunc funcR;
    funcR.setB(0.5);

    double tmp;
    vector<double> val;
    ifstream input("r19Benchmark.dat");
    while(!input.eof()){
        input>>tmp>>tmp;
        val.push_back(tmp);
    }

    double cumulError(0);
    double r=0;
    for(unsigned int i=0; i<val.size(); i++){
        cumulError+= abs(val[i] - funcR.eval(19, 0, r));
        r+= 0.01;
    }
    BOOST_CHECK_CLOSE(cumulError+1., 1., 1e-5);
}


BOOST_AUTO_TEST_SUITE_END()

#endif // SPHERICALHOFUNCTEST_H
