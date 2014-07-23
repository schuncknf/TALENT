/* 
TALENT school 2014 in Density Functional Theory and Self-Consistent Methods
ECT*, Villazzano (Tn). 
Authors: David Regnier, Simone Orioli.
*/

/*
DETAILS: 
-) Box for the setting of the 1-body matrix elements of a spherically symmetric potential.
-) Requires "Integrator.h" and "HoBasis.h"
-) It defines the value of the potential, calls the HoBasis.h to calculate the integrand and integrates it using Integrator.h
-) To make things simple, the integrand function is defined as a GSL function
*/

/*
 INCLUSIONS
*/
#include <iostream>
#include <armadillo>

#include"integratorGaussLegendreGSL.h"
#include"integratorGaussLegendre.h"
#include "sphericalhofunc.h"
#include<gsl/gsl_integration.h>

#include "VMinnesotaMatrixGenerator.h"

//#include "Integrator.h" //or whatever we will call it
//#include "HoBasis.h" //as before

/*
 FURTHER DECLARATIONS
*/
using namespace arma; //for the armadillo library

/*
 STRUCTURES
*/
typedef struct {
    //nothing inside;
} potStruct;

typedef struct {
    int n1; //main quantum number for the first function
    int n2; //main quantum number for the second function
    int n3; //angular momentum for the first function
    int n4; //angular momentum for the second function
    double r1;
    IntegratorGaussLegendreGSL integrator;
    SphericalHOFunc Rn;
    gsl_function integrandR1R2;
    gsl_function integrandR1;
} structR1;


/*
    FUNCTIONS
*/

//----------------------------------------------------------------------------------------------------------------------------------
static double kR= 1.487; //fm^-2
static double v0R= 200.; //MeV
static double kS= 0.465; //fm^-2
static double v0S= 91.85; //MeV

//----------------------------------------------------------------------------------------------------------------------------------
double integrand(double& r1, double& r2, int& n1, int& n2, int& n3, int& n4, SphericalHOFunc& Rn) {

    double RnFactor= Rn.eval(n1, 0, r1) * Rn.eval(n2, 0, r2) * Rn.eval(n3, 0, r1) * Rn.eval(n4, 0, r2);
    //double RnFactor=1.;
    double gaussR= v0R/2./kR * (exp(-kR*r1*r1 -kR*r2*r2 +2*kR*r1*r2) - exp(-kR*r1*r1 -kR*r2*r2 -2*kR*r1*r2))/2.;
    double gaussS= v0S/2./kS * (exp(-kS*r1*r1 -kS*r2*r2 +2*kS*r1*r2) - exp(-kS*r1*r1 -kS*r2*r2 -2*kS*r1*r2))/2.;

    if(std::isnan((gaussR + gaussS) * RnFactor * r1*r1 *r2*r2) ){
            cout<<"Pb: in integrand: one nan factor:"<<endl;
            cout<<"r0 "<<Rn.eval(n1, 0, r1)<<endl;
            cout<<"r1 "<<Rn.eval(n2, 0, r2)<<endl;
            cout<<"r2 "<<Rn.eval(n3, 0, r1)<<endl;
            cout<<"r3 "<<n4<<"  "<<Rn.eval(n4, 0, r2)<<endl;
            cout<<"gR "<<gaussR<<endl;
            cout<<"gS "<<gaussS<<endl;
            cout<<"Rn "<<RnFactor<<endl;

    }

    return (gaussR + gaussS) * RnFactor * r1*r1 *r2*r2;
}

//----------------------------------------------------------------------------------------------------------------------------------
double integrandR1R2(double r2, void * params){
    structR1 * p= (structR1*) params;
//    cout<<"nval= "<<p->n1<<" "<<p->n2<<" "<<p->n3<<" "<<p->n4<<endl;
    return integrand(p->r1, r2, p->n1, p->n2, p->n3, p->n4, p->Rn);
}


//----------------------------------------------------------------------------------------------------------------------------------
double integrandR1(double r1, void * params){

    structR1 * p= (structR1*) params;
    p->r1= r1;
    p->integrandR1R2.params= p;

    double aMax=500.;
    return p->integrator.integrate(p->integrandR1R2, 0., aMax);
}

//----------------------------------------------------------------------------------------------------------------------------------
double gammaElem(int i, int j, structR1* p, int dim){
    p->n2= i;
    p->n4= j;
    p->integrandR1.params= p;

    double aMax= 500.;
    double elem(0);

    // GammaHF part
    for(int n1=0; n1<dim; n1++){
        p->n1= n1;
        for(int n3=0; n3<dim; n3++){
            p->n3= n3;
            elem+= p->integrator.integrate(p->integrandR1, 0., aMax);
            cout<<n1<<" "<<i<<" "<<n3<<" "<<j<<" "<<p->integrator.integrate(p->integrandR1, 0., aMax)<<endl;
        }
    }
    return elem;
}


//----------------------------------------------------------------------------------------------------------------------------------
void VMinnesotaMatrixGenerator::calcGammaMat(mat& gamma, double b){
    // Creation of the objects to put in the structure
    IntegratorGaussLegendreGSL integrator;
    integrator.setOrder(250);

    SphericalHOFunc Rn;
    Rn.setB(b);

    gsl_function intR1R2Func;
    intR1R2Func.function= integrandR1R2;

    gsl_function intR1Func;
    intR1Func.function= integrandR1;

    structR1 * p= new structR1;
    p->Rn= Rn;
    p->integrator= integrator;
    p->r1= 0.;
    p->integrandR1R2= intR1R2Func;
    p->integrandR1= intR1Func;


    // Loop on matrix elements
    for(int i = 0; i < int(gamma.n_cols); i++ ) {
        cout<<"i "<<i<<endl;
        for (int j = i+1; j < int(gamma.n_cols); j++) {
           gamma(i,j)= gammaElem(i, j, p, gamma.n_cols);
           gamma(j,i)= gamma(i,j);
        }
        gamma(i,i)= gammaElem(i, i, p, gamma.n_cols);
    }

    // Clean: TODO !!!!
//    intR1R2Func.params= NULL;
//    intR1Func.params= NULL;
//    p->integrandR1= gsl_function();
//    p->integrandR1R2= gsl_function();
//    delete(p);

}

//----------------------------------------------------------------------------------------------------------------------------------
double hElem(int& i, int& j, mat& density, mat& gamma, double b){
    double elem(0);

    // GammaHF part
    for(unsigned int n1=0; n1<density.n_cols; n1++){
        for(unsigned int n3=0; n3<density.n_cols; n3++){
            elem+= gamma(n1,n3)* density(n3, n1);
        }
    }

    // HO part
    double mnc2= 938.90590;
    if(i == j){
        elem+= (2*i +3./2.) * mnc2/b/b;
    }

    return elem;
}



//----------------------------------------------------------------------------------------------------------------------------------
void VMinnesotaMatrixGenerator::fillHMatrix(mat& H, mat& density, mat& gamma, double b) {
    if(H.n_cols != density.n_cols || gamma.n_cols != density.n_cols){
        throw invalid_argument( (string(__FILE__)+", "+__FUNCTION__+": H matrix has not the same size than density matrix").c_str());
    }

    // Loop on matrix elements
    for(int i = 0; i < int(density.n_cols); i++ ) {
        cout<<"i "<<i<<endl;
        for (int j = 0; j < int(density.n_cols); j++ ) {
           H(i,j)= hElem(i, j, density, gamma, b);
        }
    }

}
