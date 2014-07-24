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
static double kR= 1.487; //1.487; //fm^-2
static double v0R= 200.; //200.; //MeV
static double kS= 0.465;//0.465; //fm^-2
static double v0S= 91.85;//91.85; //MeV

//----------------------------------------------------------------------------------------------------------------------------------
double integrand(double& r1, double& r2, int& n1, int& n2, int& n3, int& n4, SphericalHOFunc& Rn) {

    double gaussR= +v0R/4./kR * (exp(-kR*r1*r1 -kR*r2*r2 +2*kR*r1*r2) - exp(-kR*r1*r1 -kR*r2*r2 -2*kR*r1*r2));
    double gaussS= -v0S/4./kS * (exp(-kS*r1*r1 -kS*r2*r2 +2*kS*r1*r2) - exp(-kS*r1*r1 -kS*r2*r2 -2*kS*r1*r2));
    double Rn1Rn2= Rn.eval(n1, 0, r1) * Rn.eval(n2, 0, r2);

    double RnFactor1= Rn1Rn2 * Rn.eval(n3, 0, r1) * Rn.eval(n4, 0, r2);
    double term1= 0.5*(gaussR + gaussS) * RnFactor1 * r1*r2;

    double RnFactor2= Rn1Rn2 * Rn.eval(n4, 0, r1) * Rn.eval(n3, 0, r2);
    double term2= 0.5*(gaussR + gaussS) * RnFactor2 * r1*r2;

    if(std::isnan(term1) || std::isnan(term2) ){
            cout<<"Pb: in integrand: one nan factor:"<<endl;
            cout<<"r0 "<<Rn.eval(n1, 0, r1)<<endl;
            cout<<"r1 "<<Rn.eval(n2, 0, r2)<<endl;
            cout<<"r2 "<<Rn.eval(n3, 0, r1)<<endl;
            cout<<"r3 "<<n4<<"  "<<Rn.eval(n4, 0, r2)<<endl;
            cout<<"gR "<<gaussR<<endl;
            cout<<"gS "<<gaussS<<endl;
            cout<<"Rn1 "<<RnFactor1<<endl;
            cout<<"Rn2 "<<RnFactor2<<endl;

    }

    return term1 + term2;
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

    double aMax=200. * p->Rn.getB();
    return p->integrator.integrate(p->integrandR1R2, 0., aMax);
}

////----------------------------------------------------------------------------------------------------------------------------------
//double VabcdElem(structR1* p, int dim){
//    p->n2= i;
//    p->n4= j;
//    p->integrandR1.params= p;

//    double aMax= 500.;
//    double elem(0);

//    // GammaHF part
//    for(int n1=0; n1<dim; n1++){
//        p->n1= n1;
//        for(int n3=0; n3<dim; n3++){
//            p->n3= n3;
//            elem+= p->integrator.integrate(p->integrandR1, 0., aMax);
//            cout<<n1<<" "<<i<<" "<<n3<<" "<<j<<" "<<p->integrator.integrate(p->integrandR1, 0., aMax)<<endl;
//        }
//    }
//    return elem;
//}


//----------------------------------------------------------------------------------------------------------------------------------
void VMinnesotaMatrixGenerator::calc2BodyMat(vector<vector<vector<vector<double> > > >& Vabcd, double b){
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

    double aMax=200. * p->Rn.getB();


    // Loop on matrix elements
    int dim= Vabcd.size();
    for(int n1 = 0; n1 < dim; n1++ ) {
        p->n1= n1;
        cout<<"n1 "<<n1<<endl;
        for (int n2 = 0; n2 < dim; n2++) {
            p->n2=n2;
            for(int n3=0; n3<dim; n3++){ //n1
                p->n3=n3;
                for(int n4=n3; n4<dim; n4++){ //max(n3,n2)
                    p->n4=n4;

                    // Calc element:
                    p->integrandR1.params= p;
                    double elem= p->integrator.integrate(p->integrandR1, 0., aMax);
                    Vabcd[n1][n2][n3][n4]= elem;
                    Vabcd[n1][n2][n4][n3]= elem;

//                    Vabcd[n1][n4][n3][n2]= elem;
//                    Vabcd[n1][n4][n2][n3]= elem;

//                    Vabcd[n3][n2][n1][n4]= elem;
//                    Vabcd[n3][n2][n4][n1]= elem;
                    cout<<n1<<" "<<n2<<" "<<n3<<" "<<n4<<" "<<setprecision(12)<<elem<<endl;

                } //end n4
            } // end n3
        } // end n2
    } // end n1



    // Clean: TODO !!!!
//    intR1R2Func.params= NULL;
//    intR1Func.params= NULL;
//    p->integrandR1= gsl_function();
//    p->integrandR1R2= gsl_function();
//    delete(p);

}

//----------------------------------------------------------------------------------------------------------------------------------
double hElem(int& i, int& j, mat& density, vector<vector<vector<vector<double> > > >& Vabcd, double b){
    double elem(0);

    int sigmaI= i%2;
    int sigmaJ= j%2;

    // density dependent part
    for(unsigned int n1=0; n1<density.n_cols/2; n1++){
        for(unsigned int n3=0; n3<density.n_cols/2; n3++){
            //cout<<(sigmai == sigmaj)<<endl;
            if(sigmaI == sigmaJ){
                elem+= 0.5*Vabcd[n1][i/2][n3][j/2]* density(2*n3, 2*n1);
            }

//            if(sigma2 == sigma4){
//                elem+= Vabcd[k][i/2][l][j/2]* density(2*l + (1-sigma2), 2*k + (1-sigma2));
//            }
//            else{
//                elem+= Vabcd[k][i/2][l][j/2]* density(2*l + sigma2, 2*k + sigma4);
//            }
        }
    }

    // HO part
    double mnc2= 938.9059;
    double hbarc= 197.32891;
    double hbarOmega= hbarc*hbarc/mnc2/b/b;
    if(i == j){
        elem+= (2*(i/2) +3./2.) * hbarOmega;
    }

    return elem;
}



//----------------------------------------------------------------------------------------------------------------------------------
void VMinnesotaMatrixGenerator::fillHMatrix(mat& H, mat& density, vector<vector<vector<vector<double> > > >& Vabcd, double b) {
    if(H.n_cols != density.n_cols || Vabcd.size() != density.n_cols/2){
        throw invalid_argument( (string(__FILE__)+", "+__FUNCTION__+": H matrix has not the same size than density matrix").c_str());
    }

    // Loop on matrix elements
    for(int i = 0; i < int(H.n_cols); i++ ) {
        for (int j = 0; j < int(H.n_cols); j++ ) {
           H(i,j)= hElem(i, j, density, Vabcd, b);
        }
    }

}
