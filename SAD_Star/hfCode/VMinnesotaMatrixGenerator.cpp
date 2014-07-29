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

#include "integratorGaussLaguerre.h"
#include "integratorGaussLegendre.h"
#include "integratorGaussLegendreGSL.h"

#include "sphericalhofunc.h"
#include "VMinnesotaMatrixGenerator.h"

using namespace arma; //for the armadillo library

/*
 STRUCTURES
*/
typedef struct {
    int n1; //main quantum number for the first function
    int n2; //main quantum number for the second function
    int n3; //angular momentum for the first function
    int n4; //angular momentum for the second function
    double r1; // value of r1
    vector<double> norm;
    SphericalHOFunc Rn;
} Struct1;

typedef struct {
    Struct1* parameters;
    IntegratorGaussLegendre integrator;
    gsl_function integrandR1R2;
} Struct2;


/*
    FUNCTIONS
*/

//----------------------------------------------------------------------------------------------------------------------------------
static double kR= 1.487; //1.487; //fm^-2
static double v0R= 200.; //200.; //MeV
static double kS= 0.465;//0.465; //fm^-2
static double v0S= 91.85;//91.85; //MeV

//----------------------------------------------------------------------------------------------------------------------------------
double integrand(double& r1, double& r2, int n1, int n2, int n3, int n4, SphericalHOFunc& Rn, vector<double>& norm) {

    double rm2= (r1-r2)*(r1-r2);
    double rp2= (r1+r2)*(r1+r2);
    double gaussR= +v0R/4./kR * (exp(-kR*rm2) - exp(-kR*rp2));
    double gaussS= -v0S/4./kS * (exp(-kS*rm2) - exp(-kS*rp2));
    double Rn1Rn2= Rn.eval(n1, 0, r1, norm[n1]) * Rn.eval(n2, 0, r2, norm[n2]);

    double VD=   0.5*(gaussR+ gaussS)* Rn1Rn2* Rn.eval(n3, 0, r1, norm[n3])* Rn.eval(n4, 0, r2, norm[n4])* r1*r2;
    double VEPr= 0.5*(gaussR+ gaussS)* Rn1Rn2* Rn.eval(n4, 0, r1, norm[n4])* Rn.eval(n3, 0, r2, norm[n3])* r1*r2;

    //if(std::isnan(term1) || std::isnan(term2) ){
//            cout<<"Pb: in integrand: one nan factor:"<<endl;
//            cout<<"Rn1 "<<Rn.eval(n1, 0, r1)<<endl;
//            cout<<"Rn2 "<<Rn.eval(n2, 0, r2)<<endl;
//            cout<<"Rn3 "<<Rn.eval(n3, 0, r1)<<endl;
//            cout<<"Rn4 "<<n4<<"  "<<Rn.eval(n4, 0, r2)<<endl;
//            cout<<"gR "<<gaussR<<endl;
//            cout<<"gS "<<gaussS<<endl;
//            cout<<"Rn1 "<<RnFactor1<<endl;
//            cout<<"Rn2 "<<RnFactor2<<endl;
//            cout<<"norm n1 "<<n1<<"  "<<p->norm[n1]<<endl;
//            cout<<"integrand= "<<VD<<" "<<VEPr<<endl;

    //}

    return (VD + VEPr);
}


//----------------------------------------------------------------------------------------------------------------------------------
double integrandR1R2(double r2, void * params){
    Struct1 * p= (Struct1*) params;
//    cout<<"nval= "<<p->n1<<" "<<p->n2<<" "<<p->n3<<" "<<p->n4<<endl;
    return integrand(p->r1, r2, p->n1, p->n2, p->n3, p->n4, p->Rn, p->norm);
}


//----------------------------------------------------------------------------------------------------------------------------------
double integrandR1(double r1, void * params){

    Struct2 * p= (Struct2*) params;
    p->parameters->r1= r1;

    return p->integrator.integrate0ToInf(p->integrandR1R2);
}


//----------------------------------------------------------------------------------------------------------------------------------
void VMinnesotaMatrixGenerator::calc2BodyMat(TwoBodyMat& Vabcd, double& b, int order){
    int dim= Vabcd.size();

    // Struct1
    SphericalHOFunc Rn;
    Rn.setB(b);

    vector<double> norm;
    for(int i=0; i<dim; i++){
        norm.push_back(Rn.norm(i,0));
    }

    Struct1* param1= new Struct1;
    param1->norm= norm;
    param1->Rn= Rn;

    // Struct2
    gsl_function intR1R2Func;
    intR1R2Func.function= integrandR1R2;
    intR1R2Func.params= param1;

    IntegratorGaussLegendre integrator;
    integrator.setTableDir("../gen_legendre");
    //integrator.setTableDir("../gen_laguerre");
    integrator.setOrder(order);

    Struct2 * param2= new Struct2;
    param2->integrator= integrator;
    param2->integrandR1R2= intR1R2Func;
    param2->parameters= param1;

    // Loop on matrix elements
    gsl_function intR1Func;
    intR1Func.function= integrandR1;
    intR1Func.params= param2;

    TwoBodyMat isFilled (dim, vector<vector<vector<double> > >(dim, vector<vector<double> >(dim, vector<double>(dim, 0.))));

    for(int n1 = 0; n1 < dim; n1++ ) {
        param1->n1= n1;
        for (int n2 = n1; n2 < dim; n2++) {
           param1->n2=n2;
            for(int n3=0; n3<dim; n3++){ //n1
                param1->n3=n3;
                for(int n4=n3; n4<dim; n4++){ //max(n3,n2)
                    param1->n4=n4;

                    // Calc element:
                    if( isFilled[n1][n2][n3][n4] < 0.5){
                        double elem= integrator.integrate0ToInf(intR1Func);

                        Vabcd[n1][n2][n3][n4]= elem;
                        Vabcd[n1][n2][n4][n3]= elem;

                        Vabcd[n2][n1][n3][n4]= elem;
                        Vabcd[n2][n1][n4][n3]= elem;

                        Vabcd[n3][n4][n1][n2]= elem;
                        Vabcd[n3][n4][n2][n1]= elem;

                        Vabcd[n4][n3][n1][n2]= elem;
                        Vabcd[n4][n3][n2][n1]= elem;

                        isFilled[n1][n2][n3][n4]= 1.;
                        isFilled[n1][n2][n4][n3]= 1.;

                        isFilled[n2][n1][n3][n4]= 1.;
                        isFilled[n2][n1][n4][n3]= 1.;

                        isFilled[n3][n4][n1][n2]= 1.;
                        isFilled[n3][n4][n2][n1]= 1.;

                        isFilled[n4][n3][n1][n2]= 1.;
                        isFilled[n4][n3][n2][n1]= 1.;

//                        cout<<n1<<" "<<n2<<" "<<n3<<" "<<n4<<" "<<setprecision(25)<<elem<<endl;
                    } // if is not filled

                } //end n4
            } // end n3
        } // end n2
    } // end n1

    // Clean:
    //delete param1;
    //delete param2;
}


//----------------------------------------------------------------------------------------------------------------------------------
double hElem(int& i, int& j, mat& density, vector<vector<vector<vector<double> > > >& Vabcd, double& b){
    double elem(0.);

    int sigma2= i%2;
    int sigma4= j%2;

    // Density dependent part
    for(unsigned int n1=0; n1<density.n_cols/2; n1++){
        for(unsigned int n3=0; n3<density.n_cols/2; n3++){
            if(sigma2 == sigma4){
                elem+= Vabcd[n1][i/2][n3][j/2] * density(2*n3, 2*n1);
            }
        }
    }

    // HO part
    double hbarOmega= HBARC*HBARC/MNC2/b/b;
    if(i == j){
        elem+= (2.*(i/2) +1.5) * hbarOmega;
    }

    return elem;
}



//----------------------------------------------------------------------------------------------------------------------------------
void VMinnesotaMatrixGenerator::fillHMatrix(mat& H, mat& density, TwoBodyMat& Vabcd, double& b) {
    if(H.n_cols != density.n_cols || Vabcd.size() != density.n_cols/2){
        throw invalid_argument( (string(__FILE__)+", "+__FUNCTION__+": H matrix has not the same size than density matrix").c_str());
    }

    // Loop on matrix elements
    for(int i = 0; i < int(H.n_cols); i++ ) {
        for (int j = i+1; j < int(H.n_cols); j++ ) {
           H(i,j)= hElem(i, j, density, Vabcd, b);
           H(j,i)= H(i,j);
        }
        H(i,i)= hElem(i, i, density, Vabcd, b);
    }
}
