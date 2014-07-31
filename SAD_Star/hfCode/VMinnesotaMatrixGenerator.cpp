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
#include <fstream>
#include <armadillo>

#include "integratorGaussLaguerre.h"
#include "integratorGaussLegendre.h"
#include "integratorGaussLegendreGSL.h"

#include "sphericalhofunc.h"
#include "VMinnesotaMatrixGenerator.h"

using namespace arma; //for the armadillo library
using namespace VMinnesotaMatrixGenerator;

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
TwoBodyMat VMinnesotaMatrixGenerator::emptyMat(int nMax, int lMax){
    FourIndiceMat fourIMat( (nMax+1), vector<vector<vector<double> > >((nMax+1), vector<vector<double> >((nMax+1), vector<double>((nMax+1), 0.))));
    TwoBodyMat res(2*lMax+2,  vector<FourIndiceMat>(2*lMax+2, fourIMat));
    return res;
}

////----------------------------------------------------------------------------------------------------------------------------------
//double integrand(double& r1, double& r2, int n1, int n2, int n3, int n4, SphericalHOFunc& Rn, vector<double>& norm) {

//    double rm2= (r1-r2)*(r1-r2);
//    double rp2= (r1+r2)*(r1+r2);
//    double gaussR= +v0R/4./kR * (exp(-kR*rm2) - exp(-kR*rp2));
//    double gaussS= -v0S/4./kS * (exp(-kS*rm2) - exp(-kS*rp2));
//    double Rn1Rn2= Rn.eval(n1, 0, r1, norm[n1]) * Rn.eval(n2, 0, r2, norm[n2]);

//    double VD=   0.5*(gaussR+ gaussS)* Rn1Rn2* Rn.eval(n3, 0, r1, norm[n3])* Rn.eval(n4, 0, r2, norm[n4])* r1*r2;
//    double VEPr= 0.5*(gaussR+ gaussS)* Rn1Rn2* Rn.eval(n4, 0, r1, norm[n4])* Rn.eval(n3, 0, r2, norm[n3])* r1*r2;

//    //if(std::isnan(term1) || std::isnan(term2) ){
////            cout<<"Pb: in integrand: one nan factor:"<<endl;
////            cout<<"Rn1 "<<Rn.eval(n1, 0, r1)<<endl;
////            cout<<"Rn2 "<<Rn.eval(n2, 0, r2)<<endl;
////            cout<<"Rn3 "<<Rn.eval(n3, 0, r1)<<endl;
////            cout<<"Rn4 "<<n4<<"  "<<Rn.eval(n4, 0, r2)<<endl;
////            cout<<"gR "<<gaussR<<endl;
////            cout<<"gS "<<gaussS<<endl;
////            cout<<"Rn1 "<<RnFactor1<<endl;
////            cout<<"Rn2 "<<RnFactor2<<endl;
////            cout<<"norm n1 "<<n1<<"  "<<p->norm[n1]<<endl;
////            cout<<"integrand= "<<VD<<" "<<VEPr<<endl;

//    //}

//    return (VD + VEPr);
//}


////----------------------------------------------------------------------------------------------------------------------------------
//double integrandR1R2(double r2, void * params){
//    Struct1 * p= (Struct1*) params;
////    cout<<"nval= "<<p->n1<<" "<<p->n2<<" "<<p->n3<<" "<<p->n4<<endl;
//    return integrand(p->r1, r2, p->n1, p->n2, p->n3, p->n4, p->Rn, p->norm);
//}


////----------------------------------------------------------------------------------------------------------------------------------
//double integrandR1(double r1, void * params){

//    Struct2 * p= (Struct2*) params;
//    p->parameters->r1= r1;

//    return p->integrator.integrate0ToInf(p->integrandR1R2);
//}


//----------------------------------------------------------------------------------------------------------------------------------
void VMinnesotaMatrixGenerator::calc2BodyMat(TwoBodyMat& V, double& bParam, int order){
    // Launch the Morten's program


    // Read legend:
    cout<<"read legend"<<endl;
    string fileName="../relcom2labsystem/spM.dat";
    ifstream legInput(fileName.c_str());
    if(!legInput){
        throw logic_error("in VMinnesotaMatrixGenerator::calc2BodyMat, file missing");
    }

    vector<int> n(1,-1), l(1,-1), s(1,-1), twomJ(1,-1), twoTau(1,-1);
    string word;
    int nVal, lVal,twoJ, twomJVal,twoTauVal;

    getline(legInput, word);
    while(!legInput.eof()){
        legInput>>word>>word>>word;
        legInput>>nVal>>lVal>>twoJ>>twomJVal>>twoTauVal;

        n.push_back(nVal);
        l.push_back(lVal);
        s.push_back(twoJ-2*lVal);
        twomJ.push_back(twomJVal);
        twoTau.push_back(twoTauVal);

        //cout<<nVal<<" "<<twoTauVal<<endl;
    }


    cout<<"read values"<<endl;
    fileName="../relcom2labsystem/VM-scheme.dat";
    ifstream valInput(fileName.c_str());
    if(!valInput){
        throw logic_error("in VMinnesotaMatrixGenerator::calc2BodyMat, file missing");
    }

    getline(valInput, word);
    getline(valInput, word);
    cout<<word<<endl;
    int a,b,c,d;
    double val;
    int count=0;
    while(!valInput.eof()){
        valInput>>a>>b>>c>>d>>val;

        int i1= 2* l[a] + (1+s[a])/2;
        int i2= 2* l[b] + (1+s[b])/2;
        double twoJ= s[a]+2*l[a];
        double twoJPrime= s[b]+2*l[b];

        V[i1][i2][n[a]][n[b]][n[c]][n[d]] += val/(twoJ+1.)/(twoJPrime+1.);

        if(count%1 == 0){
            cout<<count<<endl;
            cout<<"abcd " << a<<" "<<b<<" "<<c<<" "<<d<<" "<<val<<endl;
        }
        count++;
    }
}


//----------------------------------------------------------------------------------------------------------------------------------
double hElem(int& k, int& i, int& j, vector<mat>& density, TwoBodyMat& V , double& b){
    double elem(0.);

   // Density dependent part
    for(unsigned int n2=0; n2<density[k].n_cols; n2++){
        for(unsigned int n4=0; n4<density[k].n_cols; n4++){
            for(int kPrime=0; kPrime<density.size(); kPrime++){
                elem+= V[k][kPrime][i][n2][j][n4] * density[kPrime](n4, n2);
            }
        }
    }

    // HO part
    double hbarOmega= HBARC*HBARC/MNC2/b/b;
    if(i == j){
        elem+= (2*i+ k/2 + 1.5) * hbarOmega;
    }

    return elem;
}



//----------------------------------------------------------------------------------------------------------------------------------
void VMinnesotaMatrixGenerator::fillHMatrix(vector<mat>& H, mat& density, TwoBodyMat& Vabcd, double& b) {
    if(H.n_cols != density.n_cols || Vabcd.size() != density.n_cols/2){
        throw invalid_argument( (string(__FILE__)+", "+__FUNCTION__+": H matrix has not the same size than density matrix").c_str());
    }

    // Loop on matrix elements
    for(int k= 0; k< 2*lMax+1; k++){
        for(int i = 0; i < int(H.n_cols); i++ ) {
            for (int j = i+1; j < int(H.n_cols); j++ ) {
               H[k](i,j)= hElem(k, i, j, density, Vabcd, b);
               H[k](j,i)= H[k](i,j);
            }
            H[k](i,i)= hElem(i, i, density, Vabcd, b);
        }
    }

}
