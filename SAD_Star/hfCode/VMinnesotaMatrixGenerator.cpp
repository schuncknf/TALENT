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

//------------------------------------------------------------------------------
FourIndiceMat VMinnesotaMatrixGenerator::emptyFourIndiceMat(int dim){
    return FourIndiceMat( (dim), vector<vector<vector<double> > >((dim), vector<vector<double> >((dim), vector<double>((dim), 0.))));
}


//----------------------------------------------------------------------------------------------------------------------------------
TwoBodyMat VMinnesotaMatrixGenerator::emptyTwoBodyMat(int NMax, int lMax){
    FourIndiceMat fourIMat= emptyFourIndiceMat(NMax/2 +1); // No need of all this space because nMax depends on l ...
    TwoBodyMat res(2*lMax+2,  vector<FourIndiceMat>(2*lMax+2, fourIMat));
    return res;
}


//------------------------------------------------------------------------------
void readMortenMatrix(FourIndiceMat& fourIMat, string fileName){
    ifstream valInput(fileName.c_str());
    if(!valInput){
        throw logic_error("in VMinnesotaMatrixGenerator::readMortenMatrix, file missing");
    }

    // Skip 2 first lines
    string word;
    getline(valInput, word);
    getline(valInput, word);

    int a,b,c,d;
    double val;

    while(!valInput.eof()){
        valInput>>a>>b>>c>>d>>val;

        fourIMat[a][b][c][d]= val;
//        fourIMat[a][b][d][c]= -val;

//        fourIMat[b][a][c][d]= -val;
//        fourIMat[b][a][d][c]= val;

//        fourIMat[c][d][a][b]= val;
//        fourIMat[c][d][b][a]= -val;

//        fourIMat[d][c][a][b]= -val;
//        fourIMat[d][c][b][a]= val;


//        int i1= 2* l[a] + (1+s[a])/2;
//        int i2= 2* l[b] + (1+s[b])/2;
//        double twoJ= s[a]+2*l[a];
//        double twoJPrime= s[b]+2*l[b];

    }

//    // Antisymmetry test
//    for(int a=0; a<21; a++){
//        for(int b=0; b<21; b++){
//            for(int c=0; c<21; c++){
//                for(int d=0; d<21; d++){

//                    if(fourIMat[a][b][c][d] != - 200){
//                        cout<<fourIMat[a][b][c][d]<<endl;


//                    }

//                }
//            }
//        }
//    }
}


//------------------------------------------------------------------------------
/**
 *  The int vectors must be empty
 */
void readMortenIndiceTable(vector<int>& n, vector<int>& l, vector<int>& twoJ, vector<int>& twoMj, vector<int>& twoTau, string fileName){
    ifstream legInput(fileName.c_str());
    if(!legInput){
        throw logic_error("in VMinnesotaMatrixGenerator::readMortenIndiceTable, file missing");
    }

    // Put dummy values in the 0 indice
    n.push_back(-1);
    l.push_back(0);
    twoJ.push_back(0);
    twoMj.push_back(0);
    twoTau.push_back(0);

    // Read values from file
    int nVal, lVal,twoJVal, twoMjVal,twoTauVal;
    string word;

    getline(legInput, word);
    while(!legInput.eof()){
        legInput>>word>>word>>word;
        legInput>>nVal>>lVal>>twoJVal>>twoMjVal>>twoTauVal;
//        cout<<nVal<<" "<<lVal<<endl;

        n.push_back(nVal);
        l.push_back(lVal);
        twoJ.push_back(twoJVal);
        twoMj.push_back(twoMjVal);
        twoTau.push_back(twoTauVal);
    }
}

//------------------------------------------------------------------------------
int getMortenInd(vector<int>& n, vector<int>& l, vector<int>& twoJ, vector<int>& twoMj, vector<int>& twoTau, int nVal, int lVal, int twoJVal, int twoMjVal, int twoTauVal){
    int i=0;
    while(i<int(n.size()) && (n[i]!=nVal || l[i]!= lVal || twoJ[i]!= twoJVal || twoMj[i]!= twoMjVal || twoTau[i]!=twoTauVal) ){
        i++;
    }
    if(i == int(n.size())){
        throw logic_error("In getMortenInd, quantum number not in the table");
    }
    return i;
}

//----------------------------------------------------------------------------------------------------------------------------------
void VMinnesotaMatrixGenerator::calc2BodyMat(TwoBodyMat& V, double& bParam, int order, int NMax){
    // Launch the Morten's program

    // Read Morten's matrix elements
    cout<<"read table"<<endl;
    vector<int> n,l,twoJ,twoMj, twoTau;
    readMortenIndiceTable(n, l, twoJ, twoMj, twoTau, "../relcom2labsystem/spM.dat");

    cout<<"read matrix"<<endl;
    FourIndiceMat fourIMat= emptyFourIndiceMat(n.size());
    readMortenMatrix(fourIMat, "../relcom2labsystem/VM-scheme.dat");

    // Fill the V matrix
    int kMax=2*NMax;
    int kMin=0;
    cout<<"kMax= "<<kMax<<endl;

    double sum;
    int i1, i2, i3, i4, twoJ1, twoJ2, l1, l2;
    for(int k1=kMin; k1<=kMax; k1++){
        l1= (k1+k1%2)/2;
        twoJ1= k1-k1%2+1;

        for(int k2=kMin; k2<=kMax; k2++){
            l2= (k2+k2%2)/2;
            twoJ2= k2-k2%2+1;

            for(int n1=0; n1<= (NMax-l1)/2; n1++){
                for(int n2=0; n2<= (NMax-l2)/2; n2++){
                    for(int n3=0; n3<= (NMax-l1)/2; n3++){
                        for(int n4=0; n4<= (NMax-l2)/2; n4++){

                            sum=0;
                            for(int twoMj1= -twoJ1; twoMj1<= twoJ1; twoMj1+=2){
                                for(int twoMj2= -twoJ2; twoMj2<= twoJ2; twoMj2+=2){
                                    i1= getMortenInd(n, l, twoJ, twoMj, twoTau, n1, l1, twoJ1, twoMj1, 1);
                                    i2= getMortenInd(n, l, twoJ, twoMj, twoTau, n2, l2, twoJ2, twoMj2, 1);
                                    i3= getMortenInd(n, l, twoJ, twoMj, twoTau, n3, l1, twoJ1, twoMj1, 1);
                                    i4= getMortenInd(n, l, twoJ, twoMj, twoTau, n4, l2, twoJ2, twoMj2, 1);
                                    sum+= fourIMat[i1][i2][i3][i4];
                                }
                            }// End sum over m1 m2

                            V[k1][k2][n1][n2][n3][n4]= sum/double(twoJ1+1)/double(twoJ2+1);
                        } // end n4
                    }
                }
            }// end n1
        } // end k2
    } // end k1

}


//----------------------------------------------------------------------------------------------------------------------------------
double hElem(int& k, int& i, int& j, vector<mat>& density, TwoBodyMat& V , double& b, int NMax){
    double elem(0.);

    int twoJ= k-k%2+1;
   // Density dependent part
    for(unsigned int kPrime=0; kPrime<density.size(); kPrime++){
        int lPrime= (kPrime+kPrime%2)/2;
        int twoJPrime= kPrime-kPrime%2+1;

        for(unsigned int n2=0; n2<= (NMax-lPrime)/2; n2++){
            for(unsigned int n4=0; n4<= (NMax-lPrime)/2; n4++){

                elem+= V[k][kPrime][i][n2][j][n4] * density[kPrime](n4, n2) *(twoJPrime+1);
//                cout<<V[k][kPrime][i][n2][j][n4]<<endl;
            }
        }
    }

    // HO part
    double hbarOmega= HBARC*HBARC/MNC2/b/b;
    if(i == j){
        int l= (k+k%2)/2;
        elem+= (2*i+ l + 1.5) * hbarOmega;
    }

    return elem;
}



//----------------------------------------------------------------------------------------------------------------------------------
void VMinnesotaMatrixGenerator::fillHMatrix(vector<mat>& H, vector<mat>& density, TwoBodyMat& Vabcd, double& b, int lMax) {
//    if(H[0].n_cols != density[0].n_cols || Vabcd.size() != density[0].n_cols){
//        throw invalid_argument( (string(__FILE__)+", "+__FUNCTION__+": H matrix has not the same size than density matrix").c_str());
//    }

    int NMax= lMax;

    // Loop on matrix elements
    for(int k= 0; k<= 2*lMax; k++){
        int l= (k+k%2)/2;
        for(int i = 0; i <= (NMax-l)/2; i++ ) {
            for (int j = i+1; j <= (NMax-l)/2; j++ ) {
               H[k](i,j)= hElem(k, i, j, density, Vabcd, b, NMax);
               H[k](j,i)= H[k](i,j);
            }
            H[k](i,i)= hElem(k, i, i, density, Vabcd, b, NMax);
        }
    }

}
