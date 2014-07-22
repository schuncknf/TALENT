#include <math.h>
#include<iostream>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_integration.h>

#include <armadillo>

#include "VCoulombMatrixGenerator.h"

using namespace std;
using namespace arma;




// Hydrogen
int main() {

    int dim=20;
    mat A(dim, dim, fill::zeros);
    setMatrix(A, dim, 0.5);


    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, A);

    cout<<eigval[0]<<endl;



    return 0;
}
