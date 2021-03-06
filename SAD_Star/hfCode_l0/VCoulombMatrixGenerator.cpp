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

#include "VCoulombMatrixGenerator.h"

//#include "Integrator.h" //or whatever we will call it
//#include "HoBasis.h" //as before

/*
 FURTHER DECLARATIONS
*/
using namespace arma; //for the armadillo library
using namespace VCoulombMatrixGenerator;


/*
    FUNCTIONS
*/
//----------------------------------------------------------------------------------------------------------------------------------
double VCoulombMatrixGenerator::potential(double r, void *params) {
    //Coulomb potential in natural units

   // potStruct *pS = (potStruct *)params;
    
    return -(1./r);
}

//----------------------------------------------------------------------------------------------------------------------------------
double VCoulombMatrixGenerator::Vij(double r, void *params) {
    //generates the integrand R_nl(r)*V(r)R_n'l'(r) using the function potential
    
    matElStruct *p = (matElStruct *)params; //casting of the params pointer to a prms structure
    
    int i = p->i; 
    int j = p->j;
    int l = p->l;
    double b = p->b;
    double y, z;
    
    SphericalHOFunc func;
    func.setB(b);

    //the name of the function generate_basis has to be substitute with the one we will choose
    y = func.eval(i, l, r); //obtain the value of the basis function for the different values of i and j
    z = func.eval(j, l, r);
    
    return y * z * r * r * potential(r, p->pS); //r^2 because of the jacobian of spherical coordinates
    //explicitely ignoring spherical harmonics in case of l=0, has to be fixed in different cases
}


//----------------------------------------------------------------------------------------------------------------------------------
void VCoulombMatrixGenerator::addKinMatEl(mat& A, matElStruct *mES) {

    double b = mES->b;
    int i = mES->i;
    int j = mES->j;
    int l = mES->l;

        if(i == j) {
            A(i,j) += 0.5/(b * b) * (2*i + l + 3./2.);
        } else if(i == (j - 1) ){
            A(i,j) += 0.5/(b * b) * sqrt(j * (j + l + 0.5) );
        } else if (i == (j + 1) ){
            A(i,j) += 0.5/(b * b) * sqrt(i * (i + l + 0.5) );
        }

}

//----------------------------------------------------------------------------------------------------------------------------------
void VCoulombMatrixGenerator::setMatrix(mat& A, int Mdim, double b) {
    IntegratorGaussLegendreGSL integrator;
    //integrator.readTables("../gen_legendre");
    integrator.setOrder(250);


    matElStruct *mES= new matElStruct;
    mES->b= b;
    mES->l=0;


    gsl_function F;
    F.function= Vij;
    F.params= mES;

    gsl_integration_glfixed_table *w= gsl_integration_glfixed_table_alloc(250);

    //set the matrix for the diagonalization: the gsl_function F is the integrand of the matrix element and P is its parameters structure.
    //r_max = 50 seems a good cutoff for Coulomb potential matrix elements up to 80 HO shells - no further investigation has been performed
    //angular momentum has explicitely not taken into account
    
    for( int i = 0; i < Mdim; i++ ) {
        mES->i = i; //set the first function index
        for ( int j = 0; j < Mdim; j++ ) {
            mES->j = j; //set the second function index
            
            //gsl_integration_glfixed has to be changed with the future name of the function
            A(i,j) = integrator.integrate(F, 0., 50);
            //A(i,j)= gsl_integration_glfixed(F, 0., 50., w);

            addKinMatEl(A, mES);
            //gsl_function, r_min, r_max and the workspace array, but we can put it directly inside the integration class
        }
    }

    // Clean
    delete(mES);
    gsl_integration_glfixed_table_free(w);
}
