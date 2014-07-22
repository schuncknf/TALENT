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
    potStruct *pS;
    int i; //main quantum number for the first function
    int j; //main quantum number for the second function
    int l; //angular momentum for the first function
    int l2; //angular momentum for the second function
    double b; //frequency dependent variable
} matElStruct;

/*
    FUNCTIONS
*/
//----------------------------------------------------------------------------------------------------------------------------------
double potential(double r, void *params) {
    //Coulomb potential in natural units

    potStruct *pS = (potStruct *)params;
    
    return -(1./r);
}

//----------------------------------------------------------------------------------------------------------------------------------
double Vij(double r, void *params) {
    //generates the integrand R_nl(r)*V(r)R_n'l'(r) using the function potential
    
    prms *p = (prms *)params; //casting of the params pointer to a prms structure
    
    int i = p->i; 
    int j = p->j;
    int l = p->l;
    double b = p->b;
    double y, z;
    
    //the name of the function generate_basis has to be substitute with the one we will choose
    y = generate_basis(r, i, l, b); //obtain the value of the basis function for the different values of i and j
    z = generate_basis(r, j, l, b);
    
    return y * z * r * r * potential(r, p->pS); //r^2 because of the jacobian of spherical coordinates
    //explicitely ignoring spherical harmonics in case of l=0, has to be fixed in different cases
}

//----------------------------------------------------------------------------------------------------------------------------------
void setMatrix(mat *A, int Mdim, gsl_function *F) {
    //set the matrix for the diagonalization: the gsl_function F is the integrand of the matrix element and P is its parameters structure.
    //r_max = 50 seems a good cutoff for Coulomb potential matrix elements up to 80 HO shells - no further investigation has been performed
    //angular momentum has explicitely not taken into account
    
    for( int i = 0; i < Mdim; i++ ) {
        P->i = i; //set the first function index
        for ( int j = 0; j < Mdim; j++ ) {
            P->j = j; //set the second function index
            
            //gsl_integration_glfixed has to be changed with the future name of the function
            A(i,j) = gsl_integration_glfixed (F, 0., 50, w);
            //gsl_function, r_min, r_max and the workspace array, but we can put it directly inside the integration class
        }
    }

}
