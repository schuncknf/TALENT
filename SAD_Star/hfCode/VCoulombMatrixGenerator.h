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

#include <armadillo>

using namespace arma;

// Declarations
namespace VCoulombMatrixGenerator{

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


  double potential(double r, void *params);
  double Vij(double r, void *params);
  void addKinMatEl(mat& A, matElStruct *mES);
  void setMatrix(mat& A, int Mdim, double b);
}

