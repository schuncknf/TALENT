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

#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

// Declarations
namespace VMinnesotaMatrixGenerator{
  typedef vector<vector<vector<vector<double> > > > TwoBodyMat;

  void fillHMatrix(mat& H, mat& density, vector<vector<vector<vector<double> > > >& Vabcd, double b);
  void calc2BodyMat(vector<vector<vector<vector<double> > > >& Vabcd, double b);
}

