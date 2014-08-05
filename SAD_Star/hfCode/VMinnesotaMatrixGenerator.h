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

#include "constants.h"

using namespace std;
using namespace arma;

// Declarations
namespace VMinnesotaMatrixGenerator{
  typedef vector<vector<vector<vector<double> > > > FourIndiceMat;
  typedef vector<vector< FourIndiceMat > > TwoBodyMat;


  TwoBodyMat emptyTwoBodyMat(int NMax, int lMax);
  FourIndiceMat emptyFourIndiceMat(int dim);

  void fillHMatrix(vector<mat>& H, vector<mat>& density, TwoBodyMat& Vabcd, double& b, int lMax);
  void calc2BodyMat(TwoBodyMat& Vabcd, double& b, int order, int NMax);
}

