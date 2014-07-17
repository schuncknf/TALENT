#include "sphericalhofunc.h"
#include "integratorGaussLegendre.h"
#include <armadillo>


using namespace std;
using namespace arma;


//------------------------------------------------------------------------------
double vR(double r){
    return 1./r;
}


//------------------------------------------------------------------------------
double integrand(double r, int n1, int n2) {
    SphericalHOFunc rFunc;
    double b=10;

    double res= vR(r);
    res*= rFunc.eval(n1, 0, b, r);
    res*= rFunc.eval(n2, 0, b, r);
    res*= r*r;
    return res;
}


//------------------------------------------------------------------------------
double calcElement(int n1, int n2, IntegratorGaussLegendre& integrator, int order){

    return integrator.integrate3(&integrand, 0., 1e3, order, n1, n2);
}


//------------------------------------------------------------------------------
void generateMatrix(mat &A, int nMax) {
    int order=20;
    IntegratorGaussLegendre integrator;
    integrator.readTables("lgvalues-weights.php", "lgvalues-abscissa.php");

    A.zeros(nMax,nMax);
    for(int i=0; i< nMax; i++){
        for(int j=0; j<nMax; j++){
            A(i,j)= calcElement(i, j, integrator, order);
        }
    }
}


//------------------------------------------------------------------------------
int main (int argc, char* argv[]){
  mat A;
  int nMax=5;

  generateMatrix(A, nMax);
  cout<<A<<endl;
}
