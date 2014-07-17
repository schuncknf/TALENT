#ifndef VMATRIXGENERATORC_CPP
#define VMATRIXGENERATORC_CPP

#include<functional>
#include<armadillo>

#include "sphericalhofunc.h"
#include "integratorGaussLegendre.h"



using namespace std;
using namespace arma;
using namespace std::placeholders;


//------------------------------------------------------------------------------
double vR(double r){
    return 1./r;
}


//------------------------------------------------------------------------------
double integrand(double r, int n1, int n2, SphericalHOFunc& rFunc) {
    double res= vR(r);
    res*= rFunc.eval(n1, 0, r);
    res*= rFunc.eval(n2, 0, r);
    res*= r*r;
    return res;
}


//------------------------------------------------------------------------------
double calcElement(int n1, int n2, IntegratorGaussLegendre& integrator, int order, SphericalHOFunc& rFunc){

    double elem=0;
    int l1= 0;
    int l2= 0;
    double b= rFunc.getB();

    // Potential part:
    auto integrandF= bind(integrand, _1, n1, n2, rFunc);
    elem+= integrator.integrate(integrandF, 0., 25);

    // Kinetic part:
    double kin(0);
    if(n1 == n2){
        kin= 0.5*b*b*(2*n1+l1);
    }
    else if(n1== n2 -1){
        kin= 0.5*b*b*sqrt(n2*(n2+l1+0.5));
    }
    else if (n1== n2 + 1){
        kin= 0.5*b*b*sqrt(n1*(n1+l1+0.5));
    }
    elem+= kin;

    return kin;
}


//------------------------------------------------------------------------------
void generateMatrix(mat &A, int nMax, SphericalHOFunc& rFunc) {
    int order=50;
    IntegratorGaussLegendre integrator;
    integrator.readTables("lgvalues-weights.php", "lgvalues-abscissa.php");
    integrator.setOrder(3);

    A.zeros(nMax,nMax);
    for(int i=0; i< nMax; i++){
        for(int j=0; j<nMax; j++){
            A(i,j)= calcElement(i, j, integrator, order, rFunc);
        }
    }
}


//------------------------------------------------------------------------------
int main (int argc, char* argv[]){
  mat A;
  int nMax=100;
  double b=0.5;

  SphericalHOFunc rFunc;
  rFunc.setB(b);

  generateMatrix(A, nMax, rFunc);
  cout<<A<<endl;

  cout<<"Eigen Value= "<<endl;
  vec eigenVal;
  mat eigenVec;
  eig_sym(eigenVal, eigenVec, A);
  cout<< eigenVal<<endl;

  return 0;
}


#endif
