#ifndef VMATRIXGENERATORC_CPP
#define VMATRIXGENERATORC_CPP

#include<functional>
#include<armadillo>
#include<gsl/gsl_integration.h>

#include "sphericalhofunc.h"
#include "integratorGaussLaguerre.h"
#include "integratorGaussLaguerre.h"



using namespace std;
using namespace arma;
using namespace std::placeholders;



//------------------------------------------------------------------------------
typedef struct{
  int n1;
  int n2;
  SphericalHOFunc rFunc;
} Params;


//------------------------------------------------------------------------------
double vR(double r){
    return -1./r;
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
double integrand2(double r, void* p) {
    Params* param= (Params*) p;

    double res= vR(r);
    res*= param->rFunc.eval(param->n1, 0, r);
    res*= param->rFunc.eval(param->n2, 0, r);
    res*= r*r;
    return res;
}


//------------------------------------------------------------------------------
double calcElement(int n1, int n2, IntegratorGaussLaguerre& integrator, SphericalHOFunc& rFunc){

    double elem=0;
    int l1= 0;
    int l2= 0;
    double b= rFunc.getB();

    // Potential part:
    auto integrandF= bind(integrand, _1, n1, n2, rFunc);
    elem+= integrator.integrate0ToInf(integrandF);

//    int dim= 1000;
//    double error= 1e-2;

//    Params* p= new Params;
//    p->n1= n1;
//    p->n2= n2;
//    p->rFunc= rFunc;

//    double res;
//    gsl_integration_workspace * w= gsl_integration_workspace_alloc(dim);
//    gsl_function F;
//    F.function= integrand2;
//    F.params= p;
//    elem+= gsl_integration_qagiu(&F, 0., error, error, dim, w, &res, &error);

    // Kinetic part:
    double kin(0);
    if(n1 == n2){
        kin= 0.5/b/b*(2*n1+l1 +3./2.);
    }
    else if(n1== n2 -1){
        kin= 0.5/b/b*sqrt(n2*(n2+l1+0.5));
    }
    else if (n1== n2 + 1){
        kin= 0.5/b/b*sqrt(n1*(n1+l1+0.5));
    }
    elem+= kin;

    return elem;

    //clean
    //gsl_integration_workspace_free(w);
}


//------------------------------------------------------------------------------
void generateMatrix(mat &A, int nMax, SphericalHOFunc& rFunc) {
    IntegratorGaussLaguerre integrator;
    integrator.readTables("../gen_laguerre/gauss-laguerre_n100_x.txt", "../gen_laguerre/gauss-laguerre_n100_w.txt", 100);
    integrator.readTables( "../gen_laguerre/gauss-laguerre_n1000_w.txt", "../gen_laguerre/gauss-laguerre_n1000_x.txt", 1000);
    integrator.readTables( "../gen_laguerre/gauss-laguerre_n2000_w.txt", "../gen_laguerre/gauss-laguerre_n2000_x.txt", 2000);
    integrator.setOrder(1000);

    A.zeros(nMax,nMax);
    for(int i=0; i< nMax; i++){
        for(int j=0; j< nMax; j++){
            A(i,j)= calcElement(i, j, integrator, rFunc);
        }
    }
}


//------------------------------------------------------------------------------
int main (int argc, char* argv[]){
  mat A;
  int nMax=30;

  for(double b=0.1; b<2.; b+= 0.1){
      SphericalHOFunc rFunc;
      rFunc.setB(b);
      //  for(double x=0; x<50; x+=1e-2){
      //      cout<<integrand(x, 2, 3, rFunc)<<endl;
      //  }

      generateMatrix(A, nMax, rFunc);
      //cout<<A<<endl;

      //cout<<"Eigen Value= "<<endl;
      vec eigenVal;
      mat eigenVec;
      eig_sym(eigenVal, eigenVec, A);

      double eMin;
      cout<< b<<"  "<<eigenVal(0)<<endl;
      for(int i=0; i<eigenVal.size(); i++){

      }
  }

  return 0;
}


#endif
