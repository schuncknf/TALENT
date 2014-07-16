#include "vmatrixgenerator.h"

//------------------------------------------------------------------------------
VMatrixGenerator::VMatrixGenerator() : intOrder_(5)
{
}


//------------------------------------------------------------------------------
VMatrixGenerator::~VMatrixGenerator(){

}


//------------------------------------------------------------------------------
void VMatrixGenerator::setRadialPotential(double (*func)(double)){
    vR_= func;
}


//------------------------------------------------------------------------------
void VMatrixGenerator::setIntOrder(int n){
    intOrder_= n;
}


//------------------------------------------------------------------------------
void VMatrixGenerator::generateMatrix(vector<vector<double> > &mat, int nMax) {
    cout<<"generateMatrix"<<calcElement(1, 2, 3, 4)<<endl;
}


//------------------------------------------------------------------------------
double VMatrixGenerator::calcElement(int n1, int l1, int n2, int l2){
//  double (*myF)(double) = bind(&VMatrixGenerator::integrand, n1, l1, n2, l2, placeholders::_1);
//  auto myF = bind(&VMatrixGenerator::integrand, n1, l1, n2, l2, placeholders::_1);
    n1_= n1;
    n2_= n2;
    l1_= l1;
    l2_= l2;
  return integrator_.integrate3(&VMatrixGenerator::integrand, 0., 1e3, intOrder_);
}


//------------------------------------------------------------------------------
double VMatrixGenerator::integrand(double r) const{

    SphericalHOFunc rFunc;
    double b=0.5;

    double res= vR_(r);
    res*= rFunc.eval(n1_, l1_, b, r);
    res*= rFunc.eval(n2_, l2_, b, r);
    return res;
}
