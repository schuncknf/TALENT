#ifndef VMATRIXGENERATOR_H
#define VMATRIXGENERATOR_H

#include <functional>

#include "sphericalhofunc.h"
#include "integratorGaussLegendre.h"


using namespace std;

class IntegratorGaussLegendre;

class VMatrixGenerator
{
    friend class IntegratorGaussLegendre;
public:
    VMatrixGenerator();
    ~VMatrixGenerator();

    void setRadialPotential(double (*func)(double));
    void setIntOrder(int n);
    void generateMatrix(vector<vector<double> >& mat, int nMax) ;

private:
    double (*vR_)(double);
    IntegratorGaussLegendre integrator_;
    int intOrder_;
    int n1_;
    int n2_;
    int l1_;
    int l2_;

    double calcElement(int n1, int l1, int n2, int l2) ;
    double integrand(double r) const;
};

#endif // VMATRIXGENERATOR_H
