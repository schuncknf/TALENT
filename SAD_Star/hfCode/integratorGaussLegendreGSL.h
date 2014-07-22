#ifndef INTEGRATORGAUSSLEGANDREGSL_H
#define INTEGRATORGAUSSLEGANDREGSL_H

#include<string>
#include<iostream>
#include<fstream>
#include<vector>
#include<stdexcept>

#include<gsl/gsl_integration.h>

using namespace std;


class IntegratorGaussLegendreGSL
{

public:
    IntegratorGaussLegendreGSL();
    ~IntegratorGaussLegendreGSL();

    void setOrder(int n);

    double integrate(gsl_function& func, double a, double b) const;

    template<class T>
    double integrate(T& func, double a, double b) const;


//    template<class T>
//    double integrate0ToInf(T func) const;


private:
    int order_;
    gsl_integration_glfixed_table * w_;

};


//------------------------------------------------------------------------------
template<class T>
double IntegratorGaussLegendreGSL::integrate(T& func, double a, double b) const{
    throw logic_error("Not yet implemented for this type of function...");
    return 0.;
}

////------------------------------------------------------------------------------
//template<class T>
//double integratorGaussLegendreGSL::integrate0ToInf(T func) const{
//    func_=func;
//    int dim= 1000;
//    double error= 1e-8;
//    double res;
//    gsl_workspace * w= gsl_integration_workspace_alloc(dim);
//    gsl_integration_qagiu(&gslFunc, 0., error, error, dim, w, &res, &error);

//    return res;
//}

#endif // INTEGRATOR_H
