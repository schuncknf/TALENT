#include "integratorGaussLegendreGSL.h"


//------------------------------------------------------------------------------
IntegratorGaussLegendreGSL::IntegratorGaussLegendreGSL(): order_(4){
    w_ = gsl_integration_glfixed_table_alloc(order_);
}


//------------------------------------------------------------------------------
IntegratorGaussLegendreGSL::~IntegratorGaussLegendreGSL()
{
    gsl_integration_glfixed_table_free(w_);
}


//------------------------------------------------------------------------------
void IntegratorGaussLegendreGSL::setOrder(int n){
    gsl_integration_glfixed_table_free(w_);
    order_=n;
    w_= gsl_integration_glfixed_table_alloc( (unsigned int)order_ );
}


//------------------------------------------------------------------------------
double IntegratorGaussLegendreGSL::integrate(gsl_function& func, double a, double b) const{
    return gsl_integration_glfixed (&func, a, b, w_);
}


//------------------------------------------------------------------------------
void IntegratorGaussLegendreGSL::readTables(string tabDir){
    // no need to read table
}

//------------------------------------------------------------------------------
double IntegratorGaussLegendreGSL::integrate0ToInf(gsl_function &func) const{
    //throw logic_error("method to be implemented");

    double sum=0;
    double xi;
    double wi;
    for(int i=0; i<order_; i++){
        gsl_integration_glfixed_point (-1., 1, i, &xi, &wi, w_);
        double t= xi;
        sum+= wi * func.function((1.+ t)/(1.-t), func.params) * 2./(1.-t)/(1.-t);
    }


    return sum;
}
