#include "integratorGaussLegendreGSL.h"


//------------------------------------------------------------------------------
IntegratorGaussLegendreGSL::IntegratorGaussLegendreGSL(): order_(2){
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
    w_= gsl_integration_glfixed_table_alloc(n);
}


//------------------------------------------------------------------------------
double IntegratorGaussLegendreGSL::integrate(gsl_function& func, double a, double b) const{
    return gsl_integration_glfixed (&func, a, b, w_);
}
