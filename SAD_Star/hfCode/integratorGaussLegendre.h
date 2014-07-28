#ifndef INTEGRATORGAUSSLEGENDRE_H
#define INTEGRATORGAUSSLEGENDRE_H

#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<stdexcept>
#include<map>

#include<gsl/gsl_integration.h>

using namespace std;


class IntegratorGaussLegendre
{

public:
    IntegratorGaussLegendre();
    ~IntegratorGaussLegendre();

    void readTables(string tabDir);
    void setOrder(int n);

    template<class T>
    double integrate(T& func, double a, double b) const;
    template<class T>
    double integrate0ToInf(T& func) const;
    double integrate0ToInf(gsl_function& F) const;


private:
    int order_;
    map<int, vector<double> > weights_;
    map<int, vector<double> > abscissa_;
    void readTab(const string& file, map<int, vector<double> >& data, int n);
};


//------------------------------------------------------------------------------
template<class T>
double IntegratorGaussLegendre::integrate(T& func, double a, double b) const{
    if(weights_.count(order_)<1 || abscissa_.count(order_)<1 ){
        throw logic_error( (string("in ")+__FILE__+" "+__FUNCTION__+", table unavailable for this order").c_str());
    }

    double sum=0;
    double x;
    for(int i=0; i<order_; i++){
        x= abscissa_.at(order_).at(i) * (b-a)/2. + (b+a)/2.;
        sum+= weights_.at(order_).at(i)* func(x);
    }
    sum*= (b-a)/2.;

    return sum;
}


//------------------------------------------------------------------------------
template<class T>
double IntegratorGaussLegendre::integrate0ToInf(T& func) const{
    if(weights_.count(order_)<1 || abscissa_.count(order_)<1 ){
        throw logic_error( (string("in ")+__FILE__+" "+__FUNCTION__+", table unavailable for this order").c_str());
    }

    double sum=0;
    for(int i=0; i<order_; i++){
        double t= abscissa_.at(order_).at(i);
        sum+= weights_.at(order_).at(i) * func( (1.+ t)/(1.-t) ) * 2./(1.-t)/(1.-t);
    }

    return sum;
}

#endif // INTEGRATOR_H
