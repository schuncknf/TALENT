#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include<string>
#include<iostream>
#include<fstream>
#include<vector>
#include<stdexcept>

using namespace std;


class IntegratorGaussLegendre
{

    friend class VMatrixGenerator;

public:
    IntegratorGaussLegendre();
    ~IntegratorGaussLegendre();

    void readTables(string weightFile, string abscissaFile);
    void setOrder(int n);

    template<class T>
    double integrate(T func, double a, double b) const;


private:
    int order_;
    vector<vector<double> > weights_;
    vector<vector<double> > abscissa_;
    void readTab(string file, vector<vector<double> >& data);

};


//------------------------------------------------------------------------------
template<class T>
double IntegratorGaussLegendre::integrate(T func, double a, double b) const{
    if(order_> weights_.size() || order_>abscissa_.size() || order_<2){
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

#endif // INTEGRATOR_H
