#ifndef INTEGRATORGAUSSLAGUERRE_H
#define INTEGRATORGAUSSLAGUERRE_H

#include<string>
#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<stdexcept>
#include<cmath>
#include <math.h>
#include <map>

using namespace std;

class IntegratorGaussLaguerre
{

public:
    IntegratorGaussLaguerre();
    ~IntegratorGaussLaguerre();

    void readTables(string tabDir);
    void setOrder(int n);

//    template<class T>
//    double integrate(T func, double a, double b) const;
    template<class T>
    double integrate0ToInf(T func) const;


private:
    int order_;
    map<int, vector<double> > weights_;
    map<int, vector<double> > abscissa_;
    void readTab(const string& file, map<int, vector<double> >& data, int n);

};


////------------------------------------------------------------------------------
//template<class T>
//double IntegratorGaussLaguerre::integrate(T func, double a, double b) const{
//    if(order_> weights_.size() || order_>abscissa_.size() || order_<2){
//        throw logic_error( (string("in ")+__FILE__+" "+__FUNCTION__+", table unavailable for this order").c_str());
//    }

//    double sum=0;
//    double x;
//    for(int i=0; i<order_; i++){
//        x= abscissa_.at(order_).at(i) * (b-a)/2. + (b+a)/2.;
//        sum+= weights_.at(order_).at(i)* func(x);
//    }
//    sum*= (b-a)/2.;

//    return sum;
//}


//------------------------------------------------------------------------------
template<class T>
double IntegratorGaussLaguerre::integrate0ToInf(T func) const{
    if(weights_.count(order_)<1 || abscissa_.count(order_)<1 ){
        throw logic_error( (string("in ")+__FILE__+" "+__FUNCTION__+", table unavailable for this order").c_str());
    }

    double sum=0;
    double t(0.);
    for(int i=0; i<order_; i++){
        t= abscissa_.at(order_).at(i);
        sum+= exp( log(weights_.at(order_).at(i)) + t) * func(t);
        //cout<<"exp "<< exp(t)<<"  func "<<func(t)<<"w "<<weights_.at(order_).at(i)<<endl;
//        cout<<"g-laguerre xi= "<<abscissa_.at(order_).at(i)<<endl;
//        cout<<"wi= "<<weights_.at(order_).at(i)<<endl;
        if(std::isnan(sum)){
            cout<<"STOP"<<endl;
            break;
        }

    }

    return sum;
}

#endif // INTEGRATOR_H
