#include "integratorGaussLegendre.h"


//------------------------------------------------------------------------------
IntegratorGaussLegendre::IntegratorGaussLegendre()
{
}


//------------------------------------------------------------------------------
IntegratorGaussLegendre::~IntegratorGaussLegendre()
{
}


//------------------------------------------------------------------------------
void IntegratorGaussLegendre::readTables(string weightFile, string abscissaFile){
    this->readTab(weightFile, weights_);
    this->readTab(abscissaFile, abscissa_);
}


//------------------------------------------------------------------------------
double IntegratorGaussLegendre::integrate(double (*func)(double), double a, double b, int order) const{
    if(order> weights_.size() || order>abscissa_.size() || order<2){
        throw logic_error( (string("in ")+__FILE__+" "+__FUNCTION__+", table unavailable for this order").c_str());
    }

    double sum=0;
    double x;
    for(int i=0; i<order; i++){
        x= abscissa_.at(order).at(i) * (b-a)/2. + (b+a)/2.;
        sum+= weights_.at(order).at(i)* func(x);
    }
    sum*= (b-a)/2.;

    return sum;
}


//------------------------------------------------------------------------------
double IntegratorGaussLegendre::integrate2(double (*func)(void *, double), void *object, double a, double b, int order){
    if(order> weights_.size() || order>abscissa_.size() || order<2){
        throw logic_error( (string("in ")+__FILE__+" "+__FUNCTION__+", table unavailable for this order").c_str());
    }

    double sum=0;
    double x;
    for(int i=0; i<order; i++){
        x= abscissa_.at(order).at(i) * (b-a)/2. + (b+a)/2.;
        sum+= weights_.at(order).at(i)* func(object, x);
    }
    sum*= (b-a)/2.;

    return sum;
}


//------------------------------------------------------------------------------
double IntegratorGaussLegendre::integrate3(double (*func)(double, int, int), double a, double b, int order, int n1, int n2){
    if(order> weights_.size() || order>abscissa_.size() || order<2){
        throw logic_error( (string("in ")+__FILE__+" "+__FUNCTION__+", table unavailable for this order").c_str());
    }

    double sum=0;
    double x;
    for(int i=0; i<order; i++){
        x= abscissa_.at(order).at(i) * (b-a)/2. + (b+a)/2.;
        sum+= weights_.at(order).at(i)* func(x, n1, n2);
    }
    sum*= (b-a)/2.;

    return sum;
}

//------------------------------------------------------------------------------
/**
 * @brief IntegratorGaussLegendre::readTab
 * @param file
 * @param data
 */
void IntegratorGaussLegendre::readTab(string file, vector<vector<double> > & data){
  ifstream input(file.c_str());
  if(!input){
      throw invalid_argument(string("in ")+__FILE__+", "+__FUNCTION__+", file not found");
  }

  data.resize(0);

  string word;
  int n=1;
  vector<double> valVec;
  while(! input.eof()){
      input>>word;

      if(word[0] == '$'){
          data.push_back(valVec);
          n++;
          valVec.clear();
      }
      else if(word.size() > 40 || word[0] =='0'){
          double val= std::atof(word.c_str());
          valVec.push_back(val);
      }

  } //end while
}
