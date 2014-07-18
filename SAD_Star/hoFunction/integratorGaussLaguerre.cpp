#include "integratorGaussLaguerre.h"


//------------------------------------------------------------------------------
IntegratorGaussLaguerre::IntegratorGaussLaguerre(): order_(2), weights_(1e5), abscissa_(1e5)
{
}


//------------------------------------------------------------------------------
IntegratorGaussLaguerre::~IntegratorGaussLaguerre()
{
}


//------------------------------------------------------------------------------
void IntegratorGaussLaguerre::setOrder(int n){
  order_=n;
}


//------------------------------------------------------------------------------
void IntegratorGaussLaguerre::readTables(string weightFile, string abscissaFile, int n){
    this->readTab(weightFile, weights_, n);
    this->readTab(abscissaFile, abscissa_, n);
}


//------------------------------------------------------------------------------
/**
 * @brief IntegratorGaussLaguerre::readTab
 * @param file
 * @param data
 */
void IntegratorGaussLaguerre::readTab(string file, vector<vector<double> > & data, int n){
  ifstream input(file.c_str());
  if(!input){
      throw invalid_argument(string("in ")+__FILE__+", "+__FUNCTION__+", file not found");
  }

  string word;
  vector<double> valVec;
  while(! input.eof()){
      input>>word;
      double val= std::atof(word.c_str());
      valVec.push_back(val);
  } //end while
  data.at(n)= valVec;
}
