#include "integratorGaussLegendre.h"


//------------------------------------------------------------------------------
IntegratorGaussLegendre::IntegratorGaussLegendre(): order_(2)
{
}


//------------------------------------------------------------------------------
IntegratorGaussLegendre::~IntegratorGaussLegendre()
{
}


//------------------------------------------------------------------------------
void IntegratorGaussLegendre::setOrder(int n){
  order_=n;
}


//------------------------------------------------------------------------------
void IntegratorGaussLegendre::readTables(string weightFile, string abscissaFile){
    this->readTab(weightFile, weights_);
    this->readTab(abscissaFile, abscissa_);
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
