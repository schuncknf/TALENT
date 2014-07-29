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

  if(weights_.count(n)<1){
      ostringstream weightFile;
      weightFile<<tableDir_<<"/gauss-legendre_n"<<n<<"_w.txt";
      ifstream wInput(weightFile.str().c_str());

      ostringstream abscissaFile;
      abscissaFile<<tableDir_<<"/gauss-legendre_n"<<n<<"_x.txt";
      ifstream xInput(abscissaFile.str().c_str());

      if(wInput && xInput){
          this->readTab(weightFile.str(), weights_, n);
          this->readTab(abscissaFile.str(), abscissa_, n);
      }

      wInput.close();
      xInput.close();
  }
}



//------------------------------------------------------------------------------
void IntegratorGaussLegendre::setTableDir(string dir){
    tableDir_= dir;
}

//------------------------------------------------------------------------------
void IntegratorGaussLegendre::readTables(string tableDir){
    int nMax=5e4;
    for(int n=2; n<nMax; n++){
        ostringstream weightFile;
        weightFile<<tableDir<<"/gauss-legendre_n"<<n<<"_w.txt";
        ifstream wInput(weightFile.str().c_str());

        ostringstream abscissaFile;
        abscissaFile<<tableDir<<"/gauss-legendre_n"<<n<<"_x.txt";
        ifstream xInput(abscissaFile.str().c_str());

        if(wInput && xInput){
            this->readTab(weightFile.str(), weights_, n);
            this->readTab(abscissaFile.str(), abscissa_, n);
        }

        wInput.close();
        xInput.close();
    }

}


//------------------------------------------------------------------------------
/**
 * @brief IntegratorGaussLegendre::readTab
 * @param file
 * @param data
 */
void IntegratorGaussLegendre::readTab(const string& file, map<int, vector<double> >& data, int n){
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
  data[n]= valVec;
}

//------------------------------------------------------------------------------
double IntegratorGaussLegendre::integrate0ToInf(gsl_function& func) const{
    if(weights_.count(order_)<1 || abscissa_.count(order_)<1 ){
        throw logic_error( (string("in ")+__FILE__+" "+__FUNCTION__+", table unavailable for this order").c_str());
    }

    double sum=0;
    for(int i=0; i<order_; i++){
        double t= abscissa_.at(order_).at(i);
        sum+= weights_.at(order_).at(i) * func.function((1.+ t)/(1.-t), func.params) * 2./(1.-t)/(1.-t);
    }

    return sum;
}

