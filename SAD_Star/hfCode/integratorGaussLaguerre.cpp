#include "integratorGaussLaguerre.h"


//------------------------------------------------------------------------------
IntegratorGaussLaguerre::IntegratorGaussLaguerre(): order_(2)
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
void IntegratorGaussLaguerre::readTables(string tableDir){
    int nMax=5e4;
    for(int n=2; n<nMax; n++){
        ostringstream weightFile;
        weightFile<<tableDir<<"/gauss-laguerre_n"<<n<<"_w.txt";
        ifstream wInput(weightFile.str().c_str());

        ostringstream abscissaFile;
        abscissaFile<<tableDir<<"/gauss-laguerre_n"<<n<<"_x.txt";
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
 * @brief IntegratorGaussLaguerre::readTab
 * @param file
 * @param data
 */
void IntegratorGaussLaguerre::readTab(const string& file, map<int, vector<double> >& data, int n){
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
double IntegratorGaussLaguerre::integrate0ToInf(gsl_function func) const{
    if(weights_.count(order_)<1 || abscissa_.count(order_)<1 ){
        throw logic_error( (string("in ")+__FILE__+" "+__FUNCTION__+", table unavailable for this order").c_str());
    }

    double sum=0;
    double t(0.);
    for(int i=0; i<order_; i++){
        t= abscissa_.at(order_).at(i);
        sum+= exp( log(weights_.at(order_).at(i)) + t) * func.function(t, func.params); //func(t, func.params); GSL_FN_FDF_EVAL_F(&func,t)
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
