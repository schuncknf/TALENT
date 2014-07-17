#include "integratorGaussLegendre.h"


using namespace std;

//------------------------------------------------------------------------------
double func0(double x){
    return x;
}


//------------------------------------------------------------------------------
bool testInt0(){
  IntegratorGaussLegendre integrator;
  integrator.readTables("lgvalues-weights.php", "lgvalues-abscissa.php");
  double res= integrator.integrateT(func0, 0, 10, 2);
  cout<<res<<endl;
  return abs(res - 50.)<1e-6;
}


//------------------------------------------------------------------------------
int main (int argc, char* argv[]){

    testInt0();

  return 0;
}
