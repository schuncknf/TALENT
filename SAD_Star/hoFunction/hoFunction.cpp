#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "sphericalhofunc.h"
#include "integratorGaussLegendre.h"

#define HBARRE 6.58211928

using namespace std;


//----------------------------------------------------------------------------------------
double xFunc(double x){
    return x;
}

//----------------------------------------------------------------------------------------
double hoProd(double r){
    SphericalHOFunc func;
    double omega=1.3;
    double m=2.;
    double b= sqrt(m*omega/HBARRE);

    return func.eval(2, 2, b, r) * func.eval(3, 4, b, r);
}

//----------------------------------------------------------------------------------------
/**
 * @brief main
 * @param argc
 * @param argv
 * @return
 * Evaluate the spherical HO function
 */
int main (int argc, char* argv[])
{
    cout<<"eval:"<<endl;
    SphericalHOFunc func;

    double omega=1.3;
    double m=2.;
    double b= sqrt(m*omega/HBARRE);
    double res= func.eval(100, 3, b, 0.5);
    cout<<"res= "<<res<<endl;

    cout<<"Integration test"<<endl;
    IntegratorGaussLegendre integrator;
    integrator.readTables("lgvalues-weights.php", "lgvalues-abscissa.php");

    double integral= integrator.integrate(xFunc, -1, 1., 2);
    cout<<"x=3 "<<xFunc(3.)<<endl;
    cout<<"integra= "<<integral<<endl;
    cout<<"integra= "<<integrator.integrate(xFunc, 0., 10., 2)<<endl;

    cout<<"Integral of the product of two basis functions"<<endl;
    cout<<"int= "<<integrator.integrate(hoProd, 0., 1e6, 5)<<endl;
//    for(double x=0.; x<100; x+=0.01){
//        cout<<x<<" "<<hoProd(x)<<endl;
//    }


    return 0;
}
