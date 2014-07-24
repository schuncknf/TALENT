#include <armadillo>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cassert>
#include <omp.h>
#include"gnuplot_i.hpp"
//#include kitchensink

using namespace std;
using namespace arma;

#include "gauss_legendre.h"
#include "harmonic_oscillator.hpp"
#include "functions.hpp"
#include "minnesotame.hpp"
#include "hffunctions.hpp"

void plot(Gnuplot& g1,vec X,vec Y,char* fname)
{
            std::vector<double> x = conv_to< std::vector<double> >::from(X);
            std::vector<double> y = conv_to< std::vector<double> >::from(Y);
            g1.set_style("lines").plot_xy(x,y, fname);
}



int main(int argc, char *argv[])
{


if(argc!=4)
{
cout<<"Need the nmaxoccupied,nmax, and hw(MeV)\n";
return 0;

}


int nmaxoccupied=atoi(argv[1]);
int nmax = atoi(argv[2]);
double hw = atof(argv[3]); //MeV
double b = sqrt(197.327*197.327/938.9059/hw);//fm
//This is using Nicolas' parameter
b = sqrt(197.32891*197.32891/938.9059/hw);

int NSHELL = 2*(nmax-1);
cout<<"This means Nshell = "<<NSHELL<<"\n";


double data[1];
data[0] = b;


//spbasis initialization
imat spbasis;
spinit(spbasis,nmax);
//cout<<twobody(0,1,0,1,spbasis,data)<<"\n";



//two particle initialization
imat tpbasis;
tpbasisinit(spbasis,tpbasis,nmax);

//onebody element given in terms of the orbitals referred to in spbasis
mat onebodymat;
onebodymatinit(onebodymat,spbasis,data);


//twobody elements given in terms of combination of orbitals referred to 
//in tpbasis and spbasis
mat twobodymat;
twobodymatinit(twobodymat,spbasis,tpbasis,data);



mat densities(spbasis.n_rows,spbasis.n_rows);
mat fock(spbasis.n_rows,spbasis.n_rows);
densities.fill(0);


//initializes density as diagonal
for(int i=0;i<2*nmaxoccupied;i++)
densities(i,i) = 1;


int m=0;
double enot = enotexpect(densities,onebodymat,twobodymat,spbasis,tpbasis);
double enotold = enot+.1;
cout<<m<<"th Iteration Energy: "<<enot<<"\n";
m++;
while(abs(enot-enotold)>1e-6)
{
enotold = enot;
fockinit(fock,onebodymat,twobodymat,densities,spbasis,tpbasis);
densityupdate(fock,densities,spbasis,nmaxoccupied);
enot = enotexpect(densities,onebodymat,twobodymat,spbasis,tpbasis);
cout<<setiosflags(ios::fixed);
cout<<setprecision(6);	
cout<<m<<"th Iteration Energy: "<<enot<<"\n";
m++;
}
//cout<<densities;






















return 0;
}
