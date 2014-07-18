#include <armadillo>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cassert>
#include <omp.h>





using namespace std;
using namespace arma;

#include "gauss_legendre.h"
#include "harmonic_oscillator.hpp"
#include "functions.hpp"



//This function takes
double kinetic(int nx,int ny,double b)
{
double result=0;
int l=0;
double hw = 1/b/b;
nx = nx-1;
ny = ny-1;
if(nx==ny)
result=.5*hw*(2*nx+1.5);
if(nx-ny==1)
result=.5*hw*sqrt(double(nx)*(double(nx+l)+.5));
if(nx-ny==-1)
result=.5*hw*sqrt(double(ny)*(double(ny+l)+.5));

return result;
}



double wfintegrandinf(double r, double data[])
{
double data1[3];
double data2[3];
data1[0] = data[0];
data1[1] = data[1];
data1[2] = data[2];
data2[0] = data[0+3];
data2[1] = data[1+3];
data2[2] = data[2+3];
double b = data[2];
double denom = cos(r);
denom = denom*denom;
r  = tan(r);
double result = ho_radialtalent(r,data1)*ho_radialtalent(r,data2)*r*r*b*b*b/denom;
//if(abs(result)>1e5)
//cout<<"Denom"<<denom<<"\n";
return result;
}




double meintegrandinf(double r, double data[])
{
double data1[3];
double data2[3];
data1[0] = data[0];
data1[1] = data[1];
data1[2] = data[2];
data2[0] = data[0+3];
data2[1] = data[1+3];
data2[2] = data[2+3];
double b = data[2];
double denom = cos(r);
denom = denom*denom;
r  = tan(r);

return -ho_radialtalent(r,data1)*ho_radialtalent(r,data2)*r*b*b/denom;
}

double tbmeinf(int n1,int n2,double b)
{
int l1=0;
int l2=0;
double twoparams[6];
twoparams[0] = double(n1);
twoparams[1] = double(l1);
twoparams[2] = b;
twoparams[3] = double(n2);
twoparams[4] = double(l2);
twoparams[5] = b;
double result=gauss_legendre(512, meintegrandinf, twoparams, 0., M_PI/2.01);
return result;
}



