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
double onebody(int i1,int j2,imat spbasis,double data[])
{
if(i1!=j2)
return 0;
double b = data[0];

int n1 = spbasis(i1,0);
int l1 = spbasis(i1,1);
int s1 = spbasis(i1,2);
int n2 = spbasis(j2,0);
int l2 = spbasis(j2,1);
int s2 = spbasis(j2,2);
double result=0;
double hw = 197.327*197.327/938.9059/b/b;
//This is using Nicolas' parameter
hw = 197.32891*197.32891/938.9059/b/b;
result=hw*(2*n1+1.5);
return result;
}



double tbmewfintegrandinf(double r1,double r2, double data[])
{
double data1[3];
double data2[3];
double data3[3];
double data4[3];
data1[0] = data[0]+1;
data1[1] = data[1];
data1[2] = data[2];
data2[0] = data[0+3*1]+1;
data2[1] = data[1+3*1];
data2[2] = data[2+3*1];
data3[0] = data[0+3*2]+1;
data3[1] = data[1+3*2];
data3[2] = data[2+3*2];
data4[0] = data[0+3*3]+1;
data4[1] = data[1+3*3];
data4[2] = data[2+3*3];
double b = data[2];
double denom1 = cos(r1);
denom1 = denom1*denom1;
r1  = tan(r1);
double denom2 = cos(r2);
denom2 = denom2*denom2;
r2  = tan(r2);
//denom1=1;
//denom2=1;

double V0r= 200.;
double V0s= -91.85;
//double V0t= 91.85;
double kr = 1.487;
double ks = .465;
//double kt = .465;
double Vr = V0r*(exp(-kr*(r1+r2)*(r1+r2))-exp(-kr*(r1-r2)*(r1-r2)))/kr;
double Vs = V0s*(exp(-ks*(r1+r2)*(r1+r2))-exp(-ks*(r1-r2)*(r1-r2)))/ks;

double wf = ho_radialtalent(r1,data1)*ho_radialtalent(r2,data2)*ho_radialtalent(r1,data3)*ho_radialtalent(r2,data4)*r1*r2;
wf = wf+ho_radialtalent(r1,data1)*ho_radialtalent(r2,data2)*ho_radialtalent(r2,data3)*ho_radialtalent(r1,data4)*r1*r2;
wf = wf/2/2;
return -.5*(Vr+Vs)*wf/denom1/denom2;
}




double twobody(int i1,int j2,int k3,int m4,imat spbasis,double data2[])
{
double b = data2[0];
int meshpts=64;
int n1 = spbasis(i1,0);
int l1 = spbasis(i1,1);
int s1 = spbasis(i1,2);
int n2 = spbasis(j2,0);
int l2 = spbasis(j2,1);
int s2 = spbasis(j2,2);
int n3 = spbasis(k3,0);
int l3 = spbasis(k3,1);
int s3 = spbasis(k3,2);
int n4 = spbasis(m4,0);
int l4 = spbasis(m4,1);
int s4 = spbasis(m4,2);





double data[12];
data[0]=double(n1);
data[1]=double(l1);
data[2]=b;
data[0+3*1]=double(n2);
data[1+3*1]=double(l2);
data[2+3*1]=b;
data[0+3*2]=double(n3);
data[1+3*2]=double(l3);
data[2+3*2]=b;
data[0+3*3]=double(n4);
data[1+3*3]=double(l4);
data[2+3*3]=b;



double result=0;
if(s1+s2!=s3+s4)
return result;

if(s1+s2==1||s1+s2==-1)
return result;

if(s1==s3)
result = gauss_legendre_2D_cube(meshpts, tbmewfintegrandinf, data, 0, M_PI/2.001, 0, M_PI/2.001);


if(s1==s4)
result = -gauss_legendre_2D_cube(meshpts, tbmewfintegrandinf, data, 0, M_PI/2.001, 0, M_PI/2.001);


return result;

}


