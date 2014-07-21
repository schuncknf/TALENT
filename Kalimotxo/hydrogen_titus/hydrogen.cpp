#include <armadillo>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cassert>
#include <omp.h>
#include"gnuplot_i.hpp"




using namespace std;
using namespace arma;

#include "gauss_legendre.h"
#include "harmonic_oscillator.hpp"
#include "functions.hpp"
void plot(Gnuplot& g1,vec X,vec Y,char* fname)
{
            std::vector<double> x = conv_to< std::vector<double> >::from(X);
            std::vector<double> y = conv_to< std::vector<double> >::from(Y);
            g1.set_style("lines").plot_xy(x,y, fname);
}



int main(int argc, char *argv[])
{

//oscillator paramter in bohr
double b = 1;
int nmax=20;
if(argc!=4)
{
cout<<"Need 3 inputs, b, n1, and n2 for orthogonormality check.\n";
return 0;
}
b = double(atof(argv[1]));
int n1= int(atof(argv[2])+.1);
int l1=0;
int n2= int(atof(argv[3])+.1);
int l2=0;


///////////////////checks orthogonality///////////////////////////

double twoparams[6];
twoparams[0] = double(n1);
twoparams[1] = double(l1);
twoparams[2] = b;
twoparams[3] = double(n2);
twoparams[4] = double(l2);
twoparams[5] = b;
double result=0;
int meshwf = 1024;
result += gauss_legendre(meshwf, wfintegrandinf, twoparams, 0., M_PI/8);
result += gauss_legendre(meshwf, wfintegrandinf, twoparams, M_PI/8, M_PI/4);
result += gauss_legendre(meshwf, wfintegrandinf, twoparams, M_PI/4, 3*M_PI/8);
result += gauss_legendre(meshwf, wfintegrandinf, twoparams, 3*M_PI/8, M_PI/2.02);


if(n1==n2)
result = abs(1-result);
if(n1!=n2)
result = abs(result);
cout<<"Error from expected norm: "<<result<<"\n";



/////////////////////Calculates ground state for a variety of nmaxes and oscillator parameters///////////

double blimit = 2;
int numbs = 20;
int nmaxlim = 40;


///////////////////Establishes gnuplot output via the header file/////////////////
mat results;
vec X2(numbs);
vec Y2(numbs);
Gnuplot g2("lines");
g2.set_xrange(0,blimit);
g2.set_yrange(-.6,-.3);
g2.showonscreen();
g2.cmd("set arrow from 0,-.5 to 2,-.5");
Gnuplot g3;
g3.set_xrange(0,blimit);
g3.set_yrange(-.6,-.3);
g3.cmd("set arrow from 0,-.5 to 2,-.5");
g3.savetops("gsvariation");
Gnuplot g4;
g4.set_xrange(0,blimit);
g4.set_yrange(-.2,0);
g4.cmd("set arrow from 0,-.125 to 2,-.125");
g4.savetops("excitedvariation");
vec Y3(numbs);
char graphtitle[30];
///////////////////////////////////////////////////////////////////////////////
ivec nvalues(5);
nvalues(0) =  2;
nvalues(1) =  5;
nvalues(2) = 10;
nvalues(3) = 20;
nvalues(4) = 50;


for(int m=0;m<nvalues.n_rows;m++)
{
	int l=nvalues(m);
	X2.fill(0);
	Y2.fill(0);
	Y3.fill(0);
	for(int k=1;k<=numbs;k++)
		{
		b = double(k)*blimit/double(numbs);
		X2(k-1)=b;
		mat H(l,l);
		H.fill(0);
		#pragma omp parallel for
		for(int i=0;i<l;i++)
			for(int j=0;j<l;j++)
				{
			
				H(i,j) = kinetic(i+1,j+1,b)+tbmeinf(i+1,j+1,b);
				}	
		vec eigval;
		mat eigvec;
		eig_sym( eigval, eigvec, H ); 
		Y2(k-1) = eigval(0);
		Y3(k-1) = eigval(1);
		}
	cout<<"Nmax = "<<l<<" is done\n";
	sprintf(graphtitle,"GS Nmax= %d",l);
	plot(g2,X2,Y2,graphtitle);
	plot(g3,X2,Y2,graphtitle);
	plot(g4,X2,Y3,graphtitle);
	if(m==0)
	results = join_cols(results,X2);
	results = join_rows(results,Y2);
	
}

ofstream resultfile;
resultfile.open ("bnmaxresults.dat");
resultfile<<setiosflags(ios::fixed);
cout<<setiosflags(ios::fixed);	
for(int i=0;i<results.n_rows;i++)
	for(int j=0;j<results.n_cols;j++)
		{
			
			if(j==0)
			{
				cout<<setprecision(1);
				if(i==0)
				cout<<results(i,j);
				else
				cout<<"\n"<<results(i,j);
			
			}			
			else
			{
			cout<<setprecision(9);
			cout<<" "<<results(i,j);

			}

			if(j==0)
			{
				resultfile<<setprecision(1);
				if(i==0)
				resultfile<<results(i,j);
				else
				resultfile<<"\n"<<results(i,j);
			
			}			
			else
			{
			resultfile<<setprecision(9);
			resultfile<<" "<<results(i,j);

			}
		}

  resultfile.close();



cout<<"\nGraph Saved under gsvariation.ps and excitedvariation.\n";
/////////////////////////////////////////////////////////////////////////////////



sleep(1); //gives time to get output to gnuplot before closing
return 0;
}
