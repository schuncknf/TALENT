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
#include "minnesotame.hpp"
#include "hffunctions.hpp"

void spinit(imat& spbasis,int nmax)
{

for(int i=0;i<nmax;i++)
		{
			
			for(int spin=-1;spin<=1;spin+=2)
			{
			imat temp1(1,3);
			temp1(0,0) = i;
			temp1(0,1) = 0;//this is the temp holder for l;
			temp1(0,2) = spin;
			//will add holder for occupation when moving forward
			spbasis = join_cols(spbasis,temp1);

			}

		}

}


void tpbasisinit(imat& spbasis,imat& tpbasis,int nmax)
{
int indexcount=0;
for(int i=0;i<spbasis.n_rows;i++)
	for(int j=0;j<spbasis.n_rows;j++)
			{
			//only doing spin=0 channel
			int spin1 = spbasis(i,2);
			int spin2 = spbasis(j,2);
			if(spin1+spin2==0)
				{
					imat temp2(1,2);
					temp2(0,0) = i;
					temp2(0,1) = j;
					tpbasis = join_cols(tpbasis,temp2);
					indexcount++;
				}
		

			}

}


void onebodymatinit(mat& onebodymat,imat& spbasis,double data[])
{
onebodymat.resize(spbasis.n_rows,spbasis.n_rows);
onebodymat.fill(0);
for(int i=0;i<spbasis.n_rows;i++)
	for(int j=0;j<spbasis.n_rows;j++)
		onebodymat(i,j) = onebody(i,j,spbasis,data);
}

void twobodymatinit(mat& twobodyme,imat& spbasis, imat& tpbasis,double data[])
{
imat tpindex(spbasis.n_rows,spbasis.n_rows);
tpindex.fill(-1);
for(int i=0;i<tpbasis.n_rows;i++)
{
	int i1 = tpbasis(i,0);
	int j2 = tpbasis(i,1);
	tpindex(i1,j2) = i;
}


twobodyme.resize(tpbasis.n_rows,tpbasis.n_rows);
twobodyme.fill(0);
#pragma omp parallel for
for(int i=0;i<tpbasis.n_rows;i++)
	for(int j=i;j<tpbasis.n_rows;j++)
		{
		int i1 = tpbasis(i,0);
		int j2 = tpbasis(i,1);
		int k3 = tpbasis(j,0);
		int l4 = tpbasis(j,1);
		int tp1 = tpindex(j2,i1);
		int tp2 = tpindex(l4,k3);

		if(i1<j2&&k3<l4)
		{
			twobodyme(i,j) = twobody(i1,j2,k3,l4,spbasis,data);
			twobodyme(tp1,tp2) = twobodyme(i,j);
			twobodyme(tp1,j) = -twobodyme(i,j);
			twobodyme(i,tp2) = -twobodyme(i,j);
			twobodyme(j,i) = twobodyme(i,j);
			twobodyme(tp2,tp1) = twobodyme(i,j);
			twobodyme(j,tp1) = -twobodyme(i,j);
			twobodyme(tp2,i) = -twobodyme(i,j);
		
		}
		}
}
double enotexpect(mat& densities,mat& onebodymat,mat& twobodymat,imat& spbasis,imat& tpbasis)
{
double result = 0;
for(int i=0;i<spbasis.n_rows;i++)
	for(int j=0;j<spbasis.n_rows;j++)
		result+= onebodymat(i,j)*densities(j,i);

for(int i=0;i<tpbasis.n_rows;i++)
	for(int j=0;j<tpbasis.n_rows;j++)
		{
		int i1 = tpbasis(i,0);
		int j2 = tpbasis(i,1);
		int k3 = tpbasis(j,0);
		int l4 = tpbasis(j,1);
		result+= .5*twobodymat(i,j)*densities(j2,l4)*densities(i1,k3);
		
		}
return result;

}
void fockinit(mat& fock,mat& onebodymat,mat& twobodymat,mat& densities,imat& spbasis,imat& tpbasis)
{

fock = onebodymat;
for(int i=0;i<tpbasis.n_rows;i++)
	for(int j=0;j<tpbasis.n_rows;j++)
		{
		int i1 = tpbasis(i,0);
		int j2 = tpbasis(i,1);
		int k3 = tpbasis(j,0);
		int l4 = tpbasis(j,1);
		fock(i1,k3) += twobodymat(i,j)*densities(j2,l4);
		}


}
void densityupdate(mat& fock,mat& densities,imat& spbasis,int nmaxoccupied)
{
vec eigval;
mat eigvec;
spindiagonalizeascending(eigval,eigvec,fock,spbasis);

mat temp1 = eigvec;
temp1.resize(spbasis.n_rows,2*nmaxoccupied);

//densities = (.7*densities + .3*temp1*temp1.t());
densities = temp1*temp1.t();

}

