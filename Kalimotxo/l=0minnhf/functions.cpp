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

void spindiagonalizeascending(vec& eigval, mat& eigvec,mat& obdens,imat& spbasis)
{
mat temp(spbasis.n_rows/2,spbasis.n_rows/2);
eigvec.resize(spbasis.n_rows,spbasis.n_rows);
eigval.resize(spbasis.n_rows);
temp.fill(0);
for(int i=0;i<spbasis.n_rows/2;i++)
for(int j=0;j<spbasis.n_rows/2;j++)
temp(i,j) = obdens(2*i,2*j);

vec eigvaltemp;
mat eigvectemp;
eig_sym( eigvaltemp, eigvectemp, temp );


for(int i=0;i<spbasis.n_rows/2;i++)
{
eigval(2*i) = eigvaltemp(i);
eigval(2*i+1) = eigvaltemp(i);
for(int j=0;j<spbasis.n_rows/2;j++)
{
eigvec(2*i,2*j) = eigvectemp(i,j);
eigvec(2*i+1,2*j+1) = eigvectemp(i,j);
}


}




}
void spindiagonalizedescending(vec& eigval, mat& eigvec,mat& obdens,imat& spbasis)
{
mat temp(spbasis.n_rows/2,spbasis.n_rows/2);
eigvec.resize(spbasis.n_rows,spbasis.n_rows);
eigval.resize(spbasis.n_rows);
temp.fill(0);
for(int i=0;i<spbasis.n_rows/2;i++)
for(int j=0;j<spbasis.n_rows/2;j++)
temp(i,j) = obdens(2*i,2*j);

vec eigvaltemp;
mat eigvectemp;
eig_sym( eigvaltemp, eigvectemp, temp );


for(int i=0;i<spbasis.n_rows/2;i++)
{
eigval(spbasis.n_rows - 2*i) = eigvaltemp(i);
eigval(spbasis.n_rows-(2*i+1)) = eigvaltemp(i);
for(int j=0;j<spbasis.n_rows/2;j++)
{
eigvec(spbasis.n_rows-2*i,spbasis.n_rows-2*j) = eigvectemp(i,j);
eigvec(spbasis.n_rows-(2*i+1),spbasis.n_rows-(2*j+1)) = eigvectemp(i,j);
}


}




}
