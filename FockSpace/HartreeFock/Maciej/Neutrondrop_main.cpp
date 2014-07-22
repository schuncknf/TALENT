#include <iostream>
#include <fstream>
#include <cmath>
#include <armadillo>

#include "HFsystem.h"
#include "HFsolver.h"

using namespace std;
using namespace arma;

int main() {
	
	int N=10;
	
	mat Density(N,N); //density matrix in HO basis
	mat SPH(N,N); //single-particle hamiltonian 
	mat eigenvec(N,N); //eigenvector matrix of diagonaized SPH
	vec eigenval(N); //vector of eigenvalues 
	
	
	return 0;
}
