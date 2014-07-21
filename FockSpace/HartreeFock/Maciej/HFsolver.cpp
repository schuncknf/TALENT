//HF solver

#include <iostream>
#include <fstream>
#include <cmath>
#include <armadillo>

#include "HFsystem.h"

using namespace std;
using namespace arma;

void HFsolver(mat Density, mat SPH, double precision, int N, int w){
	//w is number of iterations which were done after having reached fixed precision
	
	int k,v=0;
	vec eigenvalold(N);
	
	ofstream plik;
	plik.open("iteration_process.dat");
	
	do{	
		for(int i=0;i<N;i++){
			for(int j=i;j<N;j++){
				SPH(i,j)=0.5*(SPH(i,j)+SPH(j,i)); //symmetrization of single particle hamiltonian
			}
		}
	
		mat eigenvec(N,N); //eigenvector matrix of diagonaized SPH
		vec eigenval(N); //vector of eigenvalues 
		
		eig_sym(eigenval,eigenvec,SPH); //diagonalization of sp hamiltonian
		
		double sum=0;
		for(int i=0;i<N;i++){
			sum+=fabs(eigenvalold(i)-eigenval(i));		
		}
		
		double lambda=0;
		lambda=sum/(double)N;
		
		cout<<"No of interation "<<k<<"   HF energy "<<eigenval(0)<<"  stability parameter is "<<lambda<<endl;
		plik<<k<<"  "<<eigenval(0)<<"   "<<lambda<<endl;
	
		eigenvalold=eigenval;
		
		DensityMatrix(Density,eigenvec,N);
		SingleParticleHamiltonian(SPH,eigenvec,N);
	
		//check of the condition of presicion with the given no of further iterations
		if(lambda<precision) {
			v++;
		}
		k++;
		
	} while(v==w);
	plik.close();

	return;
}
