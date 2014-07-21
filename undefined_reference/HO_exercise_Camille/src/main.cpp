/*
 * Compile with -lgsl -lblas -larmadillo plz
 *
 */

#include "ho.hpp"
#include "quadrature.hpp" /**< don't even go there, very archaic coding style. **/
#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <cassert>
#include <iomanip>

double f_coulomb(double,void*); // the HO_coulomb function to integrate
double T(const int na, const int la, const int nb, const int lb,const double b); // get the analytical kinetic matrix elements
double V(const int na, const int la, const int nb, const int lb,const double b,Quad& q); // get a potential matrix elements

struct V_params{
	int n1,l1,n2,l2;
	double b;
};

double f_coulomb(double r, void* params){
	struct V_params p = *(struct V_params*) params;
	return -r*HO::wfn_radial(p.n1,p.l1,r,p.b)*HO::wfn_radial(p.n2,p.l2,r,p.b);
}
// there is a factor hbar omega in front of every factor, which in my convention becomes, hbar^2 b^2/m, with hbar = m = 1 this becomes just b^2
double T(const int na, const int la, const int nb, const int lb,const double b){
	assert(la==lb);
	if (na==nb){
		return b*b*0.5*(2*na+la+1.5); // N = 2*n + l 
	} else if (na == nb-1) {
		return b*b*0.5*sqrt(nb*(nb+la+0.5));
	} else if (na == nb+1) {
		return b*b*0.5*sqrt(na*(na+la+0.5));
	} else {
		return 0.;
	}
}
/**
 * @param q A quadrature instance to do the integration, see quadrature.hpp
 */
double V(const int na, const int la, const int nb, const int lb,const double b,Quad& q){
	assert(la==lb);
	struct V_params v;
	v.n1 = na;
	v.l1 = la;
	v.n2 = nb;
	v.l2 = lb;
	v.b  = b;
	return q.integrate(f_coulomb,&v); // integrate the f_coulomb function using quadrature
}


int main(int argc, char* argv[]){
	if (argc!=2){
		std::cerr << "[Error] Expected input file of quadrature knots & weights to be passed on the command line \n" << std::endl;
		exit(-1);
	}
	Quad q(argv[1]); // quadrature instance, call this with q.integrate( double f(double,void*) ) to integrate		
	const unsigned NMAX[] = {2,5,10,20,50};
	int l=0;
	for (double b=0.1; b<=4.0; b+=0.1){ // b = sqrt(m omega / hbar )
		std::cout << std::fixed << std::setprecision(2) <<b << " ";
		for (unsigned int n=0; n<5; n++) { // loop over number of basis functions
			unsigned int N = NMAX[n];
			arma::Mat<double> H(N,N); //Hamiltonian matrix, I prefer using the template explicitly
			for (unsigned n1=0; n1<N; n1++){
				for (unsigned n2=n1; n2<N; n2++){
					H(n1,n2) = T(n1,l,n2,l,b) + V(n1,l,n2,l,b,q);
					H(n2,n1) = H(n1,n2); // because all mat elements are real no conjugate needed here
				}
			}
			arma::Col<double> eigval;
			arma::Mat<double> eigvec;
			arma::eig_sym(eigval,eigvec,H);
			std::cout << std::setprecision(9) << eigval(0) << " ";
		}
		std::cout << std::endl;
	}
	return 0;
}

