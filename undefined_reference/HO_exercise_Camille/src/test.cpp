#include "quadrature.hpp"
#include "ho.hpp"
#include <iostream>
#include <iomanip>

void test_various_integrands(Quad& q);
void test_orthonormality(Quad& q);
double f_ho_radial(double,void*); // integrand to check orthonormality of the wfn's
struct f_ho_radial_params {
	int n1,l1,n2,l2;
	double b;
};

int main(int argc, char* argv[] ){
	if (argc!=2){
		std::cerr << "[Error] Expected input file of quadrature knots & weights to be passed on the command line \n" << std::endl;
		exit(-1);
	}
	Quad q(argv[1]); // quadrature instance, call this with q.integrate( double f(double,void*) ) to integrate	
	test_various_integrands(q);
	test_orthonormality(q);
}


/** here we make use of the nice lambda functions introduced in c++11 */
void test_various_integrands(Quad& q){
	std::cout << " (0,inf) exp(-x)               = " << std::setw(10) << std::setprecision(6) << q.integrate( [](double x,void* p){ return exp(-x); },NULL) << " exact = 1 " << std::endl;
	std::cout << " (0,inf) x*exp(-x)             = " << std::setw(10) << std::setprecision(6) << q.integrate( [](double x,void* p){ return x*exp(-x);},NULL) << " exact = 1 " << std::endl;
	std::cout << " (0,inf) 1/6*x^3*exp(-x)       = " << std::setw(10) << std::setprecision(6) << q.integrate( [](double x,void* p){ return 1./6.*x*x*x*exp(-x); },NULL) << " exact = 1 " << std::endl;
	std::cout << " (0,inf) 2./sqrt(pi)*exp(-x^2) = " << std::setw(10) << std::setprecision(6) << q.integrate( [](double x,void* p){ return 2./sqrt(M_PI)*exp(-x*x); },NULL) << " exact = 1 " << std::endl;
	
}

void test_orthonormality(Quad& q){
	struct f_ho_radial_params p;
	p.l1 = 0;
	p.l2 = 0;
	for (int n1=0; n1<10; n1++){
		for (int n2=n1; n2<10; n2++){
			p.n1 = n1;
			p.n2 = n2;
			std::cout << " (n1,l1,n2,l2) = (" << p.n1 << "," << p.l1 << "," << p.n2 << "," << p.l2 << ") = " << q.integrate(f_ho_radial,&p) << std::endl;
		}
	}	
}

double f_ho_radial(double r,void* params){
	struct f_ho_radial_params p = *(struct f_ho_radial_params*) params;
	return r*r*HO::wfn_radial(p.n1,p.l1,r,p.b)*HO::wfn_radial(p.n2,p.l2,r,p.b);
}
