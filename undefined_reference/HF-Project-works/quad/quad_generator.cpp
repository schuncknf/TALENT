#include "quadrature.hpp"
#include <iostream>

int main(){
	int kind,nt;
	double alpha,beta;
	char file[128]; // I don't expect a file name longer than ... characters...
	std::cout << " Kind of quadrature? (5= Gauss Laguerre) " << std::endl;
	std::cin >> kind;
	std::cout << " alpha? (0 for plain Gauss Laguerre) " << std::endl;
	std::cin >> alpha;
	std::cout << " beta? (0 if you don't need this) " << std::endl;
	std::cin >> beta;
	std::cout << " Number of knots? " << std::endl;
	std::cin >> nt;
	std::cout << " Output file? " << std::endl;
	std::cin >> file;
	Quad q(kind,alpha,beta,nt); // make a Quad instance	
	q.write(file); // write weights and knots
	return 0;
}
