#include "hamiltonian.hpp"

/**
 * this file contains all the physical details
 * of the Neutrondrop_Swave model.
 * 
 * we will probably want to add a constructor to set
 * all the parameters of the model
*/

std::complex<double> Neutrondrop_Swave::oneBodyPart(unsigned a,unsigned b){
	/* return the one body hamiltonian matrix elements here */
	return std::complex<double>(0.,0.);
}

std::complex<double> Neutrondrop_Swave::twoBodyPart(unsigned a,unsigned b,unsigned c,unsigned d){
	/* return the two body hamiltonian matrix elements here */
	return std::complex<double>(0.,0.);
}	


/** main function here for now.
  * 
  * to decide:
  *  -- each physical system has a main
  * or
  *  -- one file with a main that accesses
  *     the different physical systems...
*/

int main(){
	return 0;
}
