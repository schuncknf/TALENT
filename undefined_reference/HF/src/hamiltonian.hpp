#ifndef HAMILTONIAN_HPP
#define HAMILTONIAN_HPP
#include <complex>
/**
 *
 * the rough idea is that all that the HF
 * solver needs to have access to is 
 * one body part and two body matrix
 * elements so we make a virtual base 
 * class with these methods.
 * Implementation wise the HF solver 
 * sees nothing of the physical system.
 */
class Hamiltonian{
	virtual std::complex<double> oneBodyPart(unsigned,unsigned) = 0; // pure virtual functions, any non-virtual child classes should implement these
	virtual std::complex<double> twoBodyPart(unsigned,unsigned,unsigned,unsigned) = 0;	
};



/***********************************************
 *  List of all the physical systems           *
 ***********************************************/

class Neutrondrop_Swave : public Hamiltonian {
	std::complex<double> oneBodyPart(unsigned,unsigned);
	std::complex<double> twoBodyPart(unsigned,unsigned,unsigned,unsigned);	
};
#endif // HAMILTONIAN_HPP
