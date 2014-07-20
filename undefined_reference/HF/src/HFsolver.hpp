/**
 * This file will contain the Hartree Fock solver routine
 */
#ifndef HARTREE_FOCK_SOLVER
#define HARTREE_FOCK_SOLVER
#include <complex>
#include <armadillo>
#include "hamiltonian.hpp"
namespace HartreeFock{
	void run( std::complex<double> (Hamiltonian::*oneBodyPart)(unsigned,unsigned),
                  std::complex<double> (Hamiltonian::*twoBodyPart)(unsigned,unsigned,unsigned,unsigned));

} // namespace HartreeFock
#endif // HARTREE_FOCK_SOLVER
