/*
  compile with:
  g++ -o HF.exe HF.cpp -larmadillo

  A Hartree-Fock Solver for neutron drops
*/

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main()
{
  int nStates = 2;
  double singleParticle;
  mat hamiltonian = zeros(nStates, nStates);
  mat Dmtx = eye(nStates, nStates);
  vec eigenvalues = zeros(nStates, 1);
  vec eigenPrevious = zeros(nStates, 1);
  vec eigenDiff = zeros(nStates, 1);

  cout << "Beep Boop!" << endl;

  return 0;
}
