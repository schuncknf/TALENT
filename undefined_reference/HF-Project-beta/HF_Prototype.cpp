#include <iostream>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

#include"HF_solver.cpp"
#include"Physics.cpp"


int main()
{
double parameters[4] = {1,2,3,4};
Physical_world HF(10, V_coulomb, T_cinetik, Random_rho, parameters);

// in this function you should put the potentials that you want to youse in this code!


while (convergence_check(&HF)) 
{
PHYSICAL_rho_to_h(&HF);
SYSTEM_h_to_rho(&HF);
}

cout << HF.e << endl;

return 0;
}








