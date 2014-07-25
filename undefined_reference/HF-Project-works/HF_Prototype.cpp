#include <iostream>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

#include "quadrature.hpp"
#include "HF_solver.h"
#include "Physics.cpp"
#include "States.h"
#include "HF_Phys.h"


/*
 * 
 * Its time to clean up the code and start to implement the potential with spin/ isospin and angular momentum in the same way of the V_neutron_drop
 * 
 * Then start to change physics and states to take into account the new symmetries (maybe a envarment variable to choose the used symmetry?)
 * 
 */



int main()
{
    int N_basis = 10;
    int N_particles = 2;
    int system = 5;
    double b  = 0.49104389631014383; // remember b = sqrt(m omega /hbar) hbar omega = 10 MeV, m=938.9059
    double hw = 10;
    b = sqrt(.5/20.73 * hw);
    cout <<b << endl;

    unsigned int max_iterations = 50;
    double conv_treshold = 1e-10;


    Quad glqi = Quad("QLQI_N128.bin"); // global integrator
    
    States states(N_basis, system);

    physical_world Whatever(N_basis, N_particles, N_drops_V, T_full_HO, &hw, states, glqi, b);
    solver HF_solver(&Whatever, max_iterations, conv_treshold);
    cout << endl<<endl<<" Total energy: " <<HF_solver.Get_energy() << endl;

    return 0;
}