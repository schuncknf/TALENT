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
    int N_basis; //= 8;
    int N_particles; //= 2;
    int system; //= 5;
    double b;  //= 0.49104389631014383; // remember b = sqrt(m omega /hbar) hbar omega = 10 MeV, m=938.9059
    double hw; //= 10;
    int simmetry;
//     b = sqrt(.5/20.73 * hw);
//     cout <<b << endl;

    unsigned int max_iterations; //= 2;
    double conv_treshold; //= 1e-10;

    cin >> N_basis >> N_particles >> system >> b >> hw >> max_iterations >> conv_treshold >> simmetry;
    cout << "basis:                      "<< N_basis << endl << "Particles:                  "<< N_particles <<endl 
         << "state codific:              " << system <<endl << "b:                          " << b <<endl 
         << "hw:                         "<< hw <<endl << "Max_iterations:             " << max_iterations << endl
         << "treshold:                   "<< conv_treshold << endl << "simmetry:                   "<<simmetry << endl;
         
         
    Quad glqi = Quad("QLQI_N128.bin"); // global integrator

    States states(N_basis, system);

    physical_world      Whatever(N_basis, N_particles, N_drops_V, T_full_HO, &hw, states, glqi, b, simmetry);                 //0 - no symmetry;   1 - very basic simmetry;    3 - still not imlemented;
    solver              HF_solver(&Whatever, max_iterations, conv_treshold);
    cout << endl<<endl<<" Total energy: " <<HF_solver.Get_energy() << endl;

    return 0;
}
