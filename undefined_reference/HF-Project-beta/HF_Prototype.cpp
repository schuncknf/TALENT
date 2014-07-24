#include <iostream>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

#include "HF_solver.h"
#include "Physics.cpp"
#include "quadrature.hpp"
#include "States.h"
#include "HF_Phys.h"

int main()
{
    int N_basis = 8;
    int N_particles = 5;
    int system = 5;
    double b  = 0.491043896; // remember b = sqrt(m omega /hbar) hbar omega = 10 MeV, m=938.9059
    double hw = 10;
    
    unsigned int max_iterations = 50;
    double conv_treshold = 1e-10;

    Quad quad("quad/glqi_32.bin");

    States states(N_basis, system);

    physical_world Whatever(N_basis, N_particles, V_neutrondrop, T_full_HO, &hw, states, quad, b);
    solver HF_solver(&Whatever, max_iterations, conv_treshold);

    return 0;
}