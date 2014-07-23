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
    int N_basis = 4;
    int N_particles = 2;
    int system = 0;

    unsigned int max_iterations = 50;
    double conv_treshold = 1e-10;

    Quad quad("quad/glqi_32.bin");

    States states(N_basis, system);

    physical_world Whatever(N_basis, N_particles, V_neutrondrop, T, NULL, states, quad);
    solver HF_solver(&Whatever, max_iterations, conv_treshold);

    return 0;
}








