#ifndef HF_PHYS_H
#define HF_PHYS_H

#include <armadillo>
#include <cmath>
#include "quadrature.hpp"
#include "States.h"

class physical_world
{
    private:
int simmetry=0;
        struct SPQN{
            unsigned n;
            int sigma; //-1 is spin down, 1 is spin up
            void print() { std::cout<< "[spqn] n=" << n << " sigma = " << sigma << std::endl; }
        };

        // Those are the matrix stored (i dont what the main to see those, they can be stored in different way)
        arma::mat V_internal;
        arma::mat T_internal;
        arma::vec T_internal_diag;
        bool complex_flag; // 1 for complex stuff

        double Test;

    public:
        int N_max;      // Base dimension
        int N_part;     // Number of particles

        // Store the actual density and the hamiltonian elements
        arma::mat rho_old;                // needed for better convergence?
        arma::mat rho;
        arma::mat h;
        arma::vec e_old;
        arma::vec e;

        // This is the interface that the main will see in order to use V and T
        double V(int i1, int i2, int j1, int j2);
        double T(int i1,int j1);
        
        // this function will fill T and V then make a guess for the initial density
        physical_world(int N_basis, int N_particles, double (*fill_V)(int, int, int, int, States&, Quad&, double, void*),
                       double (*fill_T)(int, int, States&, void*), void* parameters, States& states, Quad& quad, double b);
};


#endif
