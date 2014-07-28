#ifndef HF_SOLVER_H
#define HF_SOLVER_H

#include <armadillo>
#include <cmath>
#include "HF_Phys.h"

class solver
{
    private:
        // Those are the matrix stored (i dont what the main to see those, they can be stored in different way)
        bool complex_flag; // 1 for complex stuff
        physical_world * physical_object;
        double conv_par;
        unsigned int max_it;
        double conv_treshold;


    public:

        // Store the actual density and the hamiltonian elements
        arma::mat rho_old;                // needed for better convergence?
        arma::mat rho;
        arma::mat h;
        arma::vec e_old;
        arma::vec e;

\
        // this function will fill T and V then make a guess for the initial density
        solver(physical_world * pass_object, unsigned int max_iterations, double conv_treshold);

//        ~solver();
        void Init_rho(int, int);
        void PHYSICAL_rho_to_h();
        void SYSTEM_h_to_rho();
        double Get_energy();

        double Calc_conv_par();
        void alternative_rho();
        void run_iteration();
};

#endif
