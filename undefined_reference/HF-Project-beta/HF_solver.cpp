#include <armadillo>
#include <cmath>
#include <iostream>

#include "HF_solver.h"
#include "HF_Phys.h"

solver::solver(physical_world * pass_object, unsigned int max_iterations, double treshold)
{
    // physical_object = new physical_world(pass_object);
    physical_object = pass_object;
    max_it = max_iterations;
    conv_treshold = treshold;

    int Nbasis = physical_object->N_max;

    h = arma::zeros<mat>(Nbasis,Nbasis);
    rho = arma::eye<mat>(Nbasis,Nbasis);
    e = arma::ones<vec>(Nbasis);

    run_iteration();
}

/*
solver::~solver()
{

}
*/

void solver::PHYSICAL_rho_to_h()
{
    int Nbasis = physical_object->N_max;

    for(int i1 = 0; i1 < Nbasis; i1++)
    for(int i2 = 0; i2 < Nbasis; i2++)
    {
        h(i1,i2) = physical_object->T(i1,i2);

        for(int j1 = 0; j1 < Nbasis; j1++)
        for(int j2 = 0; j2 < Nbasis; j2++)
            h(i1,i2) += physical_object->V(i1,j2,i2,j1) * rho(j1,j2);
    }
}

void solver::SYSTEM_h_to_rho()
{
    //if (hf->complex_flag)
    //else
    arma::mat D;

    int Nparticles = physical_object->N_part;
    int Nbasis = physical_object->N_max;

    e_old = e;
    rho_old = rho;
    arma::eig_sym(e, D, h);

    D = D.st();

    e.print("These are the eigenvalues:");

    D.print("eigenvectors:");

    // You need to keep the N lowest eigenvectors:
    // (are they automatically ordered using the eig_sym function?)
    D.resize(Nbasis, Nparticles);

    /* I think it's better to use transpose (= hermitian conjugation if complex
     * and then multiply
     * 
     * alternative code:
     * 

    arma::mat Dstar = conj(D);

    for(int i1; i1 < N; i1++)
    for(int i2; i2 < N; i2++)
    for(int i3; i3 < N; i3++)
        rho(i1,i2) = D(i1,i3)*Dstar(i2,i3);

     *
     *
     */

    rho = D*D.t();
}

double solver::Calc_conv_par()
{
    arma::mat sum_squared = (e + e_old)%(e + e_old);
    arma::mat diff_squared = (e - e_old)%(e - e_old);

    conv_par = diff_squared.max() / sum_squared.max();

    return conv_par;
}

void solver::alternative_rho()
{
    rho = (1. - conv_par) * rho + conv_par * rho_old;
}

void solver::run_iteration()
{
    for(unsigned int i = 0; i < max_it; i++)
    {
        std::cout << std::endl << i;

        PHYSICAL_rho_to_h();
        SYSTEM_h_to_rho();

        if(Calc_conv_par() < conv_treshold)
        {
            break;
        }

//        alternative_rho();
    }
}
