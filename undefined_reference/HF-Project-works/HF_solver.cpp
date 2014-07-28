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
    int Nparticles = physical_object->N_part;
    h = arma::zeros<mat>(Nbasis,Nbasis);
    Init_rho(Nparticles,Nbasis);
//    rho = (double) Nparticles / (double) Nbasis * arma::eye<mat>(Nbasis,Nbasis);


    e = arma::ones<vec>(Nbasis);

    run_iteration();
}

/*
solver::~solver()
{

}
*/
void solver::Init_rho(int N_p,int Nbasis) { 
  rho = arma::zeros<mat>(Nbasis,Nbasis);
  for(int i=0;i< N_p; i++) 
    rho(i,i)=1;
}

double solver::Get_energy()
{
double Energy = 0;
int    Nbasis = physical_object->N_max;
        for (int a=0;a<Nbasis;a++)
        for (int b=0;b<Nbasis;b++){
                Energy += physical_object->T(a,b)*rho(b,a);
        for (int c=0;c<Nbasis;c++)
        for (int d=0;d<Nbasis;d++)
                Energy += 0.5*physical_object->V(a,b,c,d)*rho(d,b)*rho(c,a); }
return Energy;
}


void solver::PHYSICAL_rho_to_h()
{
    int Nbasis = physical_object->N_max;

    rho.print("rho");

    for(int i1 = 0; i1 < Nbasis; i1++)
    for(int i2 = 0; i2 < Nbasis; i2++)
    {
        h(i1,i2) = physical_object->T(i1,i2);

        for(int j1 = 0; j1 < Nbasis; j1++)
        for(int j2 = 0; j2 < Nbasis; j2++)
            h(i1,i2) += physical_object->V(i1,j2,i2,j1) * rho(j1,j2);
    }

    h.print("h:");
}

void solver::SYSTEM_h_to_rho()
{
    //if (hf->complex_flag)
    //else
    mat D;
    int Nparticles = physical_object->N_part;
    int Nbasis = physical_object->N_max;

    e_old = e;
    rho_old = rho;
    arma::eig_sym(e, D, h);

//     D = D.st();

    cout.precision(11);
    e.print("These are the eigenvalues:");

    D.print("eigenvectors:");

    // You need to keep the N lowest eigenvectors:
    D.resize(Nbasis, Nparticles);
    rho = D*D.t();
}

double solver::Calc_conv_par()
{

/*
 * Old norm:
 * **********
 *
    arma::mat sum_squared = (e + e_old)%(e + e_old);
    arma::mat diff_squared = (e - e_old)%(e - e_old);

    conv_par = diff_squared.max() / sum_squared.max();

 *
 *
 */

    arma::mat sum_rho = rho - rho_old;
    arma::mat diff_rho = rho + rho_old;

    return arma::norm(diff_rho) / arma::norm(sum_rho);
}

void solver::alternative_rho()
{
    rho = (1. - conv_par) * rho + conv_par * rho_old;
}

void solver::run_iteration()
{
    for(unsigned int i = 0; i < max_it; i++)
    {
        std::cout << std::endl << endl<<endl<<endl<< "ITERATION #"<< i << "  "<<endl;

        PHYSICAL_rho_to_h();
        SYSTEM_h_to_rho();

        if(Calc_conv_par() < conv_treshold)
        {
            break;
        }

//        alternative_rho();
    }
}
