#include <armadillo>
#include <cmath>
#include <iostream>

#include "HF_solver.h"
#include "HF_Phys.h"

#include "Flags.h"

solver::solver(physical_world * pass_object, unsigned int max_iterations, double treshold,int N_n, int N_p)
{
    // physical_object = new physical_world(pass_object);
    physical_object = pass_object;
    max_it = max_iterations;
    conv_treshold = treshold;
    Number_of_neutrons = N_n;
    Number_of_protons = N_p;
    Nbasis = physical_object->N_max;
    Nparticles = physical_object->N_part;
    h = arma::zeros<mat>(Nbasis,Nbasis);
    Init_rho();
    e = arma::ones<vec>(Nbasis);

    run_iteration();
}

/*
solver::~solver()
{

}
*/
void solver::Init_rho() {
  rho = arma::zeros<mat>(Nbasis,Nbasis);
  for(int i=0;i< Nparticles; i++)
    rho(i,i)=1;
}

void solver::PHYSICAL_rho_to_h()
{
//     Nbasis //= physical_object->N_max;
    for(int i1 = 0; i1 < Nbasis; i1++)
    for(int i2 = 0; i2 < Nbasis; i2++)
    {
        h(i1,i2) = physical_object->T(i1,i2);

        for(int j1 = 0; j1 < Nbasis; j1++)
        for(int j2 = 0; j2 < Nbasis; j2++)
            h(i1,i2) += physical_object->V(i1,j2,i2,j1) * rho(j1,j2);
    }

      if (Debug>2) h.print("h:");
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

        // You need to keep the N lowest eigenvectors:
//      D.resize(Nbasis, Nparticles);



//     Particle resizer:

    int Count_proton=0, Count_neutron=0;
    double Isospin;
    for (int i=0; i<Nbasis; i++){
      Isospin=physical_object->Phy_states->state_matrix(i,ST_t);
      (Isospin<0)?Count_proton++:Count_neutron++;
      if (((Isospin<0)&&(Count_proton > Number_of_protons))  ||
          ((Isospin>0)&&(Count_neutron> Number_of_neutrons)))
             D.col(i)=zeros(Nbasis);
    if (Debug>2) cout << physical_object->Phy_states->state_matrix(i,ST_t) << endl;
    if (Debug>2)cout << "protons:  " << Count_proton<< " || " << Number_of_protons<<endl;
    if (Debug>2)cout << "neutrons: " << Count_neutron<< " || " << Number_of_neutrons<<endl;
    }
    if (Debug>2)getchar();


    for (int i=0;i<Nbasis;i++)
      for (int j=0;j<Nbasis;j++)
      if (abs(D(i,j)) < 0.0000001) D(i,j)=0;

      if (Debug>1) e.print("These are the eigenvalues:");
      if (Debug>1) D.print("eigenvectors:");
      if (Debug>1) getchar();


    rho = D*D.t();
    cout.precision(9);
         /**/  cout <<" Total energy: " <<Get_energy() << endl; /**/

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
    for(unsigned int i = 0; (i < max_it); i++)
    {
        std::cout <<endl<<"O=========================================================O"<< endl << "(  "<<i<<"  )"<<endl;

        PHYSICAL_rho_to_h();
        SYSTEM_h_to_rho();

        if((Calc_conv_par() < conv_treshold)&&(i>2))
        {
          cout << "The system converged! "<<endl<<endl<<endl;
            break;
        }

//        alternative_rho();
    }
}

double solver::Get_energy()
{
double Energy = 0;
double Kinetic=0, Direct=0, Excange=0;
int    Nbasis = physical_object->N_max;
        for (int a=0;a<Nbasis;a++)
        for (int b=0;b<Nbasis;b++){
                Kinetic += physical_object->T(a,b)*rho(b,a);
        for (int c=0;c<Nbasis;c++)
        for (int d=0;d<Nbasis;d++){
                Direct += 0.25*physical_object->V(a,b,c,d)*rho(d,b)*rho(c,a);
                Excange -= 0.25*physical_object->V(a,b,c,d)* rho(c,b)*rho(d,a) ;
        }}
        Energy=Direct+Excange+Kinetic;
        if (Debug) cout <<endl<< "Energy components:" <<endl;
        if (Debug) cout << "Kinetic:       " <<Kinetic/2.<<"  MeV"<<endl;
        if (Debug) cout << "Harmonic trap: " <<Kinetic/2.<<"  MeV" <<endl;
        if (Debug) cout << "Potential:     " <<Direct+Excange<<"  MeV"<<endl;
        if (Debug) cout << "  Direct:      "<<Direct<<"  MeV"<<endl;
        if (Debug) cout << "  Excange:     "<<Excange <<"  MeV"<<endl;
        if (Debug) cout << "-------------------------------"<<endl;
        if (Debug) cout << "ENERGY:        "<< Energy<<endl;

return Energy;
}
