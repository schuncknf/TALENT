#include <iostream>
#include <armadillo>
#include <cmath>
// #include <string>
#include "Flags.h"

string in_state_file;
string in_potential_file;

using namespace std;
using namespace arma;
#include "quadrature.hpp"
#include "HF_solver.h"
#include "Physics.cpp"
#include "States.h"
#include "HF_Phys.h"
#include "HF_solver.h"

int Debug; 
double (*Potential)(double ,double ,void* );
double (*Kinetics)(int ,int ,States& ,void* );
int contarighe(string the_file);



int main(int argc, char *argv[])
{
  if (argc>3)
  {cout << "press a button to start ("<<argv[1]<<" , "<<argv[2]<<" , "<<argv[3]<<") : "<<endl; getchar();}
  else
  {cout << "ERROR  ::  You need to specify input files (I.HF, spM, VM-scheme) "<<endl;return(1);}
    //////////////////////////SYSTEM//////////////////////////////////////
    // 0  -> 1fermione L=0  NOspin                                      //
    // 1  -> 1fermione L!=0 NOspin                                      //
    // 2  -> 1fermione L!=0 SIspin                                      //
    // 3  -> 2fermioni L!=0 SIspin                                      //
    // 4  -> 1fermione L=0  SIspin                                      //
    // 5  -> 1fermione L=0  SIspin (parametrizzazionea blocchi)         //
    // 7  -> 2fermioni L!=0 SIspin (parametrizzazionea blocchi)         //
    // -1 -> leggi da file l'ordine degli stati                         //
    //////////////////////////////////////////////////////////////////////

    //////////////////////////SIMMETRY////////////////////////////////////
    // This how the potential matrix is filled:                         //
    // 1  -> <ab|V|cd> = - <ba|V|cd> = <cd|V|ab>                        //
    // 2  -> from file                                                  //
    //////////////////////////////////////////////////////////////////////
  

   int    N_basis, N_particles, FILE;
   double hw;
   int system, simmetry;  
   unsigned int        max_iterations;
   double              conv_treshold;
   double              Total_energy, Energy_Protons, Energy_Neutrons;
   int N_n = 2;
   int N_p = 2;
   
   string simulation;
   ifstream inputfile;
   inputfile.open (argv[1]);
   inputfile >> simulation >> N_basis >> N_n>>N_p >> Debug >> hw;
   inputfile.close();
   N_particles=N_n+N_p;



      if      (simulation == "File")            {             FILE=1; simmetry=2;}
      if      (simulation == "L=0")             {system = 0;  FILE=0; simmetry=1;}





/********************************SOME INPUT***************************************/
    max_iterations = 5600;
    conv_treshold = 1e-15;
    Potential = Minnesota;
    Kinetics  = T_full_HO;

    Quad glqi = Quad("LG_64.int");   // global integrator used only if not reading from file



    if (FILE){
    in_state_file     = argv[2];
    in_potential_file = argv[3];
    N_basis = contarighe(in_state_file)-1;
    }
/********************************END INPUT***************************************/





//      double b = 0.49104389631014383; // remember b = sqrt(m omega /hbar) hbar omega = 10 MeV, m=938.9059
     double b = sqrt(.5/20.73 * hw);


    cout << "basis:                      " << N_basis << endl;
    cout << "Particles:                  " << N_particles << endl;
    cout << "  Protons:                  " << N_p << endl;
    cout << "  Neutrons:                 " << N_n << endl;
    cout << "b:                          " << b << endl;
    cout << "hw:                         " << hw << endl;
    cout << "Max_iterations:             " << max_iterations << endl;
    cout << "treshold:                   " << conv_treshold << endl;


    Total_energy = Energy_Neutrons = Energy_Protons= 0;


    //******************Sistema di neutroni********************************//
    if (N_n)
    {
    cout << "    COMPUTING NEUTRONS    "<<endl;
    if (FILE) system = 1;
    States              stati_Neutrons(N_basis, system);
    physical_world      Neutrons(N_basis, N_n, General_potential_V,Kinetics, &hw, stati_Neutrons, glqi, b, simmetry);
    solver              HF_solver_Neutrons(&Neutrons, max_iterations, conv_treshold, N_n, N_p);
    Energy_Neutrons = HF_solver_Neutrons.Get_energy();
    Total_energy += Energy_Neutrons;
    }


    //******************Sistema di protoni********************************//
    if (N_p && (FILE))
    {
    cout << "    COMPUTING PROTONS    "<<endl;
    system = 2;
    States              stati_Protons(N_basis, system);
    physical_world      Protons(N_basis, N_p, General_potential_V,Kinetics, &hw, stati_Protons, glqi, b, simmetry);
    solver              HF_solver_Protons(&Protons, max_iterations, conv_treshold, N_n, N_p);
    Energy_Protons = HF_solver_Protons.Get_energy();
    Total_energy += Energy_Protons;
    }


       //***********************Output*************************************************//
      /**/  cout.precision(9) ;                                                     /**/
     /**/  cout << endl<<endl<<" Total energy: " <<Total_energy << endl;           /**/
    /**/  cout  <<"    Proton  : " <<Energy_Protons << endl;                      /**/
   /**/  cout  <<"    Neutrons: " <<Energy_Neutrons << endl;                     /**/
  //______________________________________________________________________________//

    return 0;}








int contarighe(string the_file){
    ifstream file;
    file.open (the_file);
    string s;
    int contarighe=0;
        while(!file.eof()){
                           getline(file,s,'\n');
                           contarighe++;
        }
        if (Debug>1) cout << "the file has "<<contarighe-1<< " rows"<<endl;
        return contarighe-1; // In quanto conta anche l'ultima riga vuota
}