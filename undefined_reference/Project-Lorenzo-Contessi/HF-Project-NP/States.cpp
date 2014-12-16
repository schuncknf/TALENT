#include <armadillo>
#include "States.h"

#include "Flags.h"

using namespace std;
using namespace arma;



States::States(int N_basis, int sys)     // System:  0: n; 1: n,l,m; 2: n,l,m,s; 3: n,l,m,s,t, 4: n,s: -1 construct the matrix by hands
        {
        system = sys;
        int n,l,m,s,t,N=0,J=0, N_max, counter;
	double j;
        state_matrix = zeros<mat>(N_basis,6);
         int n_maxy=1;
         int l_maxy=1;
         int /*N,*/N_custom,/*n,l,m,*/j2,mj2,t2,spin2;
         char useless[100];
         double DOG;
         ifstream T_file;
        
        switch(system)
        {

         /*************************************************/
         /*     CASE default:  2 fermion L=0                    */
         /*************************************************/
        default:
          N=0;
          n=0;
          J=0;
          cout << "         N  n  L  m  s      t" << endl;
          while(N<N_basis)
          {
              state_matrix(N,ST_n) = n;
              state_matrix(N,ST_s) = -0.5;
              cout << "States:  " << N <<"  "<< n<<"  "<<0<<"  "<<0<<"  "<<state_matrix(N,ST_s)<<"  "<<state_matrix(N,ST_t)<< endl;
              N++;
              state_matrix(N,ST_n) = n;
              state_matrix(N,ST_s) = +0.5;
              cout << "States:  " << N <<"  "<< n<<"  "<<0<<"  "<<0<<"  "<<state_matrix(N,ST_s)<<"  "<<state_matrix(N,ST_t)<< endl;
              N++;
              n++;
          }
          getchar(); cout<<endl;
          break;



          //////////////////////////////////////////////////////////////////////////////
          //     read the Single particle enegy file and construct the map of states  //
          //////////////////////////////////////////////////////////////////////////////

        case 1:                 //        NEUTRONS
         if (Debug>1) cout << "Mapping from a file... "<<endl;
         custom_interface = zeros<vec>(N_basis+1);  //# of state can't be more than the # of convention file rows
         T_file.open(in_state_file);
         if (T_file.is_open()) {
          T_file >>useless>>useless>> useless>>useless>>useless>>useless;
          N=0;
          N_custom=0;
          do{
            if (N_custom<99) 
              T_file >>useless>>useless>> N_custom>> n>>l>>j2>>mj2>>t2;
            else{
              T_file >>useless>>useless>> n>>l>>j2>>mj2>>t2; N_custom++;}
            if (t2==1) //is a neutron
            {
              if (Debug>1)cout << std::showpos<<"Mapping:   "<<N_custom<<"  "<<N<<"  "<<n<<" "<<l<<" "<<j2<<" "<<mj2<<" "<<t2<<" "<<endl;
              custom_interface(N_custom) = N;
              state_matrix(N,ST_n) = (double)n;
              state_matrix(N,ST_l) = (double)l;
              state_matrix(N,ST_j) = (double)j2;
              state_matrix(N,ST_m) = (double)mj2;
              state_matrix(N,ST_t) = (double)t2/2.;
              N++;
            }
            else
            {
              custom_interface(N_custom) = -1;
              if (Debug>1) cout << std::showpos<<"Mapping:   "<<N_custom<<"  "<<N<<"  "<<n<<" "<<l<<" "<<j2<<" "<<mj2<<" "<<t2<<" REJECT"<<endl;
            }
          }while (!T_file.eof() && (N_custom<N_basis)) ;}
         T_file.close();
         N_max_custom_but_real = N;    //This is the real number of states;
         if (Debug>1)cout << ">> End mappature Neutrons...  "<<endl;
         getchar();
         break;



        case 2:                 //        PROTONS
         if (Debug)cout << "Mapping from a file... "<<endl;
         custom_interface = zeros<vec>(N_basis+1);  //# of state can't be more than the # of convention file rows
         T_file.open(in_state_file);
         if (T_file.is_open()) {
          T_file >>useless>>useless>> useless>>useless>>useless>>useless;
          N=0;
          N_custom=0;
          do{
            if (N_custom<99) 
              T_file >>useless>>useless>> N_custom>> n>>l>>j2>>mj2>>t2;
            else{
              T_file >>useless>>useless>> n>>l>>j2>>mj2>>t2; N_custom++;}
            if (t2==-1) //is a proton
            {
              if (Debug) cout << std::showpos<<"Mapping:   "<<N_custom<<"  "<<N<<"  "<<n<<" "<<l<<" "<<j2<<" "<<mj2<<" "<<t2<<" "<<endl;
              custom_interface(N_custom) = N;
              state_matrix(N,ST_n) = (double)n;
              state_matrix(N,ST_l) = (double)l;
              state_matrix(N,ST_j) = (double)j2;
              state_matrix(N,ST_m) = (double)mj2;
              state_matrix(N,ST_t) = (double)t2/2.;
              N++;
            }
            else
            {
              custom_interface(N_custom) = -1;
              if (Debug) cout << std::showpos<<"Mapping:   "<<N_custom<<"  "<<N<<"  "<<n<<" "<<l<<" "<<j2<<" "<<mj2<<" "<<t2<<" REJECT"<<endl;
            }
          }while (!T_file.eof() && (N_custom<N_basis)) ;}
         T_file.close();
         N_max_custom_but_real = N;    //This is the real number of states;
         if (Debug) cout << ">> End mappature Protons...  "<<endl;
         getchar();
         break;

        }}
