#include <armadillo>
#include <cmath>
#include "HF_Phys.h"
#include "States.h"
#include "quadrature.hpp"

#include <iostream>

#include "Flags.h"


// this function will fill T and V then make a guess for the initial density
physical_world::physical_world(int N_basis, int N_particles, double (*fill_V)(int, int, int, int, States&, Quad&, double, void*),
                               double (*fill_T)(int, int, States&, void*), void* parameters, States& states, Quad& quad, double b, int sym = 0)
{
    double       Value;
    N_part     = N_particles;
    N_max      = N_basis;                           // Set all to zero
    simmetry   = sym;
    Phy_states = &states;


    switch (simmetry) {

      case 1:
      {

          ofstream Potential_out;
          Potential_out.open ("V_matrix.dat");
          T_internal = arma::zeros<mat>(N_max,N_max);
          V_internal = arma::zeros<mat>(N_max*N_max,N_max*N_max);
          for(int i1 = 0; i1<N_max; i1++)                                       //here i take into account spin degeneracy
          {
            T_internal(i1,i1) = fill_T(i1, i1,states, parameters);  //   <---- This will account that: <ab|V|cd> = - <ba|V|cd> = <cd|V|ab>
            for(int i2 = i1; i2<N_max; i2++)
            for(int j1 = 0 ; j1<N_max; j1++)
            for(int j2 = j1; j2<N_max; j2++)
            if   ( ((states.state_matrix(i1,ST_s) == states.state_matrix(j1,ST_s))&&(states.state_matrix(i2,ST_s) == states.state_matrix(j2,ST_s)))
              || ( ( states.state_matrix(i1,ST_s) == states.state_matrix(j2,ST_s))&&(states.state_matrix(i2,ST_s) == states.state_matrix(j1,ST_s))))
              {
              Value =          fill_V(i1, i2, j1, j2, states, quad, b, parameters);           //fill 2b matrix
              if (Debug && (Value!=0)) cout << "Coefficient:  "<< i1<< "  "<<i2<< "  "<<j1<< "  "<<j2 << "   "<<         Value  <<"    =     V_internal ( "<<i1*N_max+i2<< ", "<< j1*N_max+j2<<" )  =  "<<Value<<endl;
              V_internal(i1*N_max+i2, j1*N_max+j2) =
              V_internal(j1*N_max+j2, i1*N_max+i2) =
              V_internal(i2*N_max+i1, j2*N_max+j1) =
              V_internal(j2*N_max+j1, i2*N_max+i1) = Value;
              V_internal(j2*N_max+j1, i1*N_max+i2) =
              V_internal(j1*N_max+j2, i2*N_max+i1) =
              V_internal(i2*N_max+i1, j1*N_max+j2) =
              V_internal(i1*N_max+i2, j2*N_max+j1) = -Value;
          }} if (Debug) {cout<< ">> ended to calculate coefficients..." <<endl ; getchar(); }
          Potential_out.close();
          break;}







      case 2:       /*            Neutrons          */
      {


        if (Debug) cout <<">> Reading potential files:  "<<endl; 
          N_max = states.N_max_custom_but_real;                   //states read the correct number of states!
          T_internal = arma::zeros<mat>(N_max,N_max);
          V_internal = arma::zeros<mat>(N_max*N_max,N_max*N_max);

         int n1,n2,m1,m2;
         char useless[100];
         double DOG;
         double* hw=(double*)parameters;



             for (int j=0;j<N_max;j++)
              T_internal(j,j) = fill_T(j,j,states,parameters);

          //////////////////////////////////////////////////////////////////////////////
          //     construction of the V matrix                                         //
          //////////////////////////////////////////////////////////////////////////////
         int count=0;
         ifstream V_file;
         V_file.open (in_potential_file);
         if (V_file.is_open()) {
         V_file >>useless>>useless>>useless>>useless>>useless>>useless>>useless>>useless>>useless>>useless>>useless>>useless>>
                  useless>>useless>>useless>>useless>>useless>>useless>>useless>>useless>>useless>>useless;
         do{
           V_file >> n1>>n2>>m1>>m2>>DOG;
//            states.custom_interface.print();
//            cout << n1 << " " << n1 << " " << m1 << " " << m2 << " " <<endl; 
           if((states.custom_interface(n1)>=0)&&(states.custom_interface(n2)>=0)&&(states.custom_interface(m1)>=0)&&(states.custom_interface(m2)>=0)) {
             V_internal(states.custom_interface(n1)*N_max+states.custom_interface(n2),states.custom_interface(m1)*N_max+states.custom_interface(m2)) = (abs(DOG)>1e-10)?DOG:0;
             if ((Debug>1) && (DOG!=0)) {
               cout <<count++<<":    morten labels ( "<<n1<<"  "<<n2<<"  "<<m1<<"  "<<m2<<")  My corrispondences: "<<states.custom_interface(n1)<<" "<<states.custom_interface(n2)<<" "<<states.custom_interface(m1)<<" "<<states.custom_interface(m2)<<")    =     V_internal( "<<states.custom_interface(n1)*N_max+states.custom_interface(n2)<<" , "<<states.custom_interface(m1)*N_max+states.custom_interface(m2)<<" )  =  "<< DOG << endl;
             }
           }}while (!V_file.eof()) ;}
         V_file.close();
         if (Debug) cout << ">> File coefficients loaded sucesfully! "<<endl;

         break;}





    }
    /************************************************************************************/


}




      double physical_world::T(int i1,int j1)
      {return T_internal(i1, j1);}

      double physical_world::V(int i1, int i2, int j1, int j2)
      {return V_internal(i1* N_max + i2, j1 * N_max + j2);}
