#include <armadillo>
#include <cmath>
#include "HF_Phys.h"
#include "States.h"
#include "quadrature.hpp"

#include <iostream>

// this function will fill T and V then make a guess for the initial density
physical_world::physical_world(int N_basis, int N_particles, double (*fill_V)(int, int, int, int, States&, Quad&, double, void*),
                               double (*fill_T)(int, int, States&, void*), void* parameters, States& states, Quad& quad, double b, int sym = 0)
{
    N_part     = N_particles;
    N_max      = N_basis;                           // Set all to zero
    simmetry   = sym;

    if (simmetry==1)
    {
//      if (N_max%2) N_max++;
     T_internal = arma::zeros<mat>(N_max,N_max);
     V_internal = arma::zeros<mat>(N_max*N_max,N_max*N_max);
     for(int i1 = 0; i1<N_max; i1++)                                    //here i take into account spin degeneracy
     {      //cout << i1 << endl;
       T_internal(i1,i1) = fill_T(i1, i1,states, parameters);       //   <---- This will account that: <ab|V|cd> = - <ba|V|cd> = <cd|V|ab>
      for(int i2 = i1; i2<N_max; i2++)
      for(int j1 = 0; j1<N_max; j1++)
      for(int j2 = j1; j2<N_max; j2++){
        if  (((states.state_matrix(i1,3) == states.state_matrix(j1,3))&&(states.state_matrix(i2,3) == states.state_matrix(j2,3))) || ((states.state_matrix(i1,3) == states.state_matrix(j2,3))&&(states.state_matrix(i2,3) == states.state_matrix(j1,3))))
          // check for symmetry conservation in order to avoid overcounting.
        {
         V_internal(i1*N_max+i2, j1*N_max+j2) =
         V_internal(j1*N_max+j2, i1*N_max+i2) =
         V_internal(i2*N_max+i1, j2*N_max+j1) =
         V_internal(j2*N_max+j1, i2*N_max+i1) =
         fill_V(i1, i2, j1, j2, states, quad, b, parameters);           //fill 2b matrix
//                         cout << "  "<< i1<< "  "<<i2<< "  "<<j1<< "  "<<j2 << "   "<<         V_internal(i1*N_max+i2, j1*N_max+j2)  << "    " <<endl;

         V_internal(j2*N_max+j1, i1*N_max+i2) =
         V_internal(j1*N_max+j2, i2*N_max+i1) =
         V_internal(i2*N_max+i1, j1*N_max+j2) =
         V_internal(i1*N_max+i2, j2*N_max+j1) =
         - V_internal(i1*N_max+i2, j1*N_max+j2);
     }}}

    }
    else{
    T_internal = arma::zeros<mat>(N_max,N_max);
    V_internal = arma::zeros<mat>(N_max*N_max,N_max*N_max);
    /*** This is the most general matrix change definition if we know some symmetries****/
    for(int i1 = 0; i1<N_max; i1++)
    for(int i2 = 0; i2<N_max; i2++)
    {
        T_internal(i1,i2) = fill_T(i1, i2,states, parameters);                                  //fill 1b matrix
        for(int j1 = 0; j1<N_max; j1++)
        for(int j2 = 0; j2<N_max; j2++)
        {V_internal(i1*N_max+i2, j1*N_max+j2) = fill_V(i1, i2, j1, j2, states, quad, b, parameters);   //fill 2b matrix
                       cout << "  "<< i1<< "  "<<i2<< "  "<<j1<< "  "<<j2 << "   "<<       fill_V(i1, i2, j1, j2, states, quad, b, parameters) <<  endl;}

    }

    /************************************************************************************/
    }

}


      double physical_world::T(int i1,int j1)                      
      {
//           if (simmetry==1)
//             return (i1 == j1)? T_internal_diag(i1%(N_max/2)):0;
//             else
            return T_internal(i1, j1);
      }
      
      
      double physical_world::V(int i1, int i2, int j1, int j2)
      {return V_internal(i1* N_max + i2, j1 * N_max + j2);}
