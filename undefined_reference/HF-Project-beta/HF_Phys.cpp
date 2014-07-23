#include <armadillo>
#include <cmath>
#include "HF_Phys.h"

// this function will fill T and V then make a guess for the initial density
physical_world::physical_world(int N_basis, int N_particles, double (*fill_V)(int, int, int, int, int, void*),
                               double (*fill_T)(int, int, int, void*),
                               double (*guess_rho)(int, int, int, void*),
                               void* parameters);
{
    N_part = N_particles;
    N_max      = N_basis;                           // Set all to zero
    T_internal = zeros<mat>(N_max,N_max);
    V_internal = zeros<mat>(N_max,N_max);
    rho        = zeros<mat>(N_max,N_max);

    /*** This is the most general matrix change definition if we know some symmetries****/
    for(int i1; i1<N_max; i1++)
    for(int i2; i2<N_max; i2++)
    {
        T_internal(i1,i2) = fill_T(i1, i2, N_max, parameters);                                  //fill 1b matrix
        rho(i1,i2) = guess_rho(i1, i2, N_max, parameters);                                      //fill a starting rho

        for(int j1; j1<N_max; j1++)
        for(int j2; j2<N_max; j2++)
            V_internal(i1*N_max+i2, j1*N_max+j2) = fill_V(i1, i2, j1, j2, N_max, parameters);   //fill 2b matrix
    }
    /************************************************************************************/
}
