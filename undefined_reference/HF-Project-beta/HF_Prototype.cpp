#include <iostream>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;


class Physical_world 
{
    private:
      // Those are the matrix stored (i dont what the main to see those, they can be stored in different way)
      mat V_internal;
      mat T_internal;
      bool complex_flag; // 1 for complex stuff

    public:
    int N_max;     //Base dimension

    // Store the actual density and the hamiltonian elements
    mat rho_old;                // needed for better convergence?
    mat rho;
    mat h;


    // This is the interface that the main will see in order to use V and T
    double V(int i1, int i2, int j1, int j2)     {return *V_internal(i1*N_max+i2,j1+N_max*j2);}
    double T(int i1,int j1)                      {return *T_internal(i1, j1);}

    // this function will fill T and V then make a guess for the initial density
    void Initialize(int N_basis, double (*fill_V)(int,int,int,int,int,void*), double (*fill_T)(int,int,int,void*), double (*guess_rho)(void*), void* parameters)
        {
        N_max      = N_basis;                           // Set all to zero
        T_internal = zeros<mat>(N_max,N_max);
        V_internal = zeros<mat>(N_max,N_max);
        rho        = zeros<mat>(N_max,N_max);

        /*** This is the most general matrix change definition if we know some symmetries****/
        for(int i1; i1<N_max; i1++)
        for(int i2; i2<N_max; i2++)
          {
          T(i1,i2)   = fill_T(i1,i2,N_max,parameters);                                      //fill 1b matrix
          rho(i1,i2) = guess_rho(i1,i2,N_max,parameters);                                   //fill a starting rho
          for(int j1; j1<N_max; j1++)
          for(int j2; j2<N_max; j2++)
            V(i1*N_max+i2,j1*N_max+j2) = fill_V(i1,i2,j1,j2, N_max,parameters);             //fill 2b matrix
          }
        /************************************************************************************/
    }

    bool convergence_check()    {return(rho == rho_old)?1;0;}
    
};






int main()
{

Physical_world HF;
HF.Initialize(100, T_cinetik, V_coulomb, Random_rho, {1,2,3,4});

while (HF.convergence_check()) 
{
PHYSICAL_rho_to_h(HF);
SYSTEM_h_to_rho(HF);

}

return 0;
}


// ***************************************************************************************//
void PHYSICAL_rho_to_h(Physical_world* hf)
{
    //if (hf->complex_flag)
    //else
hf->h = ...
  
}

void SYSTEM_h_to_rho(Physical_world* hf)
{
    //if (hf->complex_flag)
    //else
mat D;
vec e;
  
}

// ***************************************************************************************//






double T_cinetik(int i,int j, int N_basis, void* parm) {return (i==j)?1:0;}
double V_coulomb(int i,int j,int k,int l, int N_basis, void* parm) {return 0;}
double Random_rho(int i,int j, int N_basis, void* parm) {return 0;}