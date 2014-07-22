


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
    vec e_old;
    vec e;

    // This is the interface that the main will see in order to use V and T
    double V(int i1, int i2, int j1, int j2)     {return V_internal(i1*N_max+i2,j1+N_max*j2);}
    double T(int i1,int j1)                      {return T_internal(i1, j1);}

    // this function will fill T and V then make a guess for the initial density
    Physical_world(int N_basis, double (*fill_V)(int,int,int,int,int,double*), double (*fill_T)(int,int,int,double*), double (*guess_rho)(int,int,int,double*), double* parameters)
        {
        N_max      = N_basis;                           // Set all to zero
        T_internal = zeros<mat>(N_max,N_max);
        V_internal = zeros<mat>(N_max,N_max);
        rho        = zeros<mat>(N_max,N_max);

        /*** This is the most general matrix change definition if we know some symmetries****/
        for(int i1; i1<N_max; i1++)
        for(int i2; i2<N_max; i2++)
          {
          T_internal(i1,i2)   = fill_T(i1,i2,N_max,parameters);                                      //fill 1b matrix
          rho(i1,i2) = guess_rho(i1,i2,N_max,parameters);                                   //fill a starting rho
          for(int j1; j1<N_max; j1++)
          for(int j2; j2<N_max; j2++)
            V_internal(i1*N_max+i2,j1*N_max+j2) = fill_V(i1,i2,j1,j2, N_max,parameters);             //fill 2b matrix
          }
        /************************************************************************************/
    }
};






// ***************************************************************************************//
void PHYSICAL_rho_to_h(Physical_world* hf)
{

        for(int i1; i1<hf->N_max; i1++)
        for(int i2; i2<hf->N_max; i2++)
        {
          hf->h(i1,i2) = hf->T(i1,i2);
          for(int j1; j1<hf->N_max; j1++)
          for(int j2; j2<hf->N_max; j2++)
            hf->h(i1,i2) += hf->V(i1,j2,i2,j1)  *  hf->rho(j1,j2);
        }
}



void SYSTEM_h_to_rho(Physical_world* hf)
{
    //if (hf->complex_flag)
    //else
mat D;

        hf->e_old=hf->e;
        hf->rho_old  =  hf->rho;
        eig_sym(hf->e,D,hf->h);

        mat Dstar=conj(D);

        for(int i1; i1<hf->N_max; i1++)
        for(int i2; i2<hf->N_max; i2++)
        for(int i3; i3<hf->N_max; i3++)
        hf->rho(i1,i2)     =  D(i1,i3)*Dstar(i2,i3);
}




    bool convergence_check(Physical_world* hf)    {return(   abs(hf->e(0)-hf->e_old(0))  <   1.e-6)?1:0;  }
// ***************************************************************************************//