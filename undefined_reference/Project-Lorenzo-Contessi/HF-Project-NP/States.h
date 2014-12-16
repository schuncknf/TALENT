#ifndef STATES_H
#define STATES_H

#include <armadillo>

using namespace std;
using namespace arma;

class States
{
    private:
    public:
      mat state_matrix;
      int system;
      
      // neeeded for custom potential
      int N_max_custom_but_real;
      vec custom_interface;
      
    States(int N_basis, int system);     // System:  0: n; 1: n,l,m; 2: n,l,m,s; 3: n,l,m,s,t
};

#endif
