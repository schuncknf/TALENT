#include <armadillo>
#include "States.h"

using namespace std;
using namespace arma;

States::States(int N_basis, int system)     // System:  0: n; 1: n,l,m; 2: n,l,m,s; 3: n,l,m,s,t
        {
        int n,l,m,s,t,N=0,J=0, j, N_max;
        state_matrix = zeros<mat>(N_basis,5);
        switch(system)
        {
          case 0: 
            for (int i=0; i < N_basis;i++)
              state_matrix(i,0) = i;
            break;
          case 1:
            N=0;
            while(N<N_basis)
            {
            for (n=0;((n<=N/2)&&(J<N_basis));n++)
              {l=N-2*n;
              for(m=-l;((m<=l)&&(J<N_basis));m++)
                {
                state_matrix(J,0) = n;
                state_matrix(J,1) = l;
                state_matrix(J,2) = m;
                J++;
                }
              }
              N++;
            }
          break;
          
          case 2:
            J=N=0;
            N_max= N_basis/2;
            for(s=0;s<=1;s++)
            {
              while(N<N_max)
              {
              for (n=0;((n<=N/2)&&(J<N_max));n++)
                {l=N-2*n;
                for(m=-l;((m<=l)&&(J<N_max));m++)
                  {
                  state_matrix(J+s*N_basis/2,0) = n;
                  state_matrix(J+s*N_basis/2,1) = l;
                  state_matrix(J+s*N_basis/2,2) = m;
                  state_matrix(J+s*N_basis/2,3) = s-.5;
                  J++;
//                    cout << J+s*N_basis/2<< "     " << n<<l<<m<<"  "<<s<<"  "<<endl;
                  }
                }
                N++;
              }
              J=N=0;
              N_max = N_basis-N_max;
             }
          break;
          
          case 3:
            J=N=0;
            N_max= N_basis/4;
            for(t=0;t<=1;t++)
            {
              for(s=0;s<=1;s++)
              {
                while(N<N_max)
                {
                for (n=0;((n<=N/2)&&(J<N_max));n++)
                  {l=N-2*n;
                  for(m=-l;((m<=l)&&(J<N_max));m++)
                    {
                    state_matrix(J+t*N_basis/2+s*N_basis/4,0) = n;
                    state_matrix(J+t*N_basis/2+s*N_basis/4,1) = l;
                    state_matrix(J+t*N_basis/2+s*N_basis/4,2) = m;
                    state_matrix(J+t*N_basis/2+s*N_basis/4,3) = s-.5;
                    state_matrix(J+t*N_basis/2+s*N_basis/4,4) = t-.5;
                    J++;
//                   cout << J+t*N_basis/2+s*N_basis/4<< "     " << n<<l<<m<<"  "<<s<<"  "<<t<<endl;
                    }
                  }
                  N++;
                }
                J=N=0;
                if (t) N_max = N_basis/4+N_basis%4-1;
              }

             }
            
          break;
          }
        }
