#include <armadillo>
#include "States.h"

using namespace std;
using namespace arma;

States::States(int N_basis, int sys)     // System:  0: n; 1: n,l,m; 2: n,l,m,s; 3: n,l,m,s,t, 4: n,s; 7: n,l,j,m_j
        {
        system = sys;
        int n,l,m,s,t,N=0,J=0, N_max, counter;
        double j;
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

        case 4:
          N=0;
          n=0;
          J=0;
          while(N<N_basis)
          {
              state_matrix(N,0) = n;
              state_matrix(N,3) = -0.5;
              N++;
              state_matrix(N,0) = n;
              state_matrix(N,3) = +0.5;
              N++;
              n++;
          }
          break;
          
        case 5:
            for (int i=0; i < N_basis;i++)
            {
              state_matrix(i,0) = i%N_basis-N_basis*(2*i/N_basis)/2;
              state_matrix(i,3) = 2*i/N_basis-.5;
            }
         break;

        case 7:

            N = 0;
            n = 0;
            l = 0;
            j = 0;
            m = 0;

            for(counter = 0; N < N_basis; counter++){
                for(l = counter; l > -1; l = l - 2){

                    n = (N - l)/2;

                    for(j = l - 0.5; j <  l + 1.; j++){
                        for(m = -j; m < j + 1.; m++){

                            state_matrix(N,0) = n;
                            state_matrix(N,1) = l;
                            state_matrix(N,2) = j;
                            state_matrix(N,3) = m;

                            N++;

                            if(N > N_basis - 1){
                                break;
                            }
                        }

                        if(N > N_basis - 1){
                            break;
                        }
                    }

                    if(N > N_basis - 1){
                        break;
                    }
                }
            }

            break;

         }
        cout <<endl<<endl<<"mapping:   "<<endl<< state_matrix << endl;
    }
