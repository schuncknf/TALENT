#include "quadrature.hpp"
#include "ho.hpp"
#include "States.h"

#include "Flags.h"
#include "Hamiltonians.cpp"



class states_proprierty
{
    private:
    int A,B,C,D;
    States* Codification;
    public:


    int nA, nB, nC, nD;
    int lA, lB, lC, lD;
    int mA, mB, mC, mD;
    int jA, jB, jC, jD;
    int spinA, spinB, spinC, spinD;
    int isospinA, isospinB, isospinC, isospinD;
    double HO_b;
      
    states_proprierty(int tA,int tB,int tC,int tD, double B_HO, States Cod)
    {
      A=tA ;B=tB ;C=tC ;D=tD;
      nA=Cod.state_matrix(A,ST_n);              nB=Cod.state_matrix(B,ST_n);      nC=Cod.state_matrix(C,ST_n);      nD=Cod.state_matrix(D,ST_n);
      lA=Cod.state_matrix(A,ST_l);              lB=Cod.state_matrix(B,ST_l);      lC=Cod.state_matrix(C,ST_l);      lD=Cod.state_matrix(D,ST_l);
      mA=Cod.state_matrix(A,ST_m);              mB=Cod.state_matrix(B,ST_m);      mC=Cod.state_matrix(C,ST_m);      mD=Cod.state_matrix(D,ST_m);
      jA=Cod.state_matrix(A,ST_j);              jB=Cod.state_matrix(B,ST_j);      jC=Cod.state_matrix(C,ST_j);      jD=Cod.state_matrix(D,ST_j);
      spinA=Cod.state_matrix(A,ST_s);           spinB=Cod.state_matrix(B,ST_s);   spinC=Cod.state_matrix(C,ST_s);   spinD=Cod.state_matrix(D,ST_s);
      isospinA=Cod.state_matrix(A,ST_t);        isospinB=Cod.state_matrix(B,ST_t);isospinC=Cod.state_matrix(C,ST_t);isospinD=Cod.state_matrix(D,ST_t);
      HO_b = B_HO;
//       cout << "eccoB!   "<< HO_b<<endl;
    }
};



// double V_coulomb(int i,int j,int k,int l, int N_basis, double* parm) {return 0;}
// double matrN_drops_VixElement(int,int,int,int, void*);
double fIntegrand(double,double,void*);


double fIntegrand(double r1,double r2,void* p){
        states_proprierty* P = (states_proprierty*)p;
        
        
        //     Minnesota potential for neutron drops //
        int p1 = P->nA;
        int p2 = P->nB;
        int p3 = P->nC;
        int p4 = P->nD;
        
        int l1 = P->lA;
        int l2 = P->lB;
        int l3 = P->lC;
        int l4 = P->lD; 
        
        double b = P->HO_b;
        double v=Potential( r1,  r2, p);

        return     HO::wfn_radial(p1,l1,r1,b)*HO::wfn_radial(p2,l2,r2,b)* 
                   r1*r2*v*  // v = 1/2(VR+VS)
                (  HO::wfn_radial(p3,l3,r1,b)*HO::wfn_radial(p4,l4,r2,b)  // direct part
                  +HO::wfn_radial(p3,l3,r2,b)*HO::wfn_radial(p4,l4,r1,b)); // exchange part, note the reversed radial coordinates
//                   return 0;
}


double General_potential_V(int a,int b,int c,int d, States& states, Quad& I_METHOD, double B, void* param){
        states_proprierty Proprerty(a,b,c,d,B,states);
        void* p = (void*) &Proprerty;
        
        int sigmaA = (int)(2.*states.state_matrix(a,ST_s));
        int sigmaB = (int)(2.*states.state_matrix(b,ST_s));
        int sigmaC = (int)(2.*states.state_matrix(c,ST_s));
        int sigmaD = (int)(2.*states.state_matrix(d,ST_s));        


       if ( sigmaA==sigmaC && sigmaB==sigmaD){
                if (sigmaA==sigmaD)
                        return 0.;
                else
                        return I_METHOD.integrate(fIntegrand,(void*) p);
        } else if ( sigmaA==sigmaD && sigmaB==sigmaC){
                return -I_METHOD.integrate(fIntegrand,(void*) p);
        } else {
                return 0.;
        }
}





// 
