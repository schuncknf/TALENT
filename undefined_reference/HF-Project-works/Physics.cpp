#include "quadrature.hpp"
#include "ho.hpp"
#include "States.h"

// double V_coulomb(int i,int j,int k,int l, int N_basis, double* parm) {return 0;}
double matrN_drops_VixElement(int,int,int,int, void*);
double fIntegrand(double,double,void*);

double T_full_HO(int i,int j,States& states,void* parm) {
    double *wb = (double*)parm;
    return (i==j)? *wb*((2.*states.state_matrix(i,0))+1.5) : 0.;
}




double fIntegrand(double r1,double r2,void* p){
        double*n = (double *) p;
        //struct SPQN p1,p2,p3,p4;
        int p1 = (int)n[0]; //*spqns[n[0]];// p1.print();
        int p2 = (int)n[1]; //spqns[n[1]];// p2.print();
        int p3 = (int)n[2]; //spqns[n[2]];// p3.print();
        int p4 = (int)n[3]; //spqns[n[3]];// p4.print();
        double kappaR=1.487;
        double kappaS=0.465;
        double b = n[4]; //0.49104389631014383; // remember b = sqrt(m omega /hbar) hbar omega = 10 MeV, m=938.9059 ,hbarc = 197.32891
        double V0R=200.0;
        double V0S=91.85;
        double VR =  0.25*V0R/kappaR*( exp( -kappaR*pow(r1-r2,2) ) - exp( -kappaR*pow(r1+r2,2)) );
        double VS = -0.25*V0S/kappaS*( exp( -kappaS*pow(r1-r2,2) ) - exp( -kappaS*pow(r1+r2,2)) );
        double v = 0.5*(VR+VS);
        return     HO::wfn_radial(p1,0,r1,b)*HO::wfn_radial(p2,0,r2,b)* 
                   r1*r2*v*  // v = 1/2(VR+VS)
                (  HO::wfn_radial(p3,0,r1,b)*HO::wfn_radial(p4,0,r2,b)  // direct part
                  +HO::wfn_radial(p3,0,r2,b)*HO::wfn_radial(p4,0,r1,b)); // exchange part, note the reversed radial coordinates
}

int delta(int a,int b){
        return a==b;
}


double N_drops_V(int a,int b,int c,int d, States& states, Quad& I_METHOD, double B, void* param){
        double n[] = {states.state_matrix(a,0),states.state_matrix(b,0),states.state_matrix(c,0),states.state_matrix(d,0), B};
        //spqns[a].print(); spqns[b].print(); spqns[c].print(); spqns[d].print();
        /*int dd = delta(spqns[a].sigma,spqns[c].sigma)*delta(spqns[b].sigma,spqns[d].sigma)-delta(spqns[a].sigma,spqns[d].sigma)*delta(spqns[b].sigma,spqns[c].sigma);
        if (dd!=0){
                return dd*I_METHOD.integrate(fIntegrand,(void*)n);          
        }*/
        int sigmaA = (int)(2.*states.state_matrix(a,3));
        int sigmaB = (int)(2.*states.state_matrix(b,3));
        int sigmaC = (int)(2.*states.state_matrix(c,3));
        int sigmaD = (int)(2.*states.state_matrix(d,3));
        
//         cout << a <<" "<<sigmaA <<  " "<<n[0]<<" || " << b<<" "<<sigmaB<<  " "<<n[1]<<endl;
        
        if ( sigmaA==sigmaC && sigmaB==sigmaD){
                if (sigmaA==sigmaD)
                        return 0.;
                else
                        return I_METHOD.integrate(fIntegrand,(void*) n);
        } else if ( sigmaA==sigmaD && sigmaB==sigmaC){
                return -I_METHOD.integrate(fIntegrand,(void*) n);
        } else {
                return 0.;
        }
}
