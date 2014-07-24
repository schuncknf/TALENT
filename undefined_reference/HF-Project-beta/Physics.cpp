#include "quadrature.hpp"
#include "ho.hpp"
#include "States.h"

double T(int i,int j,void* parm) {return (i==j)? ((i - i%2)+1.5) :0; }
double V_coulomb(int i,int j,int k,int l, int N_basis, double* parm) {return 0;}
double Random_rho(int i,int j, int N_basis, double* parm) {return 0;} 

struct fV_neutrondrop_params {
    int n1, n2, n3, n4;
    int s1, s2, s3, s4;
    double b;
	void *eparams;
};

double fV_neutrondrop(double r1,double r2, void* param){ 
	double mu=1.0;
	struct fV_neutrondrop_params p = * (struct fV_neutrondrop_params *) param;

    double kappaR=1.487;
    double kappaS=0.465;
    double V0R=200.0;
    double V0S=91.85;
    double VR =  0.25*V0R/kappaR*( exp( -kappaR*pow(r1-r2,2) ) - exp( -kappaR*pow(r1+r2,2)) );
    double VS = -0.25*V0S/kappaS*( exp( -kappaS*pow(r1-r2,2) ) - exp( -kappaS*pow(r1+r2,2)) );
    double v = 0.5*(VR+VS);
    return HO::wfn_radial(p.n1,0,r1,p.b)*HO::wfn_radial(p.n2,0,r2,p.b)*
           r1*r2*v*  // v = 1/2(VR+VS)
        (  HO::wfn_radial(p.n3,0,r1,p.b)*HO::wfn_radial(p.n4,0,r2,p.b)  // direct part
          +HO::wfn_radial(p.n3,0,r2,p.b)*HO::wfn_radial(p.n4,0,r1,p.b)); // exchange part, note the reversed radial coordinates
}

double V_neutrondrop(int a,int b,int c, int d, States& states, Quad& q, double par_b_value, void* param){
	/** i'th row of state_matrix contains QN's of particle i, 0th column is n **/
	struct fV_neutrondrop_params p;
	p.n1 = states.state_matrix(a,0);
	p.n2 = states.state_matrix(b,0);
	p.n3 = states.state_matrix(c,0);
	p.n4 = states.state_matrix(d,0);
    p.s1 = 2 * states.state_matrix(a,1);
    p.s2 = 2 * states.state_matrix(b,1);
    p.s3 = 2 * states.state_matrix(c,1);
    p.s4 = 2 * states.state_matrix(d,1);
    p.b = par_b_value;
	p.eparams = param;	


/***********************************************
 * To check the kinetic energy matrix elements:

    return 0.;

 ***********************************************/

    if ( p.s1==p.s3 && p.s2==p.s4){
        if(p.s1 == p.s2)
        {
            return 0.;
        }
        else
        {
            return q.integrate(fV_neutrondrop,(void*) &p);
        }
    } else if ( p.s1==p.s4 && p.s2==p.s3){
        return -q.integrate(fV_neutrondrop,(void*) &p);
    } else {
        return 0.;
    }
}
