#include "quadrature.hpp"
#include "ho.hpp"
#include "States.h"

double T(int i,int j,void* parm) {return (i==j)? (2.*i+1.5) :0; }
double V_coulomb(int i,int j,int k,int l, int N_basis, double* parm) {return 0;}
double Random_rho(int i,int j, int N_basis, double* parm) {return 0;} 

struct fV_neutrondrop_params {
	int n1,n2,n3,n4;
	void *eparams;
};

double fV_neutrondrop(double r1,double r2, void* param){ 
	double mu=1.0;
	double b =1.0; // double check that expression below is correct!
	struct fV_neutrondrop_params p = * (struct fV_neutrondrop_params *) param;
	return  HO::wfn_radial(p.n1,0,0,b)*HO::wfn_radial(p.n2,0,0,b)
		*r1*r2*0.5*exp( - mu*(r1*r1+r2*r2) )*sinh(2.*r1*r2)
		*HO::wfn_radial(p.n3,0,0,b)*HO::wfn_radial(p.n4,0,0,b);
}

double V_neutrondrop(int a,int b,int c, int d,States& states,Quad& q, void* param){
	/** i'th row of state_matrix contains QN's of particle i, 0th column is n **/
	struct fV_neutrondrop_params p;
	p.n1 = states.state_matrix(a,0);
	p.n2 = states.state_matrix(b,0);
	p.n3 = states.state_matrix(c,0);
	p.n4 = states.state_matrix(d,0);
	p.eparams = param;	
	return q.integrate( fV_neutrondrop,(void*) &p);
}
	
	
