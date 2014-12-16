#include <cmath>


double T_full_HO(int i,int j,States& states,void* parm) {
    double *wb = (double*)parm;
    int na = states.state_matrix(i,ST_n);
    int nb = states.state_matrix(j,ST_n);
    int la = states.state_matrix(i,ST_l);
    int lb = states.state_matrix(j,ST_l);
    int ma = states.state_matrix(i,ST_m);
    int mb = states.state_matrix(j,ST_m);
    int sa = states.state_matrix(i,ST_s);
    int sb = states.state_matrix(j,ST_s);
    int ta = states.state_matrix(i,ST_t);
    int tb = states.state_matrix(j,ST_t);
    int ja = states.state_matrix(i,ST_j);
    int jb = states.state_matrix(j,ST_j);

    if ((la!=lb)||(ma!=mb)||(sa!=sb)||(ta!=tb)||(ja!=jb)) return 0;

    switch (na-nb){
       case  1:  return *wb*sqrt(na*(na+la+0.5));  break;
       case  0:  return *wb*((2.*na+la)+1.5);      break;
       case -1:  return *wb*sqrt(nb*(nb+lb+0.5));  break;
      default:  return 0;                          break; }
}

double T_S_HO(int i,int j,States& states,void* parm) {
    double *wb = (double*)parm;
    return (i==j)? *wb*((2.*(states.state_matrix(i,ST_n))+states.state_matrix(i,ST_l))+1.5) : 0.;
}

double Minnesota(double r1,double r2,void* p)
{
   double kappaR=1.487;
        double kappaS=0.465;
        double V0R=200.0;
        double V0S=91.85;
        double VR =  0.25*V0R/kappaR*( exp( -kappaR*pow(r1-r2,2) ) - exp( -kappaR*pow(r1+r2,2)) );
        double VS = -0.25*V0S/kappaS*( exp( -kappaS*pow(r1-r2,2) ) - exp( -kappaS*pow(r1+r2,2)) );
        double v  = 0.5*(VR+VS);
  
        return v;
}


double Minnesota_coulomb(double r1,double r2,void* p)
{
   double kappaR=1.487;
        double kappaS=0.465;
        double V0R=200.0;
        double V0S=91.85;
        double VR =  0.25*V0R/kappaR*( exp( -kappaR*pow(r1-r2,2) ) - exp( -kappaR*pow(r1+r2,2)) );
        double VS = -0.25*V0S/kappaS*( exp( -kappaS*pow(r1-r2,2) ) - exp( -kappaS*pow(r1+r2,2)) );
        double v  = 0.5*(VR+VS);
        
        v+=1.44/abs(r1-r2);
        
  
        return v;
}