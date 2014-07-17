// #include "Fill_HO_Matrix.cpp"
// GSL - Include
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_eigen.h>              // gsl eigensystem routines
#include <gsl/gsl_integration.h>        // gsl integration routines
#include <boost/concept_check.hpp>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_errno.h>



//*************************INTEGRATING ROUTINE ******************
class integral_parameters
{
    private:
  double lower_limit;            /* start integral from 1 (to infinity) */
  double abs_error;              /* to avoid round-off problems */
  double rel_error;              /* the result will usually be much better */
  double *result;                /* the result from the integration */
  double *error;                 /* the estimated error from the integration */
    public:
  void *param;

       integral_parameters(double Llower_limit, double Labs_error, double Lrel_error,double *Lresult, double *Lerror)
       {
        int i;
        lower_limit = Llower_limit;
        abs_error   = Labs_error;
        rel_error   = Lrel_error;
        result      = Lresult;
        error       = Lerror;
       }
  void integrate(double (*function)(double, void*))
       {
        gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (1000);
        gsl_function function_2b_integrated;
        function_2b_integrated.function     = function;
        function_2b_integrated.params       = param;
        gsl_integration_qagiu (&function_2b_integrated,lower_limit,abs_error, rel_error, 1000, work_ptr, result, error);

       }

};





// Parameters of HO
typedef struct                  // structure holding Hij parameters 
{
  int i;                        // 1st matrix index 
  int j;                        // 2nd matrix index 
  int l1;
  int l2;
  double mass;                  // particle mass 
  double b_ho;                  // harmonic oscillator parameter 
  int potential_index;          // indicates which potential (now there is only coulomb)
  double hbar;
}
HO_parameters;


double Hij (HO_parameters ho_parameters);
double Hij_integrand (double x, void *params_ptr);
double ho_radial (int n, int l, double b, double r);
double norm (int n, int l, double b);
double ho_epsilon (int n, int l, double b, double m);
//kinetic
double T_ij(double hb, int i, int li, int j, int lj);
// potentials
double V_coulomb (double r, double C1);



double Hij (HO_parameters ho_parameters)
{
  double Parameters[8] = {(double)ho_parameters.i,(double)ho_parameters.j,ho_parameters.mass,ho_parameters.b_ho,ho_parameters.hbar,ho_parameters.potential_index,ho_parameters.l1,ho_parameters.l2};
  double Result, Error;
  double  omega  = ho_parameters.hbar / (ho_parameters.mass * ho_parameters.b_ho * ho_parameters.b_ho);   //Definition of b
  integral_parameters Hij_int(0, 1.e-8,1.e-8,&Result,&Error);
  Hij_int.param = (void*)Parameters;
  Hij_int.integrate(&Hij_integrand);
#ifdef Teasy
   Result += T_ij(ho_parameters.hbar*omega,ho_parameters.i,ho_parameters.l1,ho_parameters.j,ho_parameters.l2);
#endif
  return (Result);              // send back the result of the integration 
}

double Hij_integrand (double r, void *param)
{
  double *parameters = (double*)param;
  int potential_index;                          // index 1,2,... for potental  (nonusato)
  int l1 = (int)parameters[6];                  // orbital angular momentum 
  int l2 = (int)parameters[7];                  // orbital angular momentum 
  int Li, Lj;                                   // principal quantum number (1,2,...) 
  double mass, b_ho;                            // local ho parameters 
  double hbar;
  double omega;                                 // harmonic oscillator frequency 
  double ho_pot;                                // value of ho potential at current r 

  Li      = (int)parameters[0]+1;
  Lj      = (int)parameters[1]+1;
  mass    = parameters[2];
  b_ho    = parameters[3];
  hbar    = parameters[4];
  omega  = hbar / (mass * b_ho * b_ho);                              // definition of omega 
  ho_pot = (1. / 2.) * mass * (omega * omega) * (r * r);             // ho pot'l 


  double C1=1.;              // Coulomb factor z e^2
  double Hij_r;
  Hij_r =  V_coulomb(r, C1);
#ifndef Teasy
  Hij_r += ho_epsilon (Lj, l2, b_ho, mass) - ho_pot;
  //Smart way to calculate the kinetic energy!
#endif
  return (ho_radial (Li, l1, b_ho, r) * Hij_r * ho_radial (Lj, l2, b_ho, r));
}

double ho_epsilon (int n, int l, double b, double m)
{  return ((2. * ((double) n - 1.) + (double) l + 3. / 2.) / (m * b * b)); }
















double
ho_radial (int n, int l, double b, double r)
{
  ////////////////////////////////////////////////////////////////////
  // b= sqrt(h/(wm) )                                               //
  ////////////////////////////////////////////////////////////////////
  double q = r / b;
  double qq = q * q;
  
  ////////////////////////////////////////////////////////////////////
  // gsl_pow_int(q, (l + 1)) <- you have not to multiply again by r //
  ////////////////////////////////////////////////////////////////////
  return    norm(n, l, b) * gsl_pow_int(q, (l + 1)) * exp (-qq / 2.) * gsl_sf_laguerre_n((n - 1), l + 1. / 2., qq);
}

double norm (int n, int l, double b)
{   return sqrt (2. * gsl_sf_fact ((unsigned) (n - 1)) / (b * gsl_sf_gamma ((double) n + (double) l + 1. / 2.))); }


double V_coulomb (double r, double C1)
{  return (-C1 / r); }

double T_ij(double hw, int i, int li, int j, int lj)
{
  if (li==lj)
  {
//     i++;   this label start from 0
//     j++;
    switch (i-j)
    {
    case (-1): 
      return (.5*hw*sqrt(j*(j+lj+.5)));
        break;
    case (0): 
      return (.5*hw*(2*i+1.5));
        break;
    case (+1): 
      return (.5*hw*sqrt(i*(i+li+.5)));
        break;
    default: 
      return (0.);;
        break;
    }
  }
  else
    return (0.);
}