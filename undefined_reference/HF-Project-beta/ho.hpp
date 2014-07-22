
#ifndef HO_HPP
#define HO_HPP

#include <iostream>
#include <armadillo>
#include <complex>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_coupling.h>
#include <cmath>
#include <cassert>
#define I_UNIT std::complex<double>(0.,1.)

namespace HO {
	/** forward declarations **/
	arma::Col< std::complex<double>> wfn(const unsigned n, const unsigned l, const unsigned j, const unsigned m, const double r,const double theta,const double phi,const double b);
	double CG(const unsigned l,const int m_l,const unsigned s, const int m_s,const unsigned j, const int m_j);
	std::complex<double> Ylm(const unsigned l, const int ml, const double theta, const double phi);
	double wfn_radial(const unsigned n, const unsigned l,const double r, const double b);

	/** pass 2*j and 2*m as j and m!!! **/
	arma::Col< std::complex<double>> wfn(const unsigned n, const unsigned l, const unsigned j, const unsigned m, const double r,const double theta,const double phi,const double b){
		arma::Col<std::complex<double>> D(2);
		const int ml_up   = m-1;         // is twice ml of upper component
		const int ml_down = m+1;         // is twice ml of lower component
		D(0)     = CG(2*l,ml_up  ,1,1 ,j,m);  // everything is two times its actual value (l,ml,s,ms,j,mj)
		D(1)     = CG(2*l,ml_down,1,-1,j,m); // idem
		std::cout << " CG's : " << D(0) << " , " << D(1) << std::endl;
		if (D(0).real() != 0.) { // only calculate if CG was non-zero!
			D(0) *= Ylm(l,ml_up/2,theta,phi);
		}
		if (D(1).real() != 0.) { // only calculate if CG was non-zero!
			D(1) *= Ylm(l,ml_down/2,theta,phi) ;
		}
		return wfn_radial(n,l,r,b)*D;
	}

	/** pass everything as two times its actual value!!!**/
	double CG(const unsigned l,const int m_l,const unsigned s, const int m_s,const unsigned j, const int m_j) {
		return pow(-1,(-l+s-m_j )/2)*sqrt(j+1)*gsl_sf_coupling_3j(l,s,j,m_l,m_s,-m_j); /**< convert wigner 3j to CG coeff **/
	}
	
	std::complex<double> Ylm(const unsigned l, const int ml, const double theta, const double phi){
		if ( ml>= 0){
			return gsl_sf_legendre_sphPlm(l,ml,cos(theta))*exp(ml*phi*I_UNIT);
		} else {
			return pow(-1,ml)*std::conj(Ylm(l,-ml,theta,phi));
		}
	}
	// b = sqrt( m \omega / hbar )
	double wfn_radial(const unsigned n, const unsigned l,const double r, const double b){
		const double ksi  = r*b;
		const double A_nl = sqrt ( b*b*b / sqrt(M_PI) * exp( (n+l+2)*log(2) + gsl_sf_lnfact(n) - gsl_sf_lndoublefact(2*n+2*l+1) ) );
		return A_nl*exp(-0.5*ksi*ksi)*pow(ksi,l)*gsl_sf_laguerre_n(n,l+0.5,ksi*ksi);
	}
}; // NAMESPACE HO
#endif // HO_HPP
