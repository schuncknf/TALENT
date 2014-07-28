#include <math.h>
#include <gsl/gsl_sf_laguerre.h>

// Anton Repko, 28th July 2014
// spherical harmonic oscillator - calculation of R_nl(r) in
// psi_nlm(r) = R_nl(r) * Y_lm(theta,phi)

// corresponding formula (in fact, a function from GSL is used):
// R_nl(r) = 2 * (mw/pi)^0.25 * exp(-mw * r^2 / 2) * r^l
// * sqrt(2^l * mw^(l+1) * (2n+2l+1)!! / (2^n * n!))
// * sum(i=0..n){binom(i,n) * (-2*mw*r^2)^i / (2l+2i+1)!!}
// where mw = m * omega / hbar

double sho_wf(double r, double mw, int n, int l)
{
  int i, npl;
  double x, mwr2, cur, last, last2, factor;
  if ((r < 0) || (mw < 0) || (n < 0) || (l < 0))
    return 0.;
  if (r == 0) {
    if (l != 0)
      return 0.;
    factor = mw;
    for (i = 1; i <= n; i++)
      factor *= (i + 0.5) / i;
    return 2 * sqrt(sqrt(mw / M_PI) * factor);
  }
  x = sqrt(mw) * r;
  mwr2 = x * x;
  npl = n + l;
  factor = mw;
  for (i = 1; i <= n; i++)
    factor *= 2 * i / (double)(2 * i + 1);
  for (i = 2 * n + 3; i <= 2 * npl + 1; i += 2)
    factor *= 2 * mwr2 / i;
  return 2 * pow(mw / M_PI, 0.25) * sqrt(factor) * exp(-0.5*mwr2)
         * gsl_sf_laguerre_n(n, l + 0.5, mwr2);
}
