#include <math.h>

// Anton Repko, 16th July 2014
// spherical harmonic oscillator - calculation of R_nl(r) in
// psi_nlm(r) = R_nl(r) * Y_lm(theta,phi)

// corresponding formula (in fact, a special recursion relations are used instead):
// R_nl(r) = 2 * (mw/pi)^0.25 * exp(-mw * r^2 / 2) * r^l
// * sqrt(2^l * mw^(l+1) * (2n+2l+1)!! / (2^n * n!))
// * sum(i=0..n){binom(i,n) * (-2*mw*r^2)^i / (2l+2i+1)!!}
// where mw = m * omega / hbar

// calculation through a recursion relation among functions
// of constant (n+l), starting with n = 0
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
  for (i = 3; i <= 2 * npl + 1; i += 2)
    factor *= 2 * mwr2 / i;
  cur = 0.;
  last2 = sqrt(factor);
  if (l == npl)
    cur = last2;
  last = ((npl + 0.5) / x - x) * last2;
  if (l == npl - 1)
    cur = last;
  for (i = npl - 2; i >= l; i--) {
    cur = (((i + 1.5) / x - x) * last - sqrt(npl - i - 1) * last2) / sqrt(npl - i);
    last2 = last;
    last = cur;
  }
  return 2 * pow(mw / M_PI, 0.25) * exp(-0.5*mwr2) * cur;
}
