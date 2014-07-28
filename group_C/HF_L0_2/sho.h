#ifndef _sho_h
#define _soh_h
// spherical harmonic oscillator - calculation of R_nl(r) in
// psi_nlm(r) = R_nl(r) * Y_lm(theta,phi)
// where mw = m * omega / hbar
double sho_wf(double r, double mw, int n, int l);
#endif
