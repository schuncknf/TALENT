group_C: Anton Repko, Deni Vale, Noemi Rocco, Angela Mecca
======

Makefile: contains commands for compiling the programs
- command "make" will compile all

mtrx_test.c: diagonalization of a random matrix with LAPACK routine dsyevr
- compilation on Mac needs addition of "-llapack" in makefile

sho.c: implementation of spherical harmonic oscillator wavefunctions

sho_print.c: interface program for making a table of R_nl(r) values
- use it like "./sho_print 1 20 10 0.01 20. >harm_osc_n=20_l=10.txt"

sho_test.c: testing of orthogonality of R_nl(r) with Simpson integration
- use it like "./sho_test_gauss 1 10 2 10 2"

sho_test_gauss.c: testing of orthogonality of R_nl(r) with Gauss-Laguerre integration

gauleg.c: generation of Gauss-Laguerre abscissas and weights (for given exponent and scaling)

hydrogen.c: diagonalization of hydrogen atom in spherical HO basis
- compilation on Mac needs addition of "-llapack" in makefile

eigen.c: an attempt to separate matrix functions (in progress)

gauher.c: Gauss-Hermite quadrature (in progress)
