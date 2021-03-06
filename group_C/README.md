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

hydrogen.c: diagonalization of hydrogen atom in spherical HO basis (now modular design)
- compilation on Mac needs addition of "-llapack" in makefile

eigen.c: an attempt to separate matrix functions

gauher.c: Gauss-Hermite quadrature (to be used for the on-fly generation of weights)

sho_test_gauher.c: testing of orthogonality of R_nl(r) with Gauss-Hermite integration - much better than Gauss-Laguerre. It uses fixed number of abscissas: 64, scaling is set accoring to hw. This is good up to polynomials of degree 124, i.e. 2n+l < 125.

hydrogen-bench.c: program for creation of the benchmark file
