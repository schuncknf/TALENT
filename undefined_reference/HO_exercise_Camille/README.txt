
Hydrogen harmonic oscillator expansion exercise
================================================

just doing
$> cd src/ && make
should do the trick.
gsl, blas, armadillo are required.
Executables will be moved to run/

quad_generator
===============
will generate knots and weights files 
needed by the other executables

test
==============
tests of several integrands and 
orthonormality, uses knots and weights
created by quad_generator

main
==============
the main program of the exercise
calculates the ground state energy
of the ydrogen atom in natural units.
Also uses knots and weights created
by quad_generator




The compiler will spit out a load of warnings but 
these are mostly from some ugly header i found somewhere
floating on the interwebs to calculate the Gauss Laguerre
Quadrature knots and weights.
The knots and weights are stored in binary format.
