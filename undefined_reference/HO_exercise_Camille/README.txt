
Hydrogen harmonic oscillator expansion exercise
================================================

just doing
$> cd src/ && make
should do the trick.
gsl, blas, armadillo are required.
Executables will be moved to run/

The compiler will spit out a load of warnings but 
these are all from some ugly header i found somewhere
floating on the interwebs to calculate the Gauss Laguerre
Quadrature knots and weights.
The knots and weights are stored in binary format.
