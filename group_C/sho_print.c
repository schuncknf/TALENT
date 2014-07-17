#include <stdio.h>
#include <stdlib.h>
#include "sho.h"

// prints the values of spherical harmonic oscillator wavefunctions R_nl(r)
// where mw = m*omega/hbar, n and l are quantum numbers
int main(int argc, char *argv[])
{
    int n, l;
    double mw, r, rstep, rmax;
    if (argc != 6) {
        printf("Usage: ./sho_print mw n l rstep rmax\n");
        return 1;
    }
    mw = atof(argv[1]);
    n = atoi(argv[2]);
    l = atoi(argv[3]);
    rstep = atof(argv[4]);
    rmax = atof(argv[5]);
    for (r = 0.; r <= rmax; r += rstep)
        printf("%.5lf\t%.8le\n", r, sho_wf(r, mw, n, l));
    return 0;
}
