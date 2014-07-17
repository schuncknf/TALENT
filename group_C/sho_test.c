#include <stdio.h>
#include <stdlib.h>
#include "sho.h"
#define RSTEP 0.01
#define RMAX 50.

// checks the orthogonality of spherical harmonic oscillator radial wavefunctions
// by Simpson integration of R_(n1 l1)(r) and R_(n2 l2). The oscillator strength
// is given by mw = m*omega/hbar
int main(int argc, char *argv[])
{
    int n1, l1, n2, l2;
    double mw, r, r2, sum;
    if (argc != 6) {
        printf("Usage: ./sho_test mw n1 l1 n2 l2\n");
        return 1;
    }
    mw = atof(argv[1]);
    n1 = atoi(argv[2]);
    l1 = atoi(argv[3]);
    n2 = atoi(argv[4]);
    l2 = atoi(argv[5]);
    sum = 0.;
    for (r = 0.; r <= RMAX; r += RSTEP) {
        sum += sho_wf(r, mw, n1, l1) * sho_wf(r, mw, n2, l2) * r * r;
        r2 = r + 0.5*RSTEP;
        sum += 0.5 * sho_wf(r2, mw, n1, l1) * sho_wf(r2, mw, n2, l2) * r2 * r2;
    }
    printf("Overlap is %.10lf\n", sum * RSTEP / 1.5);
    return 0;
}
